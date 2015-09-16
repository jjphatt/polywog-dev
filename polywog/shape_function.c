// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "polywog/shape_function.h"

struct shape_function_kernel_t
{
  char* name;
  void* context;
  void (*compute)(void* context, point_t* points, real_t* extents, int num_points, point_t* x, real_t* value, vector_t* gradient);
  void (*dtor)(void* context);
};

static void shape_function_kernel_free(void* ctx, void* dummy)
{
  shape_function_kernel_t* kernel = ctx;
  string_free(kernel->name);
  if ((kernel->dtor != NULL) && (kernel->context != NULL))
    kernel->dtor(kernel->context);
}

shape_function_kernel_t* shape_function_kernel_new(const char* name,
                                                   void* context,
                                                   void (*compute)(void* context, point_t* points, real_t* extents, int num_points, point_t* x, real_t* values, vector_t* gradients),
                                                   void (*dtor)(void* context))
{
  shape_function_kernel_t* kernel = GC_MALLOC(sizeof(shape_function_kernel_t));
  kernel->name = string_dup(name);
  kernel->context = context;
  kernel->compute = compute;
  kernel->dtor = dtor;
  GC_register_finalizer(kernel, shape_function_kernel_free, kernel, NULL, NULL);
  return kernel;
}

void shape_function_kernel_compute(shape_function_kernel_t* kernel, 
                                   point_t* points, 
                                   real_t* extents, 
                                   int num_points, 
                                   point_t* x, 
                                   real_t* values,
                                   vector_t* gradients)
{
  kernel->compute(kernel->context, points, extents, num_points, x, 
                  values, gradients);
}

static void simple_compute(void* context, 
                           point_t* points, 
                           real_t* extents, 
                           int num_points, 
                           point_t* x, 
                           real_t* values, 
                           vector_t* gradients)
{
  real_t eta_max = *((real_t*)context);
  for (int i = 0; i < num_points; ++i)
  {
    point_t* xi = &points[i];
    real_t h = extents[i];
    real_t Di = point_distance(xi, x);
    values[i] = MAX(0.0, 2.0*eta_max - Di*Di/(h*h));
    real_t dWdDi = -2.0*Di/(h*h);
    gradients[i].x = dWdDi * (xi->x-x->x)/Di;
    gradients[i].y = dWdDi * (xi->y-x->y)/Di;
    gradients[i].z = dWdDi * (xi->z-x->z)/Di;
  }
}

static void spline4_compute(void* context, 
                            point_t* points, 
                            real_t* extents, 
                            int num_points, 
                            point_t* x, 
                            real_t* values, 
                            vector_t* gradients)
{
  real_t eta_max = *((real_t*)context);
  for (int i = 0; i < num_points; ++i)
  {
    point_t* xi = &points[i];
    real_t h = extents[i];
    real_t Di = point_distance(xi, x);
    real_t eta = Di/h;
    real_t X = eta/eta_max;
    if (X < 1.0)
    {
      values[i] = 1.0 - 6.0*X*X + 8.0*X*X*X - 3.0*X*X*X*X;
      real_t dWdDi = (-12.0*X + 24.0*X*X - 12.0*X*X*X)/(eta_max*h);
      gradients[i].x = dWdDi * (xi->x-x->x)/Di;
      gradients[i].y = dWdDi * (xi->y-x->y)/Di;
      gradients[i].z = dWdDi * (xi->z-x->z)/Di;
    }
    else
    {
      values[i] = 0.0;
      gradients[i].x = 0.0;
      gradients[i].y = 0.0;
      gradients[i].z = 0.0;
    }
  }
}

shape_function_kernel_t* simple_shape_function_kernel_new(real_t eta_max)
{
  char name[1024];
  snprintf(name, 1023, "Simple kernel: W(x, x0, h) = %g - (||x-x0||/h)**2", 2.0*eta_max);
  real_t* context = polymec_malloc(sizeof(real_t));
  *context = eta_max;
  return shape_function_kernel_new((char*)name, context, simple_compute, polymec_free);
}

shape_function_kernel_t* spline4_shape_function_kernel_new(real_t eta_max)
{
  char name[1024];
  snprintf(name, 1023, "Spline4 kernel");
  real_t* context = polymec_malloc(sizeof(real_t));
  *context = eta_max;
  return shape_function_kernel_new((char*)name, context, spline4_compute, polymec_free);
}

struct shape_function_t 
{
  char* name;
  void* context;
  shape_function_vtable vtable;
  int i, N;
};

shape_function_t* shape_function_new(const char* name, 
                                     void* context, 
                                     shape_function_vtable vtable)
{
  ASSERT(vtable.neighborhood_size != NULL);
  ASSERT(vtable.get_neighborhood_points != NULL);
  ASSERT(vtable.compute != NULL);

  shape_function_t* phi = polymec_malloc(sizeof(shape_function_t));
  phi->name = string_dup(name);
  phi->context = context;
  phi->vtable = vtable;
  phi->i = -1;
  phi->N = -1;
  return phi;
}

void shape_function_free(shape_function_t* phi)
{
  string_free(phi->name);
  if ((phi->vtable.dtor != NULL) && (phi->context != NULL))
    phi->vtable.dtor(phi->context);
  polymec_free(phi);
}

void shape_function_set_neighborhood(shape_function_t* phi, int point_index)
{
  ASSERT(point_index >= 0);
  phi->i = point_index;
  phi->N = phi->vtable.neighborhood_size(phi->context, point_index);
  if (phi->vtable.set_neighborhood != NULL)
    phi->vtable.set_neighborhood(phi->context, phi->i);
}

int shape_function_num_points(shape_function_t* phi)
{
  return phi->N;
}

void shape_function_get_points(shape_function_t* phi, point_t* points)
{
  phi->vtable.get_neighborhood_points(phi->context, phi->i, points);
}

void shape_function_compute(shape_function_t* phi, 
                            point_t* x,
                            real_t* values,
                            vector_t* gradients)
{
  phi->vtable.compute(phi->context, phi->i, x, values, gradients);
}

