// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/options.h"
#include "core/dense_local_matrix.h"
#include "core/linear_algebra.h"
#include "core/special_functions.h"
//#include "core/silo_file.h"
#include "polywog/gmls_matrix.h"
#include "polywog/mlpg_quadrature.h"
#include "make_mlpg_lattice.h"
#include "poisson_gmls_functional.h"
#include "elastic_gmls_functional.h"

#undef I // no imaginary stuff here!

void test_gmls_matrix_ctor(void** state)
{
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  make_mlpg_lattice(&bbox, 10, 10, 10, 3.0, &points, &extents, &stencil);
  multicomp_poly_basis_t* P = standard_multicomp_poly_basis_new(1, 2);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(P, W, points, extents, stencil);

  // Clean up.
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}

// Franke's function is a solution to Poisson's equation.
static void franke(void* context, point_t* x, real_t* u)
{
  real_t X = x->x, Y = x->y;
  real_t term1 = 0.75 * exp(-0.25 * (pow(9.0*X - 2.0, 2.0) + pow(9.0*Y - 2.0, 2.0)));
  real_t term2 = 0.75 * exp(-1.0/49.0 * pow(9.0*X + 1.0, 2.0) - 1.0/10.0 * pow(9.0*Y + 1.0, 2.0));
  real_t term3 = 0.50 * exp(-0.25 * (pow(9.0*X - 7.0, 2.0) + pow(9.0*Y - 3.0, 2.0)));
  real_t term4 = 0.20 * exp(-pow(9.0*X - 4.0, 2.0) - pow(9.0*Y - 7.0, 2.0));
  u[0] = term1 + term2 + term3 + term4;
}

// Here's its gradient.
static void grad_franke(void* context, point_t* x, real_t* grad_u)
{
  real_t X = x->x, Y = x->y;
  real_t h1_X = hermite_hn(1, X);
  real_t h1_Y = hermite_hn(1, Y);
  real_t term1 = 0.75 * exp(-0.25 * (pow(9.0*X - 2.0, 2.0) + pow(9.0*Y - 2.0, 2.0)));
  real_t ddx_term1 = 0.25 * h1_X * term1 * 9.0;
  real_t ddy_term1 = 0.25 * h1_Y * term1 * 9.0;
  real_t term2 = 0.75 * exp(-1.0/49.0 * pow(9.0*X + 1.0, 2.0) - 1.0/10.0 * pow(9.0*Y + 1.0, 2.0));
  real_t ddx_term2 = 1.0/49.0 * h1_X * term2 * 9.0;
  real_t ddy_term2 = 1.0/49.0 * h1_Y * term2 * 9.0;
  real_t term3 = 0.50 * exp(-0.25 * (pow(9.0*X - 7.0, 2.0) + pow(9.0*Y - 3.0, 2.0)));
  real_t ddx_term3 = 0.25 * h1_X * term3 * 9.0;
  real_t ddy_term3 = 0.25 * h1_Y * term3 * 9.0;
  real_t term4 = 0.20 * exp(-pow(9.0*X - 4.0, 2.0) - pow(9.0*Y - 7.0, 2.0));
  real_t ddx_term4 = h1_X * term4 * 9.0;
  real_t ddy_term4 = h1_Y * term4 * 9.0;
  grad_u[0] = ddx_term1 + ddx_term2 + ddx_term3 + ddx_term4;
  grad_u[1] = ddy_term1 + ddy_term2 + ddy_term3 + ddy_term4;
  grad_u[2] = 0.0;
}

void test_gmls_matrix_with_frankes_function(void** state)
{
  // Construct a point cloud, set of extents, and a stencil.
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  int nx = 10, ny = 10;
  bool neumann = false; // No Neumann BCs by default.
  
  // Override options if desired.
  {
    options_t* opts = options_argv();
    char* res_str = options_value(opts, "resolution");
    if (res_str != NULL)
    {
      if (string_casecmp(res_str, "low") == 0)
        nx = ny = 10;
      else if (string_casecmp(res_str, "medium") == 0)
        nx = ny = 20;
      else if (string_casecmp(res_str, "high") == 0)
        nx = ny = 40;
    }
    char* neumann_str = options_value(opts, "neumann");
    if (neumann_str != NULL)
      neumann = string_as_boolean(neumann_str);
  }

  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0/nx};
  make_mlpg_lattice(&bbox, nx, ny, 1, 3.0, &points, &extents, &stencil);
  log_debug("Point cloud has %d points and %d ghosts.", points->num_points, points->num_ghosts);
  int N = points->num_points + points->num_ghosts;
//  point_cloud_fprintf(points, stdout);

  // Now set up the GMLS machinery.
  multicomp_poly_basis_t* P = standard_multicomp_poly_basis_new(1, 2);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(P, W, points, extents, stencil);
  real_t delta = 0.5; // ratio of subdomain extent to point extent.
  gmls_functional_t* lambda = poisson_gmls_functional_new(2, points, extents, delta);
  sp_func_t* F = sp_func_from_func("Franke's function", franke, SP_FUNC_HETEROGENEOUS, 1);
  sp_func_t* grad_F = sp_func_from_func("Franke's function gradient", grad_franke, SP_FUNC_HETEROGENEOUS, 3);
  volume_integral_t* Qv = mlpg_cube_volume_integral_new(points, extents, 2, delta);

  // Set up our linear system using a dense matrix. This is inefficient but very simple.
  local_matrix_t* A = dense_local_matrix_new(N);
  real_t B[N]; // RHS vector.

  // Treat boundary nodes first.
  int_unordered_set_t* boundary_nodes = int_unordered_set_new();
  if (neumann) // Neumann BCs on -x/+x boundaries, Dirichlet on others.
  {
    tagger_unite_tag(points->tags, "neumann_boundary", "-x");
    tagger_unite_tag(points->tags, "neumann_boundary", "+x");
    tagger_unite_tag(points->tags, "dirichlet_boundary", "-y");
    tagger_unite_tag(points->tags, "dirichlet_boundary", "+y");
    tagger_unite_tag(points->tags, "dirichlet_boundary", "-z");
    tagger_unite_tag(points->tags, "dirichlet_boundary", "+z");
    int num_nbnodes, num_dbnodes; 
    int* nbnodes = point_cloud_tag(points, "neumann_boundary", &num_nbnodes);
    ASSERT(nbnodes != NULL);
    log_debug("Found %d Neumann boundary nodes in point cloud.", num_nbnodes);
    int* dbnodes = point_cloud_tag(points, "dirichlet_boundary", &num_dbnodes);
    ASSERT(dbnodes != NULL);
    log_debug("Found %d Dirichlet boundary nodes in point cloud.", num_dbnodes);

    // Boundary functionals.
    gmls_functional_t* dirichlet_bc = gmls_matrix_dirichlet_bc_new(matrix);
    sp_func_t* n = point_cloud_property(points, "normal");
    st_func_t* nt = st_func_from_sp_func(n);
    gmls_functional_t* neumann_bc = gmls_matrix_neumann_bc_new(matrix, nt);

    for (int b = 0; b < num_nbnodes; ++b) // Neumann BC nodes
    {
      int bnode = nbnodes[b];
      int num_coeffs = gmls_matrix_num_coeffs(matrix, bnode);
      int rows[num_coeffs], cols[num_coeffs];
      real_t coeffs[num_coeffs];
      gmls_matrix_compute_coeffs(matrix, bnode, neumann_bc, 0.0, NULL, 
                                 rows, cols, coeffs);

      // Make sure all the coefficients go in the same row.
      for (int j = 1; j < num_coeffs; ++j)
        assert_true(rows[j] == rows[0]);

      // Now dump the coefficients into our matrix.
      real_t row_vector[N];
      memset(row_vector, 0, sizeof(real_t) * N);
      for (int j = 0; j < num_coeffs; ++j)
        row_vector[cols[j]] = coeffs[j];
      local_matrix_add_row_vector(A, 1.0, bnode, row_vector);
      int_unordered_set_insert(boundary_nodes, bnode);

      // RHS contribution.
      point_t* xb = &points->points[bnode];
      real_t dF[3], nb[3];
      sp_func_eval(grad_F, xb, dF);
      sp_func_eval(n, xb, nb); 
      B[bnode] = nb[0]*dF[0] + nb[1]*dF[1] + nb[2]*dF[2];
    }

    for (int b = 0; b < num_dbnodes; ++b) // Dirichlet BC nodes
    {
      int bnode = dbnodes[b];
      int num_coeffs = gmls_matrix_num_coeffs(matrix, bnode);
      int rows[num_coeffs], cols[num_coeffs];
      real_t coeffs[num_coeffs];
      gmls_matrix_compute_coeffs(matrix, bnode, dirichlet_bc, 0.0, NULL, 
                                 rows, cols, coeffs);

      // Make sure all the coefficients go in the same row.
      for (int j = 1; j < num_coeffs; ++j)
        assert_true(rows[j] == rows[0]);

      // Now dump the coefficients into our matrix.
      real_t row_vector[N];
      memset(row_vector, 0, sizeof(real_t) * N);
      for (int j = 0; j < num_coeffs; ++j)
        row_vector[cols[j]] = coeffs[j];
      local_matrix_add_row_vector(A, 1.0, bnode, row_vector);
      int_unordered_set_insert(boundary_nodes, bnode);

      // RHS contribution.
      point_t* xb = &points->points[bnode];
      sp_func_eval(F, xb, &B[bnode]);
    }

    // Clean up the boundary functionals.
    gmls_functional_free(dirichlet_bc);
    gmls_functional_free(neumann_bc);
  }
  else // Dirichlet BCs all around.
  {
    int num_bnodes; 
    int* bnodes = point_cloud_tag(points, "boundary", &num_bnodes);
    ASSERT(bnodes != NULL);
    log_debug("Found %d boundary nodes in point cloud.", num_bnodes);

    // Boundary functional.
    gmls_functional_t* dirichlet_bc = gmls_matrix_dirichlet_bc_new(matrix);

    for (int b = 0; b < num_bnodes; ++b)
    {
      int bnode = bnodes[b];
      int num_coeffs = gmls_matrix_num_coeffs(matrix, bnode);
      int rows[num_coeffs], cols[num_coeffs];
      real_t coeffs[num_coeffs];
      gmls_matrix_compute_coeffs(matrix, bnode, dirichlet_bc, 0.0, NULL, 
                                 rows, cols, coeffs);

      // Make sure all the coefficients go in the same row.
      for (int j = 1; j < num_coeffs; ++j)
        assert_true(rows[j] == rows[0]);

      // Now dump the coefficients into our matrix.
      real_t row_vector[N];
      memset(row_vector, 0, sizeof(real_t) * N);
      for (int j = 0; j < num_coeffs; ++j)
        row_vector[cols[j]] = coeffs[j];
      local_matrix_add_row_vector(A, 1.0, bnode, row_vector);
      int_unordered_set_insert(boundary_nodes, bnode);

      // RHS contribution.
      point_t* xb = &points->points[bnode];
      sp_func_eval(F, xb, &B[bnode]);
    }

    // Clean up the boundary functional.
    gmls_functional_free(dirichlet_bc);
  }

  // Now interior nodes.
  for (int i = 0; i < N; ++i)
  {
    if (!int_unordered_set_contains(boundary_nodes, i))
    {
      int num_coeffs = gmls_matrix_num_coeffs(matrix, i);
      int rows[num_coeffs], cols[num_coeffs];
      real_t coeffs[num_coeffs];
      gmls_matrix_compute_coeffs(matrix, i, lambda, 0.0, NULL, 
                                 rows, cols, coeffs);

      // Make sure all the coefficients go in the same row.
      for (int j = 1; j < num_coeffs; ++j)
        assert_true(rows[j] == rows[0]);

      // Now dump the coefficients into our matrix.
      real_t row_vector[N];
      memset(row_vector, 0, sizeof(real_t) * N);
      for (int j = 0; j < num_coeffs; ++j)
        row_vector[cols[j]] = coeffs[j];
      local_matrix_add_row_vector(A, 1.0, i, row_vector);

      // Interior RHS.
      volume_integral_set_domain(Qv, i);
      volume_integral_compute(Qv, F, &B[i]);
    }
  }

  // Clean up the Poisson functional.
  gmls_functional_free(lambda);
//printf("A = ");
//local_matrix_fprintf(A, stdout);
//printf("B = ");
//vector_fprintf(B, nx*ny, stdout);

  // Solve the linear system.
  real_t U[N];
  bool solved = local_matrix_solve(A, B, U);
  assert_true(solved);

//printf("U = ");
//vector_fprintf(U, nx*ny, stdout);

  // Measure the L2 error.
  real_t L2 = 0.0;
  for (int i = 0; i < points->num_points; ++i)
  {
    point_t* xi = &points->points[i];
    real_t Usol;
    sp_func_eval(F, xi, &Usol);
    real_t err = fabs(U[i] - Usol);
    L2 += err*err;
  }

  // Scale the L2 error by the uniform volume element.
  real_t dV = 1.0/nx * 1.0/ny;
  L2 = sqrt(L2 * dV);
  log_urgent("L2(error) = %g\n", L2);

  // Super crude hard-wired success metric. :-/
  if (nx == 5)
    assert_true(L2 < 0.05);
  else if (nx == 10)
    assert_true(L2 < 0.0165);
  else if (nx == 20)
    assert_true(L2 < 0.0045);

  // Clean up.
  int_unordered_set_free(boundary_nodes);
  local_matrix_free(A);
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}

typedef struct
{
  real_t L, D, P, E, nu;
} cantileaver_t;

// The solution for a 2D cantileaver beam problem.
static void cantileaver_beam(void* context, point_t* x, real_t* u)
{
  cantileaver_t* c = context;
  real_t X = x->x, Y = x->y;
  real_t L = c->L, D = c->D, P = c->P, E = c->E, nu = c->nu;
  real_t I = D*D*D/12.0;
  u[0] = -P/(6.0*E*I) * (Y-0.5*D) * (3.0*X * (2.0*L - X) + (2.0 + nu) * Y * (Y - D));
  u[1] = -P/(6.0*E*I) * (X*X * (3.0*L - X) + 3.0*nu * (L - X) * (Y - 0.5*D) * (Y - 0.5*D) + 0.25*(4.0 + 5.0*nu) *D*D * X);
  u[2] = 0.0;
}

void test_gmls_matrix_with_cantileaver_beam(void** state)
{
  // Construct a point cloud, set of extents, and a stencil.
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  int nx = 33, ny = 5, nz = 5;
  
  // Material properties.
  cantileaver_t* cantileaver = polymec_malloc(sizeof(cantileaver_t));
  cantileaver->L = 8.0;
  cantileaver->D = 1.0;
  cantileaver->P = 1.0;
  cantileaver->E = 1.0;
  cantileaver->nu = 0.25;

  // Override options if desired.
  {
    options_t* opts = options_argv();
    char* res_str = options_value(opts, "resolution");
    if (res_str != NULL)
    {
      if (string_casecmp(res_str, "low") == 0)
      {
        nx = 16;
        ny = nz = 3;
      }
      else if (string_casecmp(res_str, "medium") == 0)
      {
        nx = 33;
        ny = nz = 5;
      }
      else if (string_casecmp(res_str, "high") == 0)
      {
        nx = 66;
        ny = nz = 10;
      }
    }
  }

  real_t L = cantileaver->L;
  bbox_t bbox = {.x1 = 0.0, .x2 = L, .y1 = 0.0, .y2 = L*ny/nx, .z1 = 0.0, .z2 = L*nz/nx};
  make_mlpg_lattice(&bbox, nx, ny, nz, 3.0, &points, &extents, &stencil);
  log_debug("Point cloud has %d points and %d ghosts.", points->num_points, points->num_ghosts);
  int N = points->num_points + points->num_ghosts;
//  point_cloud_fprintf(points, stdout);

  // Now set up the GMLS machinery.
  multicomp_poly_basis_t* P = standard_multicomp_poly_basis_new(3, 2);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(P, W, points, extents, stencil);
  real_t delta = 0.5; // ratio of subdomain extent to point extent.
  gmls_functional_t* lambda = elastic_gmls_functional_new(cantileaver->E, cantileaver->nu,
                                                          2, points, extents, delta);
  sp_func_vtable cantileaver_vtable = {.eval = cantileaver_beam, .dtor = polymec_free};
  sp_func_t* F = sp_func_new("Cantileaver beam solution", cantileaver, 
                             cantileaver_vtable, SP_FUNC_HETEROGENEOUS, 3);
  volume_integral_t* Qv = mlpg_cube_volume_integral_new(points, extents, 2, delta);

  // Set up our linear system using a dense matrix. This is inefficient but very simple.
  local_matrix_t* A = dense_local_matrix_new(3*N);
  real_t B[3*N];

  // Treat boundary nodes first.
  int num_bnodes; 
  int* bnodes = point_cloud_tag(points, "boundary", &num_bnodes);
  ASSERT(bnodes != NULL);
  gmls_functional_t* dirichlet_bc = gmls_matrix_dirichlet_bc_new(matrix);
  log_debug("Found %d boundary nodes in point cloud.", num_bnodes);
  int_unordered_set_t* boundary_nodes = int_unordered_set_new();
  for (int b = 0; b < num_bnodes; ++b)
  {
    int bnode = bnodes[b];
    int num_coeffs = gmls_matrix_num_coeffs(matrix, bnode);
    int rows[num_coeffs], cols[num_coeffs];
    real_t coeffs[num_coeffs];
    gmls_matrix_compute_coeffs(matrix, bnode, dirichlet_bc, 0.0, NULL, 
                               rows, cols, coeffs);

    // Now dump the coefficients into our matrix.
    real_t row_vectors[3][3*N];
    memset((void*)row_vectors, 0, sizeof(real_t) * 3*3*N);
    for (int j = 0; j < num_coeffs; ++j)
    {
//printf("%d: row %d\n", bnode, rows[j] % 3);
      int which = rows[j] % 3;
      row_vectors[which][cols[j]] = coeffs[j];
    }
    local_matrix_add_row_vector(A, 1.0, 3*bnode,   row_vectors[0]);
    local_matrix_add_row_vector(A, 1.0, 3*bnode+1, row_vectors[1]);
    local_matrix_add_row_vector(A, 1.0, 3*bnode+2, row_vectors[2]);

    // RHS contributions.
    point_t* xb = &points->points[bnode];
    sp_func_eval(F, xb, &B[3*bnode]);

    int_unordered_set_insert(boundary_nodes, bnode);
  }

  // Now interior nodes.
  for (int i = 0; i < points->num_points; ++i)
  {
    if (!int_unordered_set_contains(boundary_nodes, i))
    {
      int num_coeffs = gmls_matrix_num_coeffs(matrix, i);
      int rows[num_coeffs], cols[num_coeffs];
      real_t coeffs[num_coeffs];
      gmls_matrix_compute_coeffs(matrix, i, lambda, 0.0, NULL, 
                                 rows, cols, coeffs);

      // Now dump the coefficients into our matrix.
      real_t row_vectors[3][3*N];
      memset((void*)row_vectors, 0, sizeof(real_t) * 3*3*N);
      for (int j = 0; j < num_coeffs; ++j)
      {
//printf("%d: row %d\n", i, rows[j] % 3);
        int which = rows[j] % 3;
        row_vectors[which][cols[j]] = coeffs[j];
      }
      local_matrix_add_row_vector(A, 1.0, 3*i,   row_vectors[0]);
      local_matrix_add_row_vector(A, 1.0, 3*i+1, row_vectors[1]);
      local_matrix_add_row_vector(A, 1.0, 3*i+2, row_vectors[2]);

      // Interior RHS.
      volume_integral_set_domain(Qv, i);
      volume_integral_compute(Qv, F, &B[3*i]);
    }
  }
//  real_t diag_A[3*N];
//  local_matrix_get_diagonal(A, diag_A);
//printf("diag A = [");
//for (int i = 0; i < 3*N; ++i)
//printf("%d: %g\n", i, diag_A[i]);
//printf("]\n");
//vector_fprintf(diag_A, 3*N, stdout);
//printf("A = ");
//local_matrix_fprintf(A, stdout);
//printf("B = ");
//vector_fprintf(B, 3*N, stdout);

  // Clean up the functionals.
  gmls_functional_free(dirichlet_bc);
  gmls_functional_free(lambda);

  // Solve the linear system.
  real_t U[3*N];
  bool solved = local_matrix_solve(A, B, U);
  assert_true(solved);

//printf("U = ");
//vector_fprintf(U, nx*ny, stdout);

  // Measure the L2 error.
  real_t L2 = 0.0;
  for (int i = 0; i < points->num_points; ++i)
  {
    point_t* xi = &points->points[i];
    real_t Usol[3];
    sp_func_eval(F, xi, Usol);
    real_t ux = U[3*i], uy = U[3*i+1], uz = U[3*i+2];
    real_t ux_sol = Usol[0], uy_sol = Usol[1], uz_sol = Usol[2];
    real_t err2 = (ux_sol-ux)*(ux_sol-ux) + (uy_sol-uy)*(uy_sol-uy) + (uz_sol-uz)*(uz_sol-uz);
    L2 += err2;
  }
//  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "gmls", ".", 1, 0, 0, 0.0);
//  silo_file_write_point_cloud(silo, "cantileaver", points);
//  const char* Unames[] = {"ux", "uy", "uz"};
//  silo_file_write_point_field(silo, (const char**)Unames, "cantileaver", U, 3, NULL);
//  silo_file_write_vector_expression(silo, "u", "{<ux>, <uy>, <uz>}");
//  silo_file_close(silo);

  // Scale the L2 error by the uniform volume element.
  real_t dV = 1.0/nx * 1.0/ny;
  L2 = sqrt(L2 * dV);
  log_urgent("L2(error) = %g\n", L2);

  // Super crude hard-wired success metric. :-/
  if (nx == 16)
    assert_true(L2 < 204.6);
  else if (nx == 33)
    assert_true(L2 < 39.8);

  // Clean up.
  int_unordered_set_free(boundary_nodes);
  local_matrix_free(A);
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_gmls_matrix_ctor),
    cmocka_unit_test(test_gmls_matrix_with_frankes_function),
    cmocka_unit_test(test_gmls_matrix_with_cantileaver_beam)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
