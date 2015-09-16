// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_MLS_SHAPE_FUNCTION_H
#define POLYWOG_MLS_SHAPE_FUNCTION_H

#include "core/point_cloud.h"
#include "model/stencil.h"
#include "polywog/shape_function.h"

// Creates a moving-least-squares polynomial shape function defined on the 
// given point cloud with neighborhoods given by the given stencil. Here,
// smoothing_lengths is a field (array) assigning a characteristic extent, h, 
// to each point in the domain, and must have enough storage for ghost 
// points.
shape_function_t* mls_shape_function_new(int polynomial_order,
                                         shape_function_kernel_t* kernel,
                                         point_cloud_t* domain,
                                         stencil_t* neighborhoods,
                                         real_t* smoothing_lengths);

#endif
