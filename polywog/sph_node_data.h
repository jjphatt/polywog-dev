// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_SPH_NODE_DATA_H
#define POLYWOG_SPH_NODE_DATA_H

#include "core/point.h"
#include "core/tensor2.h"

// The SPH node data object is a Plain Old Datatype (POD) that can hold data
// pertaining to an SPH node and its neighborhood.
typedef struct 
{
  // Moments.
  real_t        zeroth_moment;  // kernel sum minus self contribution.
  vector_t      first_moment;
  sym_tensor2_t second_moment;
} sph_node_data_t;

#endif
