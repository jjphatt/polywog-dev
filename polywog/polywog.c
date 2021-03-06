// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file really only checks that the version of polywog matches that 
// of polymec.

#include "core/polymec_version.h"
#include "polywog/polywog.h"

#if POLYWOG_MAJOR_VERSION > 0

#if POLYMEC_MAJOR_VERSION != POLYWOG_MAJOR_VERSION
#error "The installed major version of polymec differs from that of polywog. Please make sure these versions match."
#endif

#endif

