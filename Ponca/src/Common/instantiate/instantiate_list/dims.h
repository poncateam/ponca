/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#if defined(_PONCA_INSTANTIATE_2D) || defined(_PONCA_INSTANTIATE_ALL)
#    define _D 2
#    include "points.h"
#    undef _D
#endif

#if defined(_PONCA_INSTANTIATE_3D) || defined(_PONCA_INSTANTIATE_ALL)
#    define _D 3
#    include "points.h"
#    undef _D
#endif
