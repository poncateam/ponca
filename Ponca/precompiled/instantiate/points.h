/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#if defined(_PONCA_INSTANTIATE_POINTPOSITION) || defined(_PONCA_INSTANTIATE_ALL)
#    define _P PointPosition<_S, _D>
#    include "filters.h"
#    undef _P
#endif

#if defined(_PONCA_INSTANTIATE_POINTPOSITIONNORMAL) || defined(_PONCA_INSTANTIATE_ALL)
#    define _P PointPositionNormal<_S, _D>
#    define _Normal
#    include "filters.h"
#    undef _Normal
#    undef _P
#endif
