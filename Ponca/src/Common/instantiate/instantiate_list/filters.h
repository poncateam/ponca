/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#if defined(_PONCA_INSTANTIATE_SMOOTHWEIGHT) || defined(_PONCA_INSTANTIATE_ALL) || defined(_PONCA_INSTANTIATE_DEFAULT)
#    define _NF DistWeightFilter<_P, SmoothWeightKernel<_S>>
#    include "baskets.h"
#    undef _NF
#endif

#if defined(_PONCA_INSTANTIATE_CONSTANTWEIGHT) || defined(_PONCA_INSTANTIATE_ALL)
#    define _NF DistWeightFilter<_P, ConstantWeightKernel<_S>>
#    include "baskets.h"
#    undef _NF
#endif

#if defined(_PONCA_INSTANTIATE_NOWEIGHT) || defined(_PONCA_INSTANTIATE_ALL)
#    define _NF Ponca::NoWeightFilter<_P>
#    include "baskets.h"
#    undef _NF
#endif
