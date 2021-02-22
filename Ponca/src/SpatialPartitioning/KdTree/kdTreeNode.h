/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../defines.h"

namespace Ponca {
template<typename Scalar>
struct KdTreeNode
{
    union {
        struct {
			Scalar       splitValue;
            unsigned int   firstChildId:24;
            unsigned int   dim:2;
            unsigned int   leaf:1;
        };
        struct {
            unsigned int   start;
            unsigned short size;
        };
    };
};

}   
