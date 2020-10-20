#pragma once

#include "../defines.h"

namespace Ponca {

struct KdTreeNode
{
    union {
        struct {
            SPScalar       splitValue;
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

} // namespace pca
