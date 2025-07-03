#pragma once

// The different fit types are compiled depending on if they are defined through this macro or not
#define ENUM_FITS \
ENUM_FIT(PSS) \
ENUM_FIT(ASO) \
ENUM_FIT(APSS) \
ENUM_FIT(PCA) \
ENUM_FIT(Sphere) \


////////////////////////////////////////////////////////////

namespace Estimators {

    /*! \brief Enum of every available fit type */
    enum class FitType {
#define ENUM_FIT(type) type,
        ENUM_FITS
#undef ENUM_FIT
    };

    /*! \brief Count the number of different available fit types */
    inline constexpr int NUMBER_OF_FIT_TYPES = 0
#define ENUM_FIT(type) +1
        ENUM_FITS
#undef ENUM_FIT
    ;
}