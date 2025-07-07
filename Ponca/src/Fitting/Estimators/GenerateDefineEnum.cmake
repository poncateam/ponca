set(ENUM_FITS "PSS ASO APSS PCA Sphere" CACHE STRING "The fit types that the estimator factory will make available")

message(STATUS "Compiling Estimator Factory for the following FIT TYPES : ${ENUM_FITS}")



# Transform the list into lines like: ENUM_FIT(name)
# This assumes names are space-separated
string(REPLACE " " ")\\\nENUM_FIT(" ENUM_FITS_CONTENT "${ENUM_FITS}")
set(ENUM_FITS_CONTENT "ENUM_FIT(${ENUM_FITS_CONTENT})")

# Configure header file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Ponca/src/Fitting/Estimators/defineEnum.h.in
        ${CMAKE_CURRENT_SOURCE_DIR}/Ponca/src/Fitting/Estimators/defineEnum.h @ONLY)