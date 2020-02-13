###############################################################################
# Grenaille GL Viewer                                                         #
###############################################################################


find_package(OpenGL)
find_package(Qt5Widgets)
find_package(Qt5OpenGL)

set(grenaille_glviewer_MISSING "")
if(NOT OPENGL_FOUND)
    set(grenaille_glviewer_MISSING "${grenaille_glviewer_MISSING} OpenGL")
endif()
if(NOT Qt5Widgets_FOUND)
    set(grenaille_glviewer_MISSING "${grenaille_glviewer_MISSING} Qt5Widgets")
endif()
if(NOT Qt5OpenGL_FOUND)
    set(grenaille_glviewer_MISSING "${grenaille_glviewer_MISSING} Qt5OpenGL")
endif()


if(NOT grenaille_glviewer_MISSING)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)

    # Remove a warning triggred by Qt.
    if(CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "-Wno-long-long ${CMAKE_CXX_FLAGS}")
    endif(CMAKE_COMPILER_IS_GNUCXX)

    include_directories(
        ${EIGEN3_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR}
        "../../.."
    )

    # Add the include directories for the Qt 5 Widgets module to
    # the compile lines. Mandatory on some old plateforms
    include_directories(
        ${Qt5Widgets_INCLUDE_DIRS}
        ${Qt5OpenGL_INCLUDE_DIRS}
    )

    set(CMAKE_INCLUDE_CURRENT_DIR ON)

    set(grenaille_glviewer_UIS
      ${CMAKE_CURRENT_SOURCE_DIR}/mainwindow.ui
    )

    set(grenaille_glviewer_HEADERS
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer.h
      ${CMAKE_CURRENT_SOURCE_DIR}/fittingmanager.h
      ${CMAKE_CURRENT_SOURCE_DIR}/mainwindow.h
    )

    set(grenaille_glviewer_SRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/fittingmanager.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/mainwindow.cpp
    )


    # GLSL files
    set(SHADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/shaders)
    file(GLOB SHADER_VERT_FILES ${SHADERS_DIR}/*.vert)
    file(GLOB SHADER_FRAG_FILES ${SHADERS_DIR}/*.frag)

    # Assets
    set(ASSETS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/assets)

    add_custom_target( install_res
        COMMENT "Copying shaders and other ressources" VERBATIM
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${SHADERS_DIR} ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${ASSETS_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_executable(grenaille_glviewer
                          ${grenaille_glviewer_SRCS}
                          ${grenaille_glviewer_HEADERS}
                          ${grenaille_glviewer_UIS}
                          ${SHADER_VERT_FILES}
                          ${SHADER_FRAG_FILES}
    )
    add_dependencies(grenaille_glviewer install_res)
    add_dependencies(grenaille_glviewer patate_common Grenaille)
    add_dependencies(grenaille_examples grenaille_glviewer)

    target_link_libraries(grenaille_glviewer
                          ${OPENGL_LIBRARIES}
                          Qt5::Widgets
                          Qt5::OpenGL
    )
else()
    message(WARNING "Skipping grenaille_glviewer because of the following missing dependencies:${grenaille_glviewer_MISSING}")
endif()
