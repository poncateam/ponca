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
        ${GLEW_INCLUDE_DIRS}
        "../../.."
    )

    set(CMAKE_INCLUDE_CURRENT_DIR ON)

    set(grenaille_glviewer_UIS
      mainwindow.ui
    )

    set(grenaille_glviewer_HEADERSS
      viewer.h
      mainwindow.h
    )

    set(grenaille_glviewer_SRCS
      main.cpp
      viewer.cpp
      mainwindow.cpp
    )

    add_executable(grenaille_glviewer
                          ${grenaille_glviewer_SRCS}
                          ${grenaille_glviewer_HEADERS}
                          ${grenaille_glviewer_UIS})
    add_dependencies(grenaille_examples grenaille_glviewer)
    add_dependencies(grenaille_glviewer patate_common Grenaille)

    target_link_libraries(grenaille_glviewer
                          ${OPENGL_LIBRARIES}
                          Qt5::Widgets
                          Qt5::OpenGL
    )
else()
    message(WARNING "Skipping grenaille_glviewer because of the following missing dependencies:${grenaille_glviewer_MISSING}")
endif()
