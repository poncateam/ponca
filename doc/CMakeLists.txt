
project(doc LANGUAGES CXX VERSION 0.0.1)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  set (PONCA_PROJECT_NUMBER ${GIT_CHANGESET})

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

  add_custom_target(ponca-doc SOURCES ${doc_file_list})

  set(doc_file_list
      ${CMAKE_CURRENT_SOURCE_DIR}/../CHANGELOG
      ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_SOURCE_DIR}/layout.xml
      ${CMAKE_CURRENT_SOURCE_DIR}/html/ponca.css
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ponca.bib
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ponca.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/concepts.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ponca_changelog.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ponca_user_manual.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ponca_module_fitting.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ponca_module_spatialpartitioning.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/ponca_getting_started.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/example.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/example_cxx_basic.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/example_cxx_binding.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/example_cxx_fit_plane.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/example_cxx_pcl.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/example_cu_ssc.mdoc
      ${CMAKE_CURRENT_SOURCE_DIR}/src/example_python_ssc.mdoc
      )
  add_custom_target(ponca-doc-doxygen
                    SOURCES ${doc_file_list}
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                      ${CMAKE_CURRENT_SOURCE_DIR}/images/
                      ${CMAKE_CURRENT_BINARY_DIR}/html/images
                    COMMAND ${CMAKE_COMMAND} -E copy
                      ${CMAKE_CURRENT_SOURCE_DIR}/../CHANGELOG
                      ${CMAKE_CURRENT_BINARY_DIR}/
                    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    DEPENDS
                      ${doc_file_list}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen" VERBATIM )
  add_dependencies(ponca-doc ponca-doc-doxygen)

  install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
      DESTINATION share/doc/
      OPTIONAL )
endif(DOXYGEN_FOUND)
