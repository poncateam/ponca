/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "Ponca", "index.html", [
    [ "About Ponca", "index.html#ponca_about_sec", null ],
    [ "Basic FAQs", "index.html#ponca_instructions_sec", null ],
    [ "Credits", "index.html#ponca_credits_sec", [
      [ "Developers and contributors", "index.html#ponca_credits_subsection_crew", null ],
      [ "Patate/Grenaille crew and contributors", "index.html#ponca_credits_subsection_grenaille", null ],
      [ "Citation", "index.html#ponca_credits_subsection_citation", null ],
      [ "Users", "index.html#ponca_credits_subsection_users", null ]
    ] ],
    [ "Examples", "example_page.html", [
      [ "Ponca basic CPU", "example_cxx_basic_page.html", null ],
      [ "Ponca basic plane fit", "example_cxx_fit_plane_page.html", null ],
      [ "Ponca basic line fit", "example_cxx_fit_line_page.html", null ],
      [ "Ponca data-structure binding", "example_cxx_binding_page.html", null ],
      [ "Using Ponca to compute surface curvature in PCL", "example_cxx_pcl_page.html", [
        [ "Introduction", "example_cxx_pcl_page.html#pcl_intro_sec", [
          [ "Usage", "example_cxx_pcl_page.html#pcl_usage_sec", null ]
        ] ],
        [ "PCL wrapper", "example_cxx_pcl_page.html#pcl_wrapper_sec", [
          [ "Datastructures declaration", "example_cxx_pcl_page.html#pcl_wrapper_sec_pcl_wrapper_sec_h", null ],
          [ "Datastructures implementation", "example_cxx_pcl_page.html#pcl_wrapper_sec_pcl_wrapper_sec_hpp", null ],
          [ "Datastructures instanciation", "example_cxx_pcl_page.html#pcl_wrapper_sec_pcl_wrapper_sec_cpp", null ],
          [ "Main program", "example_cxx_pcl_page.html#pcl_main_sec", null ]
        ] ]
      ] ],
      [ "Screen Space Curvature using Cuda/C++", "example_cu_ssc_page.html", [
        [ "Introduction", "example_cu_ssc_page.html#cu_ssgl_intro_sec", [
          [ "Installation and usage", "example_cu_ssc_page.html#cu_ssgl_sec_dep_subsec", null ]
        ] ],
        [ "Cuda programming", "example_cu_ssc_page.html#cu_ssgl_cuda_sec", [
          [ "Define fitting data structure", "example_cu_ssc_page.html#cu_ssgl_cuda_mypoint_sec", null ],
          [ "Define weighting functions", "example_cu_ssc_page.html#cu_ssgl_cuda_weight_sec", null ],
          [ "Define fitting primitive", "example_cu_ssc_page.html#cu_ssgl_cuda_fit_sec", null ],
          [ "Kernel", "example_cu_ssc_page.html#cu_ssgl_cuda_kernel_sec", null ],
          [ "Memory access", "example_cu_ssc_page.html#cu_ssgl_cuda_access_sec", null ]
        ] ],
        [ "The whole code", "example_cu_ssc_page.html#cu_ssgl_sec", null ]
      ] ],
      [ "Screen Space Curvature using Cuda and Python", "example_python_ssc_page.html", [
        [ "Introduction", "example_python_ssc_page.html#pyssgl_intro_sec", [
          [ "Installation and usage", "example_python_ssc_page.html#pyssgl_intro_sec_dep_subsec", null ]
        ] ],
        [ "Cuda programming", "example_python_ssc_page.html#pyssgl_cuda_sec", [
          [ "Define fitting data structure", "example_python_ssc_page.html#pyssgl_cuda_mypoint_sec", null ],
          [ "Define weighting functions", "example_python_ssc_page.html#pyssgl_cuda_weight_sec", null ],
          [ "Define fitting primitive", "example_python_ssc_page.html#pyssgl_cuda_fit_sec", null ],
          [ "Kernel", "example_python_ssc_page.html#pyssgl_cuda_kernel_sec", null ],
          [ "Memory access", "example_python_ssc_page.html#pyssgl_cuda_access_sec", null ]
        ] ],
        [ "Python script", "example_python_ssc_page.html#pyssgl_python_sec", null ]
      ] ],
      [ "Ponca::KdTree neighbor searches", "example_cxx_neighbor_search.html", null ],
      [ "Comparison between Nanoflann and Ponca KdTree APIs", "example_cxx_nanoflann_page.html", [
        [ "Introduction", "example_cxx_nanoflann_page.html#nanoflann_intro_sec", [
          [ "Compilation", "example_cxx_nanoflann_page.html#nanoflann_compilation_sec", null ]
        ] ],
        [ "API comparisons", "example_cxx_nanoflann_page.html#nanoflann_comparison_sec", null ],
        [ "Timings", "example_cxx_nanoflann_page.html#nanoflann_timings_sec", null ],
        [ "Example source code", "example_cxx_nanoflann_page.html#nanoflann_sourcecode_sec", null ]
      ] ]
    ] ],
    [ "Releases overview", "ponca_changelog.html", [
      [ "Complete Changelog", "ponca_changelog.html#ponca_changelog_sec", null ]
    ] ],
    [ "Getting started", "ponca_getting_started_page.html", [
      [ "Requirements", "ponca_getting_started_page.html#getting_started_requirements_sec", null ],
      [ "Download", "ponca_getting_started_page.html#getting_started_download_sec", null ],
      [ "Installation", "ponca_getting_started_page.html#getting_started_installation_sec", [
        [ "Use as cmake subdirectory", "ponca_getting_started_page.html#getting_started_installation_cmakesubdirectory_subsec", null ],
        [ "Compilation", "ponca_getting_started_page.html#getting_started_installation_compilation_subsec", null ],
        [ "Cmake package", "ponca_getting_started_page.html#getting_started_installation_cmake_subsec", null ]
      ] ],
      [ "First steps", "ponca_getting_started_page.html#getting_started_first_step_sec", null ]
    ] ],
    [ "User Manual", "user_manual_page.html", [
      [ "Ponca Concepts", "ponca_concepts.html", null ],
      [ "Fitting Module: User Manual", "fitting.html", [
        [ "Introduction", "fitting.html#fitting_intro", [
          [ "Design choices and programing techniques", "fitting.html#fitting_design", null ],
          [ "Fitting primitives and compute neighborhood properties", "fitting.html#fitting_primitiveOverview", null ],
          [ "Fitting techniques overview", "fitting.html#fitting_availableFunctionalities", null ],
          [ "Structure of the documentation", "fitting.html#fitting_dicStructure", null ]
        ] ],
        [ "First Steps", "fitting.html#fitting_firstSteps", [
          [ "Include directives", "fitting.html#fitting_codeStructure", null ],
          [ "Data Samples", "fitting.html#fitting_datas", null ],
          [ "Definition of the Fitting object", "fitting.html#fitting_Define", null ],
          [ "Fitting Process", "fitting.html#fitting_Fitting", null ],
          [ "Check fitting status", "fitting.html#fitting_Checkstatus", null ],
          [ "Basic Outputs", "fitting.html#fitting_outputs", null ]
        ] ],
        [ "Advanced usage", "fitting.html#fitting_advanced", [
          [ "Computing derivatives", "fitting.html#fitting_derivatives", null ],
          [ "Computational objets, basket and CRTP", "fitting.html#fitting_extensions_deps", null ],
          [ "Computational objets capabilities and requirements", "fitting.html#fitting_capabilities", null ],
          [ "Sharing computations between fits", "fitting.html#fitting_multiprimitive", null ],
          [ "Computing Curvatures", "fitting.html#fitting_cuvature", null ],
          [ "Cuda", "fitting.html#fitting_cuda", null ]
        ] ],
        [ "Fitting module: Concepts", "fitting_concepts.html", [
          [ "API of Computational Objects", "fitting_concepts.html#concepts_computObject", [
            [ "Objects used in Basket", "fitting_concepts.html#concepts_computObjectBasket", null ],
            [ "Objects used in BasketDiff", "fitting_concepts.html#concepts_computObjectBasketDiff", null ]
          ] ],
          [ "Concepts related to weighting functions", "fitting_concepts.html#concepts_weighting", null ]
        ] ]
      ] ],
      [ "Fitting module: Concepts", "fitting_concepts.html", [
        [ "API of Computational Objects", "fitting_concepts.html#concepts_computObject", [
          [ "Objects used in Basket", "fitting_concepts.html#concepts_computObjectBasket", null ],
          [ "Objects used in BasketDiff", "fitting_concepts.html#concepts_computObjectBasketDiff", null ]
        ] ],
        [ "Concepts related to weighting functions", "fitting_concepts.html#concepts_weighting", null ]
      ] ],
      [ "Spatial Partitioning: User Manual", "spatialpartitioning.html", [
        [ "Introduction", "spatialpartitioning.html#spatialpartitioning_intro", [
          [ "Datastructures", "spatialpartitioning.html#spatialpartitioning_intro_structures", null ],
          [ "Queries", "spatialpartitioning.html#spatialpartitioning_intro_queries", null ]
        ] ],
        [ "KdTree", "spatialpartitioning.html#spatialpartitioning_kdtree", [
          [ "Specifications", "spatialpartitioning.html#spatialpartitioning_kdtree_implementation", null ],
          [ "Basic usage", "spatialpartitioning.html#spatialpartitioning_kdtree_usage", [
            [ "Construction", "spatialpartitioning.html#spatialpartitioning_kdtree_usage_construction", null ],
            [ "Queries", "spatialpartitioning.html#spatialpartitioning_kdtree_usage_queries", null ],
            [ "Samples and indexing", "spatialpartitioning.html#spatialpartitioning_kdtree_usage_samples_and_indexing", null ]
          ] ],
          [ "Extending KdTree", "spatialpartitioning.html#spatialpartitioning_kdtree_extending", null ],
          [ "Usage of the convenience classes KdTree and KdTreeBase", "spatialpartitioning.html#spatialpartitioning_kdtree_usage_which_class", null ]
        ] ],
        [ "KnnGraph", "spatialpartitioning.html#spatialpartitioning_knngraph", [
          [ "Basic usage", "spatialpartitioning.html#spatialpartitioning_knngraph_usage", [
            [ "Construction", "spatialpartitioning.html#spatialpartitioning_knngraph_usage_construction", null ],
            [ "Queries", "spatialpartitioning.html#spatialpartitioning_knngraph_usage_queries", null ]
          ] ]
        ] ]
      ] ],
      [ "Common module: User Manual", "common.html", [
        [ "Overview of the proposed tools", "common.html#common_intro", null ]
      ] ]
    ] ],
    [ "Bibliography", "citelist.html", null ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", null ],
        [ "Typedefs", "functions_type.html", null ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Enumerator", "functions_eval.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"classPonca_1_1KdTreeNearestIterator.html",
"classPonca_1_1PolynomialSmoothWeightKernel.html#a57638efcd9b72fa08be0b3a9ff385e47",
"pages.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';