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
      [ "Object-space Curvature using Cuda/C++ and KdTree queries", "example_cu_kdtree_page.html", [
        [ "Introduction", "example_cu_kdtree_page.html#cu_kdtree_intro_sec", null ],
        [ "Source code", "example_cu_kdtree_page.html#cu_kdtree_code", null ]
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
    [ "The ponca ecosystem", "poncaecosysrem.html", null ],
    [ "Getting started", "ponca_getting_started_page.html", [
      [ "Requirements", "ponca_getting_started_page.html#getting_started_requirements_sec", null ],
      [ "Download", "ponca_getting_started_page.html#getting_started_download_sec", null ],
      [ "Installation", "ponca_getting_started_page.html#getting_started_installation_sec", [
        [ "Use as cmake subdirectory", "ponca_getting_started_page.html#getting_started_installation_cmakesubdirectory_subsec", null ],
        [ "Compilation", "ponca_getting_started_page.html#getting_started_installation_compilation_subsec", null ],
        [ "Cmake package", "ponca_getting_started_page.html#getting_started_installation_cmake_subsec", null ],
        [ "CPM", "ponca_getting_started_page.html#getting_started_installation_cmake_cpm", null ]
      ] ],
      [ "First steps", "ponca_getting_started_page.html#getting_started_first_step_sec", null ]
    ] ],
    [ "User Manual", "user_manual_page.html", [
      [ "Defining Points in Ponca", "ponca_points.html", [
        [ "Points", "ponca_points.html#ponca_points_intro", null ],
        [ "Utility classes and functions", "ponca_points.html#points_utility", null ]
      ] ],
      [ "Fitting Module: User Manual", "fitting.html", [
        [ "Introduction", "fitting.html#fitting_intro", null ],
        [ "First Steps", "fitting.html#fitting_firstSteps", [
          [ "Include directives", "fitting.html#fitting_codeStructure", null ],
          [ "Definition of the Fitting object", "fitting.html#fitting_Define", null ],
          [ "Fitting Process", "fitting.html#fitting_Fitting", null ],
          [ "Check fitting status", "fitting.html#fitting_Checkstatus", null ],
          [ "Basic Outputs", "fitting.html#fitting_outputs", null ]
        ] ],
        [ "Evaluation schemes and projection", "fitting.html#evaluation_schemes", null ],
        [ "Computing derivatives", "fitting.html#fitting_derivatives", null ],
        [ "Changing fitting process", "fitting.html#fitting_custom", null ],
        [ "Defining a new neighbor filter", "fitting.html#fitting_newfilter", [
          [ "Defining a new distance kernel", "fitting.html#fitting_newkernel", null ],
          [ "Filter API", "fitting.html#fitting_newfilterapi", null ]
        ] ],
        [ "Defining a new Estimator", "fitting.html#fitting_newestimator", [
          [ "Understanding CRTP in Ponca", "fitting.html#crtpestimator", null ],
          [ "Estimator API", "fitting.html#fitting_newestim", [
            [ "Minimal requirements", "fitting.html#fitting_newestimbase", null ],
            [ "Computational objets capabilities and requirements", "fitting.html#fitting_newcapabilities", null ],
            [ "Providing cast operations", "fitting.html#fitting_cast", null ]
          ] ]
        ] ]
      ] ],
      [ "Fitting Module: Reference Manual", "fittingreference.html", [
        [ "Fitting techniques Overview", "fittingreference.html#fittingreference_primitiveoverview", null ],
        [ "Capabilities of Fitting tools", "fittingreference.html#fittingreference_outputandcapabilities", null ],
        [ "Computing Curvatures", "fittingreference.html#fitting_curvature", null ],
        [ "Available filters", "fittingreference.html#fitting_filters", null ],
        [ "Available Weighting kernel", "fittingreference.html#fitting_kernel", null ]
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
          [ "Usage in Cuda kernels", "spatialpartitioning.html#spatialpartitioning_kdtree_cuda", null ],
          [ "Customizing the KdTree using <tt>Traits</tt>", "spatialpartitioning.html#spatialpartitioning_kdtree_extending", null ],
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
        [ "Typedefs", "functions_type.html", "functions_type" ],
        [ "Enumerator", "functions_eval.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"classPonca_1_1GLSDer.html#adbd5c9acc1515c79390d6c1962ef5d63",
"classPonca_1_1NormalDerivativeWeingartenEstimator.html#a4cf8ebc202d9d0aa6ac0ce7f0100e8fe",
"classPonca_1_1UnorientedSphereFitImpl.html#a36ee49b0dc6a829720f57c2d7962b3d2",
"structPonca_1_1ComputeObject.html#af4490afa97b2845366ee339c56661b5d"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';