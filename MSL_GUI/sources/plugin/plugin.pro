#  +-----------------------------------------------------------------------+
#  | MSL GUI - A Graphical User Interface for the Marburg Software Library |
#  |                                                                       |
#  | Copyright (C) 2018 Henning Zickermann                                 |
#  | Contact: <zickermann@mathematik.uni-marburg.de>                       |
#  +-----------------------------------------------------------------------+
#
#    This file is part of MSL GUI.
#
#    MSL GUI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MSL GUI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MSL GUI.  If not, see <https://www.gnu.org/licenses/>.


TEMPLATE        = lib
CONFIG         += plugin
QT             += widgets

CONFIG         += c++11
QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -include ../../throw_assert.h
#QMAKE_CXXFLAGS += -Wno-ignored-qualifiers

TARGET          = mslgui_plugin

win32 {
    CONFIG(debug, release|debug):DESTDIR = ../debug/
    CONFIG(release, release|debug):DESTDIR = ../release/
} else {
    DESTDIR     = ../
}

EXAMPLE_FILES   = main/MslGuiPlugin.json

CONFIG         += object_parallel_to_source

DEFINES        += MSL_GUI


SOURCES += main/mslgui_plugin.cpp \
    main/gui_communicator.cpp \
    methods/parameter_widgetmanager.cpp \
    misc/gui_convergence_logger.cpp \
    parsed/parsedfunction_base.cpp \
    parsed/parsed_1d_function.cpp \
    parsed/parsed_sturm_bvp.cpp \
    raw/sturm_bvp_containers.cpp \
    muparser-2.2.5/src/muParser.cpp \
    muparser-2.2.5/src/muParserBase.cpp \
    muparser-2.2.5/src/muParserBytecode.cpp \
    muparser-2.2.5/src/muParserCallback.cpp \
    muparser-2.2.5/src/muParserError.cpp \
    muparser-2.2.5/src/muParserTokenReader.cpp \
    methods/parameter_widget.cpp \
    main/interval/sturm_bvp/problemtype.cpp \
    main/interval/sturm_bvp/dm_basis.cpp \
    main/interval/poisson_bvp/problemtype.cpp \
    main/cube/elliptic_bvp/problemtype.cpp \
    main/cube/elliptic_bvp/dm_basis.cpp \
    main/cube/poisson_bvp/problemtype.cpp \
    main/ldomain/poisson_bvp/problemtype.cpp \
    main/recring/poisson_bvp/problemtype.cpp \
    main/ldomain/poisson_bvp/dm_qframe.cpp \
    main/recring/poisson_bvp/dm_qframe.cpp \
    main/ldomain/poisson_bvp/dm_aggframe.cpp \
    main/recring/poisson_bvp/dm_aggframe.cpp \
    main/ldomain/biharmonic_bvp/problemtype.cpp \
    main/ldomain/biharmonic_bvp/dm_aggframe.cpp \
    main/interval/sturm_bvp/dm_qframe.cpp \
    main/slitdomain/poisson_bvp/problemtype.cpp \
    main/slitdomain/poisson_bvp/dm_qframe.cpp \
    main/recring/biharmonic_bvp/problemtype.cpp \
    main/recring/biharmonic_bvp/dm_aggframe.cpp \
    main/interval/poisson_bvp/dm_aggframe.cpp \
    main/interval/sturm_bvp/dm_aggframe.cpp \
    main/ldomain/elliptic_bvp/problemtype.cpp \
    main/ldomain/elliptic_bvp/dm_aggframe.cpp \
    main/interval/biharmonic_bvp/problemtype.cpp \
    main/interval/biharmonic_bvp/dm_basis.cpp \
    main/interval/biharmonic_bvp/dm_aggframe.cpp \
    main/cube/elliptic_bvp/dm_qframe.cpp


HEADERS  += discr/abstract_discretized_problem.h \
    discr/discretization_module_base.h \
    discr/discretization_module_guidata.h \
    discr/generic_discretized_problem.h \
    discr/generic_basis_discretization_module.h \
    main/generic_problemtype_module.h \
    main/gui_communicator.h \
    main/mslgui_plugin.h \
    methods/parameter_widgetmanager.h \
    misc/cached_frame_problem_helper.h \
    misc/convert_basis_to_qstring.h \
    misc/gui_convergence_logger.h \
    misc/string_conversion.h \
    misc/typelist.h \
    parsed/parsed_1d_function.h \
    parsed/parsed_expression_checker.h \
    parsed/parsed_function.h \
    parsed/parsed_sturm_bvp.h \
    parsed/parsedfunction_base.h \
    raw/elliptic_bvp_containers.h \
    raw/raw_problem_container.h \
    raw/raw_problem_module.h \
    raw/sturm_bvp_containers.h \
    solution/abstract_solution.h \
    solution/generic_solution.h \
    solution/solution_tools.h \
    methods/basis/cdd1.h \
    methods/method_base.h \
    methods/parameter_widget.h \
    methods/basis/list.h \
    misc/instantiate.h \
    main/interval/sturm_bvp/problemtype.h \
    main/interval/sturm_bvp/dm_basis.h \
    main/interval/poisson_bvp/problemtype.h \
    main/interval/poisson_bvp/dm_basis.h \
    main/cube/elliptic_bvp/problemtype.h \
    main/cube/elliptic_bvp/dm_basis.h \
    main/cube/poisson_bvp/problemtype.h \
    main/cube/poisson_bvp/dm_basis.h \
    discr/generic_qframe_discretization_module.h \
    main/ldomain/poisson_bvp/dm_qframe.h \
    methods/qframe/cdd2_quarklet.h \
    methods/qframe/list.h \
    main/ldomain/poisson_bvp/problemtype.h \
    main/recring/poisson_bvp/problemtype.h \
    main/recring/poisson_bvp/dm_qframe.h \
    misc/generate_has_member_type.h \
    misc/basis_or_frame.h \
    main/ldomain/poisson_bvp/dm_aggframe.h \
    discr/generic_aggframe_discretization_module.h \
    methods/aggframe/list.h \
    methods/aggframe/steepest_descent_solver.h \
    methods/aggframe/implementation_steepest_descent.h \
    methods/qframe/implementation_cdd2_quarklet.h \
    methods/aggframe/additive_schwarz_solver.h \
    methods/aggframe/implementation_additive_schwarz.h \
    main/recring/poisson_bvp/dm_aggframe.h \
    methods/aggframe/implementation_multiplicative_schwarz.h \
    methods/aggframe/multiplicative_schwarz_solver.h \
    main/ldomain/biharmonic_bvp/problemtype.h \
    raw/biharmonic_bvp_containers.h \
    main/ldomain/biharmonic_bvp/dm_aggframe.h \
    main/interval/sturm_bvp/dm_qframe.h \
    main/interval/poisson_bvp/dm_qframe.h \
    methods/qframe/implementation_steepest_descent_ks_quarklet.h \
    methods/qframe/steepest_descent_ks_quarklet.h \
    methods/basis/stevenson_awgm.h \
    main/interval/sturm_bvp/sturm_qframe_equation.h \
    main/slitdomain/poisson_bvp/problemtype.h \
    main/slitdomain/poisson_bvp/dm_qframe.h \
    main/recring/biharmonic_bvp/problemtype.h \
    main/recring/biharmonic_bvp/dm_aggframe.h \
    main/interval/poisson_bvp/dm_aggframe.h \
    methods/basis/implementation_cdd2.h \
    methods/basis/cdd2.h \
    methods/aggframe/implementation_richardson.h \
    methods/aggframe/richardson_solver.h \
    main/interval/sturm_bvp/dm_aggframe.h \
    main/ldomain/elliptic_bvp/problemtype.h \
    main/ldomain/elliptic_bvp/dm_aggframe.h \
    main/interval/biharmonic_bvp/problemtype.h \
    main/interval/biharmonic_bvp/dm_basis.h \
    main/interval/biharmonic_bvp/dm_aggframe.h \
    main/cube/elliptic_bvp/dm_qframe.h \
    main/cube/poisson_bvp/dm_qframe.h \
    misc/interval_bases_list.h \
    misc/interval_qframes_list.h


INCLUDEPATH += ../ \
               ../../../ \
               ../../../MathTL \
               ../../../WaveletTL \
               ../../../FrameTL \
               muparser-2.2.5/include
