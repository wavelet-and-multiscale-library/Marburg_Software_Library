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


QT       += core gui charts datavisualization xml

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11
QMAKE_CXXFLAGS += -std=gnu++11

TARGET = msl_gui

win32 {
    CONFIG(debug, release|debug):DESTDIR = ../debug/
    CONFIG(release, release|debug):DESTDIR = ../release/
} else {
    DESTDIR    = ../
}

OBJECTS_DIR = obj
MOC_DIR = moc
RCC_DIR = rcc
UI_DIR = ui

TEMPLATE = app


SOURCES += src/main.cpp \
           src/mainwindow.cpp \
           src/epsilon_spinbox.cpp \
           src/convergence_chart.cpp \
           src/computation_table.cpp \
           src/computation_manager.cpp \
           src/export_combobox.cpp \
           src/qoutput_stream.cpp \
           src/plot_manager.cpp \
           src/convergence_chartview.cpp \
           src/matrixnorm_dialog.cpp \
           src/problem_definition_manager.cpp \
           src/problem_definition_dialog.cpp \
           src/computation_status_box.cpp \
           src/info_defining_dialog.cpp \
           src/info_controls_dialog.cpp \
    src/info_about_dialog.cpp

HEADERS  += src/mainwindow.h \
            interfaces/abstract_problemtype_module.h \
            interfaces/abstract_gui_communicator.h \
            interfaces/gui_inputdata.h \
            interfaces/mslgui_plugin_interface.h \
            interfaces/abstract_parameterwidget_manager.h \
            src/epsilon_spinbox.h \
            src/convergence_chart.h \
            src/computation_table.h \
            src/computation_manager.h \
            src/export_combobox.h \
            src/qoutput_stream.h \
            interfaces/computation_state_enums.h \
            src/plot_manager.h \
            src/convergence_chartview.h \
            src/matrixnorm_dialog.h \
            src/problem_definition_manager.h \
            src/problem_definition_dialog.h \
            src/computation_status_box.h \
            src/info_defining_dialog.h \
            src/info_controls_dialog.h \
    src/info_about_dialog.h

FORMS    += ui/mainwindow.ui \
            ui/matrixnorm_dialog.ui \
            ui/problem_definition_dialog.ui \
            ui/info_defining_dialog.ui \
            ui/info_controls_dialog.ui \
            ui/info_about_dialog.ui

RESOURCES += resources/resources.qrc
