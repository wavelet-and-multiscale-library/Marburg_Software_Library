/*  -*- c++ -*-

   +-----------------------------------------------------------------------+
   | MSL GUI - A Graphical User Interface for the Marburg Software Library |
   |                                                                       |
   | Copyright (C) 2018 Henning Zickermann                                 |
   | Contact: <zickermann@mathematik.uni-marburg.de>                       |
   +-----------------------------------------------------------------------+

     This file is part of MSL GUI.

     MSL GUI is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     MSL GUI is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with MSL GUI.  If not, see <https://www.gnu.org/licenses/>.
*/


#ifndef ABSTRACT_GUI_COMMUNICATOR_H
#define ABSTRACT_GUI_COMMUNICATOR_H

#include <QObject>

#include "computation_state_enums.h"

struct GuiInputData;
class AbstractProblemTypeModule;



class AbstractGuiCommunicator : public QObject
{
public:
    virtual ~AbstractGuiCommunicator() {}

    virtual bool gotAbortRequest() const = 0;

signals:
    virtual void coeffComputationEnded(int computationNo, CoeffComputationState::Enum endState) const = 0;
    virtual void plotComputationEnded(int computationNo, PlotComputationState::Enum endState, int resolution) const = 0;

    virtual void convergencePlotDataPairComputed(double supportsize, double error) const = 0;
    virtual void convergenceTimePlotDataPairComputed(double seconds, double error) const = 0;

    virtual void solutionSaved(int computationNo, const QStringList& availableOptionalLogs) const = 0;

    virtual void matrixNormsComputed(int computationNo, double norm_A, double norm_Ainv) const = 0;

    /*
     * xValues and yValues are C-arrays which hold the x- and y-Values of a solution plot part on a single patch.
     * samplingPointCount is the number of (x,y)-data pairs, hence the length of the C-arrays xValues and yValues.
     * newPlot indicates wether the data belongs to (the first patch of) a new solution plot for the given
     * computation number or to a plot for that number which in parts has already been delivered to the GUI.
     */
    virtual void dataFor1DSolutionPlotComputed(int computationNo, const double* xValues,
                                               const double* yValues, int samplingPointCount, bool newPlot) const = 0;
    /*
     * Data for a 2D solution plot patch is supposed to be laid out on a rectangular grid.
     * Then xValues, yValues and zValues are C-arrays which hold the x-, y- and z-Values
     * in column major ordering on that grid (rows run in x-direction, columns in y-direction).
     * sampleCountX and sampleCountY are the numbers of different(!) x-values and y-Values on the grid,
     * hence sampleCountX is the number of columns and sampleCountY is the number of rows of that grid.
     * newPlot indicates wether the data belongs to (the first patch of) a new solution plot for the given
     * computation number or to a plot for that number which in parts has already been delivered to the GUI.
     */
    virtual void dataFor2DSolutionPlotComputed(int computationNo, const double* xValues, const double* yValues,
                                               const double* zValues, int sampleCountX, int sampleCountY, bool newPlot) const = 0;

    virtual void statusMessageGenerated(const QString& str) const = 0;
    virtual void errorOccured(const QString& errorMessage) const = 0;

public slots:
    virtual void computeSolutionCoeffs(AbstractProblemTypeModule* problemType, const GuiInputData& input) = 0;
    virtual void requestAbort() = 0;

    virtual void computeSolutionPlot(int computationNo, int resolution) = 0;
    virtual void exportSolutionPlotToMatlabFile(int computationNo, const QString& filename) const = 0;
    virtual void exportIndexPlotToMatlabFile(int computationNo, const QString& filename) const = 0;

    virtual void exportConvergenceLogsToMatlabFile(int computationNo, const QString& filename) const = 0;
    virtual void exportOptionalLogToMatlabFile(int computationNo, int logIndex,
                                               const QString& filename) const = 0;
    virtual void deleteSolution(int computationNo) = 0;
};


#endif // ABSTRACT_GUI_COMMUNICATOR_H
