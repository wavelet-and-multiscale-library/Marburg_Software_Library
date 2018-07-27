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


#ifndef GUI_COMMUNICATOR_H
#define GUI_COMMUNICATOR_H

#include <memory>

#include "GUI/interfaces/abstract_gui_communicator.h"
#include "MathTL/geometry/sampled_mapping.h"
#include "discr/abstract_discretized_problem.h"
#include "solution/abstract_solution.h"


class AbstractProblemTypeModule;


class GuiCommunicator : public AbstractGuiCommunicator
{
    Q_OBJECT
public:
    GuiCommunicator();

    void sendPlotDataToGui(const MathTL::SampledMapping<1>& plotData, int computationNo, bool newPlot = true) const;
    void sendPlotDataToGui(const MathTL::SampledMapping<2>& plotData, int computationNo, bool newPlot = true) const;
    void sendPlotDataToGui(const MathTL::Array1D< MathTL::SampledMapping<1> >& plotData, int computationNo) const;
    void sendPlotDataToGui(const MathTL::Array1D< MathTL::SampledMapping<2> >& plotData, int computationNo) const;

// Prescribed by AbstractGuiCommunicator:

    bool gotAbortRequest() const override;

signals:
    void coeffComputationEnded(int computationNo, CoeffComputationState::Enum endState) const override;
    void plotComputationEnded(int computationNo, PlotComputationState::Enum endState, int resolution) const override;

    void convergencePlotDataPairComputed(double supportsize, double error) const override;
    void convergenceTimePlotDataPairComputed(double seconds, double error) const override;

    void solutionSaved(int computationNo, const QStringList& availableOptionalLogs) const override;

    void matrixNormsComputed(int computationNo, double norm_A, double norm_Ainv) const override;

    void dataFor1DSolutionPlotComputed(int computationNo, const double* xValues, const double* yValues,
                                       int sampleCount, bool newPlot) const override;
    void dataFor2DSolutionPlotComputed(int computationNo, const double* xValues, const double* yValues,
                                       const double* zValues, int sampleCountX, int sampleCountY, bool newPlot) const override;

    void statusMessageGenerated(const QString& str) const override;
    void errorOccured(const QString& errorMessage) const override;

public slots:
    void computeSolutionCoeffs(AbstractProblemTypeModule* problemType, const GuiInputData& input) override;
    void requestAbort() override;

    void computeSolutionPlot(int computationNo, int resolution) override;
    void exportSolutionPlotToMatlabFile(int computationNo, const QString& filename) const override;
    void exportIndexPlotToMatlabFile(int computationNo, const QString& filename) const override;

    void exportConvergenceLogsToMatlabFile(int computationNo, const QString& filename) const override;
    void exportOptionalLogToMatlabFile(int computationNo, int logIndex,
                                       const QString& filename) const override;
    void deleteSolution(int computationNo) override;

private:
    bool abortRequested_;

    std::unique_ptr<AbstractDiscretizedProblem> lastDiscretizedProblem_;
    std::map<int, std::unique_ptr<AbstractSolution> > solutions_;   // mapped by computationNumber
};

#endif // GUI_COMMUNICATOR_H
