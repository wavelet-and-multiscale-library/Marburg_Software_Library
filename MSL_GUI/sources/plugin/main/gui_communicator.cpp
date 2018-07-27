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


#include "gui_communicator.h"
#include "GUI/interfaces/abstract_problemtype_module.h"
#include "GUI/interfaces/gui_inputdata.h"
#include "misc/string_conversion.h"



GuiCommunicator::GuiCommunicator() :
    abortRequested_(false),
    lastDiscretizedProblem_()
{

}



void GuiCommunicator::sendPlotDataToGui(const MathTL::SampledMapping<1>& plotData, int computationNo, bool newPlot) const
{
    const double* xValues = plotData.points().begin();
    const double* yValues = plotData.values().begin();
    int sampleCount = int (plotData.size());

    emit dataFor1DSolutionPlotComputed(computationNo, xValues, yValues, sampleCount, newPlot);
}



void GuiCommunicator::sendPlotDataToGui(const MathTL::SampledMapping<2>& plotData, int computationNo, bool newPlot) const
{
    const double* xValues = plotData.gridx().entries_vector().begin();
    const double* yValues = plotData.gridy().entries_vector().begin();
    const double* zValues = plotData.values().entries_vector().begin();

    int sampleCountX = int (plotData.gridx().column_dimension());
    int sampleCountY = int (plotData.gridx().row_dimension());

    emit dataFor2DSolutionPlotComputed(computationNo, xValues, yValues, zValues, sampleCountX, sampleCountY, newPlot);

}



void GuiCommunicator::sendPlotDataToGui(const MathTL::Array1D< MathTL::SampledMapping<1> >& plotData, int computationNo) const
{
    auto begin = plotData.begin();
    for (auto it = plotData.begin(); it != plotData.end(); ++it)
    {
        if (it == begin)
            sendPlotDataToGui(*it, computationNo, true);
        else
            sendPlotDataToGui(*it, computationNo, false);
    }
}



void GuiCommunicator::sendPlotDataToGui(const MathTL::Array1D< MathTL::SampledMapping<2> >& plotData, int computationNo) const
{
    auto begin = plotData.begin();
    for (auto it = plotData.begin(); it != plotData.end(); ++it)
    {
        if (it == begin)
            sendPlotDataToGui(*it, computationNo, true);
        else
            sendPlotDataToGui(*it, computationNo, false);
    }
}



bool GuiCommunicator::gotAbortRequest() const
{
    return abortRequested_;
}



void GuiCommunicator::computeSolutionCoeffs(AbstractProblemTypeModule* problemType,
                                            const GuiInputData& input)
{
    abortRequested_ = false;

    bool canReuseLastProblem = input.reuse_if_possible &&
                               (lastDiscretizedProblem_ != nullptr) &&
                               lastDiscretizedProblem_->canBeReusedFor(input);

    if (canReuseLastProblem) {
        try
        {
            lastDiscretizedProblem_->updateTo(input);
            emit statusMessageGenerated("Reusing last discretized problem.");
        }
        catch(const std::exception& theException)
        {
            QString errorMessage("Error while trying to update last "
                                 "discretized problem to new input:\n");
            errorMessage.append(QString(theException.what()));
            emit statusMessageGenerated(errorMessage);
            lastDiscretizedProblem_.reset();
            canReuseLastProblem = false;
        }
        catch(...)
        {
            QString errorMessage("Unknown error while trying to update last "
                                 "discretized problem to new input.");
            emit statusMessageGenerated(errorMessage);
            lastDiscretizedProblem_.reset();
            canReuseLastProblem = false;
        }
    }

    if (!canReuseLastProblem) {

        AbstractDiscretizedProblem* newProblem;

        try
        {
            newProblem = problemType->createDiscretizedProblem(input);
        }
        catch(const std::exception& theException)
        {
            QString errorMessage("Error during creation of discretized problem:\n\n");
            errorMessage.append(QString(theException.what()));
            emit errorOccured(errorMessage);
            emit coeffComputationEnded(input.computationNumber, CoeffComputationState::ERROR_PROBLEM_CREATION);
            return;
        }
        catch(...)
        {
            QString errorMessage("Unknown error during creation of discretized problem.");
            emit errorOccured(errorMessage);
            emit coeffComputationEnded(input.computationNumber, CoeffComputationState::ERROR_PROBLEM_CREATION);
            return;
        }

        lastDiscretizedProblem_.reset(newProblem);
    }

    AbstractSolution* newSolution;
    try
    {
        newSolution = lastDiscretizedProblem_->computeSolution(input.method, input.epsilon, input.jmax, this);
    }
    catch(...)
    {
        QString errorMessage("Error: Could not create solution object!");
        emit errorOccured(errorMessage);
        emit coeffComputationEnded(input.computationNumber, CoeffComputationState::ERROR_NO_SOLUTION_CREATED);
        return;
    }

    solutions_[input.computationNumber] = std::unique_ptr<AbstractSolution>(newSolution);
    emit solutionSaved(input.computationNumber, newSolution->getOptionalLogNames());

    CoeffComputationState::Enum endState;

    if (newSolution->isIncompleteDueToAbort())
    {
        endState = CoeffComputationState::ABORTED;
    }
    else
    {
        if (newSolution->isIncompleteDueToError())
        {
            endState = CoeffComputationState::ERROR_COEFF_COMPUTATION;
        }
        else
        {
            endState = CoeffComputationState::COMPLETE;
            if (!input.normEstimatesProvided)
            {
                double norm_A = lastDiscretizedProblem_->norm_A();
                double norm_Ainv = lastDiscretizedProblem_->norm_Ainv();

                emit matrixNormsComputed(input.computationNumber, norm_A, norm_Ainv);
            }
        }
    }

    emit coeffComputationEnded(input.computationNumber, endState);
}



void GuiCommunicator::requestAbort()
{
    abortRequested_ = true;
}



void GuiCommunicator::computeSolutionPlot(int computationNo, int resolution)
{
    abortRequested_ = false;
    try
    {
        solutions_.at(computationNo)->computePlotData(resolution);
    }
    catch(const std::exception& theException)
    {
        QString errorMessage("Error during computation of plot data:\n\n");
        errorMessage.append(QString(theException.what()));
        emit errorOccured(errorMessage);
        emit plotComputationEnded(computationNo, PlotComputationState::ERROR, -1);
        return;
    }
    catch(...)
    {
        QString errorMessage("Unknown error during computation of plot data.");
        emit errorOccured(errorMessage);
        emit plotComputationEnded(computationNo, PlotComputationState::ERROR, -1);
        return;
    }

    if (gotAbortRequest()) {
        emit plotComputationEnded(computationNo, PlotComputationState::ABORTED, -1);
    }
    else {
        solutions_.at(computationNo)->sendPlotDataToGuiVia(*this);
        emit plotComputationEnded(computationNo, PlotComputationState::COMPLETE, resolution);
    }
}



void GuiCommunicator::exportSolutionPlotToMatlabFile(int computationNo, const QString& filename) const
{
    try
    {
        solutions_.at(computationNo)->writeSolutionPlotToMatlabFile(cStringFromQString(filename));
    }
    catch(const std::exception& theException)
    {
        QString errorMessage = QString("Error during export of plot data "
                                       "to Matlab file \'%1\':\n\n").arg(filename);
        errorMessage.append(QString(theException.what()));
        emit errorOccured(errorMessage);
        return;
    }
    catch(...)
    {
        QString errorMessage = QString("Unknown error during export of plot data "
                                       "to Matlab file \'%1\'.").arg(filename);
        emit errorOccured(errorMessage);
        return;
    }

    QString message = QString("Plot data successfully exported to file \'%1\'.").arg(filename);
    emit statusMessageGenerated(message);
}



void GuiCommunicator::exportIndexPlotToMatlabFile(int computationNo, const QString& filename) const
{
    try
    {
        solutions_.at(computationNo)->writeIndexPlotToMatlabFile(cStringFromQString(filename));
    }
    catch(const std::exception& theException)
    {
        QString errorMessage = QString("Error during export of index plot data "
                                       "to Matlab file \'%1\':\n\n").arg(filename);
        errorMessage.append(QString(theException.what()));
        emit errorOccured(errorMessage);
        return;
    }
    catch(...)
    {
        QString errorMessage = QString("Unknown error during export of index plot data "
                                       "to Matlab file \'%1\'.").arg(filename);
        emit errorOccured(errorMessage);
        return;
    }

    QString message = QString("Index plot data successfully exported "
                              "to file \'%1\'.").arg(filename);
    emit statusMessageGenerated(message);
}



void GuiCommunicator::exportConvergenceLogsToMatlabFile(int computationNo, const QString& filename) const
{
    try
    {
        solutions_.at(computationNo)->writeConvergenceLogsToMatlabFile(cStringFromQString(filename));
    }
    catch(const std::exception& theException)
    {
        QString errorMessage = QString("Error during export of convergence log data "
                                       "to Matlab file \'%1\':\n\n").arg(filename);
        errorMessage.append(QString(theException.what()));
        emit errorOccured(errorMessage);
        return;
    }
    catch(...)
    {
        QString errorMessage = QString("Unknown error during export of convergence log data "
                                       "to Matlab file \'%1\'.").arg(filename);
        emit errorOccured(errorMessage);
        return;
    }

    QString message = QString("Convergence log data successfully exported "
                              "to file \'%1\'.").arg(filename);
    emit statusMessageGenerated(message);
}



void GuiCommunicator::exportOptionalLogToMatlabFile(int computationNo, int logIndex,
                                                    const QString& filename) const
{
    try
    {
        solutions_.at(computationNo)->writeOptionalLogToMatlabFile(logIndex, cStringFromQString(filename));
    }
    catch(const std::exception& theException)
    {
        QString errorMessage = QString("Error during export of optional log "
                                       "to Matlab file \'%1\':\n\n").arg(filename);
        errorMessage.append(QString(theException.what()));
        emit errorOccured(errorMessage);
        return;
    }
    catch(...)
    {
        QString errorMessage = QString("Unknown error during export of optional log "
                                       "to Matlab file \'%1\'.").arg(filename);
        emit errorOccured(errorMessage);
        return;
    }

    QString message = QString("Optional log successfully exported "
                              "to file \'%1\'.").arg(filename);
    emit statusMessageGenerated(message);
}



void GuiCommunicator::deleteSolution(int computationNo)
{
    solutions_.erase(computationNo);
}
