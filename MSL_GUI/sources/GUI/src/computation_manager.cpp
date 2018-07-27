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


#include <QStatusBar>
#include <QProgressBar>
#include <QLabel>
#include <QTextEdit>

#include "computation_table.h"
#include "computation_status_box.h"
#include "interfaces/gui_inputdata.h"

#include "computation_manager.h"



ComputationManager::ComputationManager(QObject* parent, ComputationTable* table,
                                       ComputationStatusBox* statusBox,
                                       QStatusBar* statusBar,
                                       QTextEdit* textLogEdit)
    : QObject(parent),
      computationTable_(table),
      statusBox_(statusBox),
      statusBar_(statusBar),
      textLogEdit_(textLogEdit)
{
    progressBar_ = new QProgressBar;
    progressBar_->setMaximum(0);
    progressBar_->setMaximumWidth(490);
    progressBar_->setVisible(false);
    statusBar_->addWidget(progressBar_);

    connect(statusBox_, &ComputationStatusBox::abortRequested, this, &ComputationManager::abortRequested);
    connect(statusBox_, &ComputationStatusBox::coeffComputationTimeChanged, computationTable_, &ComputationTable::setTimeEntry);
    connect(computationTable_, &ComputationTable::selectedComputationChanged, this, &ComputationManager::handleSelectedComputationChanged);
}



void ComputationManager::startNewCoeffComputation(AbstractProblemTypeModule* problemType,
                                                  const GuiInputData& input,
                                                  const QString& summary)
{
    int computationNo = input.computationNumber;

    ComputationLogData logData;
    logData.input = input;
    logData.inputSummary = QString("<span style=\" font-style:italic; text-decoration: underline;\">Computation %1</span><br><br> ").arg(input.computationNumber) + summary;
    logData.coeffState = CoeffComputationState::RUNNING;
    logData.plotState = PlotComputationState::NOT_STARTED_YET;
    computationLogs_[computationNo] = logData;

    computationTable_->addEntry(computationNo);

    statusBox_->showCoefficientsComputing(computationNo);
    progressBar_->setVisible(true);

    emit coeffComputationRequested(problemType, input);
}



void ComputationManager::startNewPlotComputation(int computationNo, int resolution, bool subsequent)
{
    ComputationLogData& log = computationLogs_.at(computationNo);

    QString coeffEndState;

    switch (log.coeffState) {
    case CoeffComputationState::ERROR_COEFF_COMPUTATION:
        coeffEndState = QStringLiteral("Computing coefficients: Ended with error");
        break;
    case CoeffComputationState::ABORTED:
        coeffEndState = QStringLiteral("Computing coefficients: Aborted");
        break;
    case CoeffComputationState::COMPLETE:
        coeffEndState = QStringLiteral("Computing coefficients: Done ✔ ");
        break;
    default:
        coeffEndState = QStringLiteral("Computing coefficients: Unknown status");
        break;
    }

    statusBox_->showPlotComputing(computationNo, resolution, subsequent, coeffEndState);

    if (!subsequent)
    {
        progressBar_->setVisible(true);
    }

    log.plotState = PlotComputationState::RUNNING;

    if (getSelectedComputationNumber() == computationNo)
    {
        emit computationSelected(&log);
    }
    emit plotComputationRequested(computationNo, resolution);
}



void ComputationManager::deleteComputation(int computationNo)
{
    computationLogs_.erase(computationNo);
    computationTable_->removeEntry(computationNo);
}



int ComputationManager::getSelectedComputationNumber() const
{
    return computationTable_->getSelectedComputationNo();
}



QList<int> ComputationManager::getComputationNumbers(int& runningComputationNo) const
{
    QList<int> list;
    for (const auto& pair : computationLogs_)
    {
        list.append(pair.first);
        if ((pair.second.coeffState == CoeffComputationState::RUNNING)
             || (pair.second.plotState == PlotComputationState::RUNNING))
            runningComputationNo = pair.first;
    }
    return list;
}



void ComputationManager::handleEndOfCoeffComputation(int computationNo, CoeffComputationState::Enum endState)
{
    QString statusTableEntry;
    QString coeffEndState;

    switch (endState)
    {
    case CoeffComputationState::ERROR_PROBLEM_CREATION:
        statusTableEntry = coeffEndState = QStringLiteral("Error: problem setup failed");
        break;
    case CoeffComputationState::ERROR_NO_SOLUTION_CREATED:
        statusTableEntry = coeffEndState = QStringLiteral("Error: solution setup failed");
        break;
    case CoeffComputationState::ERROR_COEFF_COMPUTATION:
        statusTableEntry = QStringLiteral("Incomplete (error)");
        coeffEndState = QStringLiteral("Computing coefficients: Ended with error");
        break;
    case CoeffComputationState::ABORTED:
        statusTableEntry = QStringLiteral("Incomplete (aborted)");
        coeffEndState = QStringLiteral("Computing coefficients: Aborted");
        break;
    case CoeffComputationState::COMPLETE:
        statusTableEntry = QStringLiteral("Complete");
        coeffEndState = QStringLiteral("Computing coefficients: Done ✔ ");
        break;
    default:
        statusTableEntry = QStringLiteral("Unknown");
        coeffEndState = QStringLiteral("Computing coefficients: Unknown status");
        break;
    }
    statusBox_->showCoefficientsComputingEnded(coeffEndState);
    computationTable_->setStatusEntry(computationNo, statusTableEntry);
    //computationTable_->setTimeEntry(computationNo, statusBox_->getCoeffTimeString());

    ComputationLogData& log = computationLogs_.at(computationNo);

    log.coeffState = endState;
    log.textLog = textLogEdit_->toHtml();

    if (endState == CoeffComputationState::COMPLETE && log.input.plotSolutionRequested)
    {
        startNewPlotComputation(computationNo, log.input.resolution, true);
    }
    else
    {
        progressBar_->setVisible(false);
        statusBar_->showMessage(coeffEndState, 8000);
        emit computationEnded();
    }
}



void ComputationManager::handleEndOfPlotComputation(int computationNo, PlotComputationState::Enum endState, int resolution)
{
    ComputationLogData& log = computationLogs_.at(computationNo);
    log.plotState = endState;

    QString plotTableEntry;
    QString plotEndState;

    switch (endState)
    {
    case PlotComputationState::ERROR:
        plotTableEntry = QStringLiteral("Failed");
        plotEndState = QStringLiteral("Computing plot samples: Ended with error");
        break;
    case PlotComputationState::ABORTED:
        plotTableEntry = QStringLiteral("Aborted");
        plotEndState = QStringLiteral("Computing plot samples: Aborted");
        break;
    case PlotComputationState::COMPLETE:
        plotTableEntry = QString("Available (%1)").arg(resolution);
        plotEndState = QStringLiteral("Computing plot samples: Done ✔ ");
        break;
    default:
        plotTableEntry = QStringLiteral("Unknown");
        plotEndState = QStringLiteral("Computing plot samples: Unknown status");
        break;
    }

    computationTable_->setPlotEntry(computationNo, plotTableEntry);
    statusBox_->showPlotComputingEnded(plotEndState);

    progressBar_->setVisible(false);
    statusBar_->showMessage(plotEndState, 8000);
    emit computationEnded();
}



void ComputationManager::saveOptionalLogsList(int computationNo, const QStringList& availableOptionalLogs)
{
    computationLogs_.at(computationNo).optionalLogs = availableOptionalLogs;
}



void ComputationManager::saveComputedMatrixNorms(int computationNo, double norm_A, double norm_Ainv)
{
    computationLogs_.at(computationNo).inputSummary +=
            QString("<br><br>Matrix norms (computed):<br>||A|| = %3, ||A%1%2|| = %4").arg(QChar(0x207B), QChar(0x00B9)).arg(norm_A).arg(norm_Ainv);
}



void ComputationManager::handleSelectedComputationChanged(int computationNo)
{
    if (computationNo >= 0)
    {
        emit computationSelected(&computationLogs_.at(computationNo));
    }
    else
    {
        emit computationSelected(nullptr);
    }
}
