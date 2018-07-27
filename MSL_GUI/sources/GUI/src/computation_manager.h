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


#ifndef COMPUTATION_MANAGER_H
#define COMPUTATION_MANAGER_H

#include <QObject>

#include "interfaces/computation_state_enums.h"
#include "interfaces/gui_inputdata.h"

class ComputationTable;
class ComputationStatusBox;
class AbstractProblemTypeModule;
class QStatusBar;
class QLabel;
class QProgressBar;
class QTextEdit;


struct ComputationLogData
{
    GuiInputData input;
    QString inputSummary;
    QString textLog;
    QStringList optionalLogs;
    CoeffComputationState::Enum coeffState;
    PlotComputationState::Enum plotState;
};



class ComputationManager : public QObject
{
    Q_OBJECT
public:
    explicit ComputationManager(QObject* parent, ComputationTable* table,
                                ComputationStatusBox* statusBox, QStatusBar* statusBar, QTextEdit* textLogEdit);

    void startNewCoeffComputation(AbstractProblemTypeModule* problemType, const GuiInputData& input,
                                  const QString& summary);
    void startNewPlotComputation(int computationNo, int resolution, bool subsequent = false);

    void deleteComputation(int computationNo);
    int getSelectedComputationNumber() const;

    QList<int> getComputationNumbers(int& runningComputationNo) const;

signals:
    void coeffComputationRequested(AbstractProblemTypeModule* problemType, const GuiInputData& input) const;
    void abortRequested() const;
    void plotComputationRequested(int computationNo, int resolution) const;

    void computationSelected(const ComputationLogData* log) const;
    void computationEnded() const;

public slots:
    void handleEndOfCoeffComputation(int computationNo, CoeffComputationState::Enum endState);
    void handleEndOfPlotComputation(int computationNo, PlotComputationState::Enum endState, int resolution);

    void saveOptionalLogsList(int computationNo, const QStringList& availableOptionalLogs);
    void saveComputedMatrixNorms(int computationNo, double norm_A, double norm_Ainv);

private slots:
    void handleSelectedComputationChanged(int computationNo);

private:
    ComputationTable* computationTable_;
    ComputationStatusBox* statusBox_;

    QStatusBar* statusBar_;
    QProgressBar* progressBar_;
    QTextEdit* textLogEdit_;

    // maps computation numbers to computation logs:
    std::map<int, ComputationLogData> computationLogs_;
};

#endif // COMPUTATION_MANAGER_H
