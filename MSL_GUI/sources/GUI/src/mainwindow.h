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


#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>

namespace Ui {
class MainWindow;
}
class ComputationManager;
struct ComputationLogData;
struct GuiInputData;
class MatrixNormDialog;
class InfoDefiningDialog;
class InfoControlsDialog;
class InfoAboutDialog;
class ProblemDefinitionManager;
class PlotManager;
class MslGuiPluginInterface;
class AbstractGuiCommunicator;
class AbstractParameterWidgetManager;
class AbstractProblemTypeModule;


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

protected:
    void closeEvent(QCloseEvent *event) override;

signals:

private slots:
    void handleSelectedProblemChanged(const QStringList& fullFunctionStrings, bool isExampleProblem);

    void handleNoProblemSelected();

    void handleEndOfComputation();

    void handleComputationSelected(const ComputationLogData* log);

    void adaptToSelectedComputation();

    void showErrorMessage(const QString& message);

    void showStatusMessage(const QString& message);

    void showDockWidgetInWindow(bool topLevel);

    void on_pushButton_start_clicked();

    void on_comboBox_domain_currentIndexChanged(const QString& domainName);

    void on_comboBox_problemClass_currentIndexChanged(int index);

    void on_comboBox_discretizationType_currentIndexChanged(int index);

    void on_comboBox_Basis1D_currentTextChanged(const QString& text);

    void on_comboBox_method_currentIndexChanged(int index);

    void on_pushButton_computeSolutionPlotSamples_clicked();

    void on_pushButton_export_clicked();

    void on_pushButton_deleteComputationEntry_clicked();

    void on_spinBox_jmax_valueChanged(int value);

    void on_checkBox_plotSolution_toggled(bool checked);

    void on_actionClear_computation_history_triggered();

    void on_actionDefining_custom_problems_triggered();

    void on_actionControls_triggered();

    void on_actionAbout_MSL_GUI_triggered();

private:
    bool loadPlugin();
    void setupAdditionalConnections();
    void setupDockWidgets();

    QString inputSummary();
    QString inputSummaryMatlab();
    void setGuiInput(GuiInputData& input);

    int getDefaultResolution(int jmax, const QString& basis1D) const;

    Ui::MainWindow* ui_;
    MatrixNormDialog* matrixNormDialog_;
    InfoDefiningDialog* infoDefiningDialog_;
    InfoControlsDialog* infoControlsDialog_;
    InfoAboutDialog* infoAboutDialog_;

    ComputationManager* computationManager_;
    ProblemDefinitionManager* problemDefinitionManager_;
    PlotManager* plotManager_;

    MslGuiPluginInterface* plugin_;
    AbstractGuiCommunicator* guiCommunicator_;

    AbstractParameterWidgetManager* parameterWidgetManager_;
    AbstractProblemTypeModule* selectedProblemType_;

    int computationNumber_;
    QThread workerThread_;

    bool computationIsRunning_;
    const ComputationLogData* selectedComputation_;
};

#endif // MAINWINDOW_H
