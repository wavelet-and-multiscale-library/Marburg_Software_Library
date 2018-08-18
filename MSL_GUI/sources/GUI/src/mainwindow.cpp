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


#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <iostream>

#include <QtWidgets>

#include "interfaces/mslgui_plugin_interface.h"
#include "interfaces/abstract_gui_communicator.h"
#include "interfaces/abstract_problemtype_module.h"
#include "interfaces/abstract_parameterwidget_manager.h"
#include "computation_manager.h"
#include "matrixnorm_dialog.h"
#include "info_defining_dialog.h"
#include "info_controls_dialog.h"
#include "info_about_dialog.h"
#include "problem_definition_manager.h"
#include "plot_manager.h"
#include "qoutput_stream.h"



MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui_(new Ui::MainWindow),
    infoDefiningDialog_(nullptr),
    infoControlsDialog_(nullptr),
    infoAboutDialog_(nullptr),
    computationNumber_(1),
    computationIsRunning_(false),
    selectedComputation_(nullptr)
{
    if (!loadPlugin())
    {
        throw false;
    }

    ui_->setupUi(this);
    setupDockWidgets();
    resize(1280, 500);
    setWindowIcon(QIcon(":/icons/Logo-WMSL.png"));

    // Redirecting console output to textEdit_log:
    new QOutputStream(std::cout, ui_->textEdit_runningComputationTextLog, this);

    ui_->computationStatusBox->setWidgets(ui_->widget_coeffStatus, ui_->widget_plotStatus);
    ui_->computationStatusBox->setAbortButton(ui_->pushButton_abort);
    ui_->computationStatusBox->setLabels(ui_->label_coeffStatus, ui_->label_coeffTime,
                                         ui_->label_plotStatus, ui_->label_plotTime);

    guiCommunicator_ = &(plugin_->getGuiCommunicator());
    matrixNormDialog_ = new MatrixNormDialog(this);
    computationManager_ = new ComputationManager(this, ui_->computationTable, ui_->computationStatusBox, statusBar(), ui_->textEdit_runningComputationTextLog);
    problemDefinitionManager_ = new ProblemDefinitionManager(this, ui_->comboBox_problem);
    plotManager_ = new PlotManager(this, ui_->convergenceView, ui_->convergenceTimeView, ui_->stackedWidget_solutionPlot, ui_->label_samplePointsCount);
    ui_->convergenceView->setRubberBand(QChartView::RectangleRubberBand);
    ui_->convergenceView->setRenderHint(QPainter::Antialiasing);
    ui_->convergenceTimeView->setRubberBand(QChartView::RectangleRubberBand);
    ui_->convergenceTimeView->setRenderHint(QPainter::Antialiasing);


    parameterWidgetManager_ = &(plugin_->getParameterWidgetManager());
    parameterWidgetManager_->addAllWidgetsTo(ui_->stackedWidget_methodParameters);

    qRegisterMetaType<CoeffComputationState::Enum>("CoeffComputationState::Enum");
    qRegisterMetaType<PlotComputationState::Enum>("PlotComputationState::Enum");
    qRegisterMetaType<GuiInputData>("GuiInputData");
    setupAdditionalConnections();

    plugin_->initProblemTypeModules();
    ui_->comboBox_domain->addItems(plugin_->getDomainNameList());
    ui_->comboBox_domain->setCurrentIndex(0);

    adaptToSelectedComputation();
}



MainWindow::~MainWindow()
{
    workerThread_.quit();
//    guiCommunicator_->requestAbort();
//    workerThread_.wait();
    delete ui_;
    delete plugin_;
}



void MainWindow::closeEvent(QCloseEvent* event)
{
    if (!computationIsRunning_)
        event->accept();
    else
    {
        int answer = QMessageBox::question(this, QStringLiteral("Stop computation and quit?"),
                              QStringLiteral("There is still a computation running. Are you sure you want "
                                             "to stop it and quit?"), QMessageBox::Yes, QMessageBox::No);
        if (answer == QMessageBox::No)
        {
            event->ignore();
        }
        else
        {
            event->accept();
            hide();
            workerThread_.quit();
            guiCommunicator_->requestAbort();
            if (!workerThread_.wait(10000))
                workerThread_.terminate();
        }

    }
}



bool MainWindow::loadPlugin()
{
    QDir pluginDir(qApp->applicationDirPath());

    QString pluginName = "mslgui_plugin";

    QPluginLoader pluginLoader(pluginDir.absoluteFilePath(pluginName));
    QObject* pluginObject = pluginLoader.instance();

    QString loadErrors;

    if(pluginObject) {
        plugin_ = qobject_cast<MslGuiPluginInterface*>(pluginObject);
        if (plugin_) {
            return true;
        }
        else {
            loadErrors = QString("The shared library \'%1\' does not implement a compatible "
                                 "version of MslGuiPluginInterface!").arg(pluginName);
        }
    }
    else {
        loadErrors = QString("Could not load the shared library \'%1\'! Error message:\n\n%2")
                            .arg(pluginName, pluginLoader.errorString());
    }

    showErrorMessage(loadErrors);
    return false;
}



void MainWindow::setupAdditionalConnections()
{
    bool success =

       connect(computationManager_, &ComputationManager::abortRequested,
               guiCommunicator_, &AbstractGuiCommunicator::requestAbort, Qt::DirectConnection)

    && connect(computationManager_, &ComputationManager::computationSelected,
               this, &MainWindow::handleComputationSelected)

    && connect(computationManager_, &ComputationManager::coeffComputationRequested,
               guiCommunicator_, &AbstractGuiCommunicator::computeSolutionCoeffs)

    && connect(computationManager_, &ComputationManager::plotComputationRequested,
               guiCommunicator_, &AbstractGuiCommunicator::computeSolutionPlot)

    && connect(guiCommunicator_, SIGNAL(coeffComputationEnded(int,CoeffComputationState::Enum)),
               computationManager_, SLOT(handleEndOfCoeffComputation(int,CoeffComputationState::Enum)))

    && connect(guiCommunicator_, SIGNAL(plotComputationEnded(int,PlotComputationState::Enum,int)),
               computationManager_, SLOT(handleEndOfPlotComputation(int,PlotComputationState::Enum,int)))

    && connect(guiCommunicator_, SIGNAL(solutionSaved(int,QStringList)),
               computationManager_, SLOT(saveOptionalLogsList(int,QStringList)))

    && connect(guiCommunicator_, SIGNAL(matrixNormsComputed(int,double,double)),
               computationManager_, SLOT(saveComputedMatrixNorms(int,double,double)))

    && connect(computationManager_, &ComputationManager::computationEnded,
               this, &MainWindow::handleEndOfComputation)

    && connect(guiCommunicator_, SIGNAL(convergencePlotDataPairComputed(double,double)),
               plotManager_, SLOT(addToCurrentConvPlot(double,double)))

    && connect(guiCommunicator_, SIGNAL(convergenceTimePlotDataPairComputed(double,double)),
               plotManager_, SLOT(addToCurrentConvTimePlot(double,double)))

    && connect(ui_->computationTable, &ComputationTable::checkboxToogled,
               plotManager_, &PlotManager::toggleConvPlotsVisibility)

    && connect(guiCommunicator_, SIGNAL(errorOccured(const QString&)),
            this, SLOT(showErrorMessage(const QString&)))

    && connect(guiCommunicator_, SIGNAL(statusMessageGenerated(const QString&)),
            this, SLOT(showStatusMessage(const QString&)))

    && connect(guiCommunicator_, SIGNAL(dataFor1DSolutionPlotComputed(int,const double*,const double*,int,bool)),
               plotManager_, SLOT(setup1DsolutionPlotPatch(int,const double*,const double*,int,bool)))

    && connect(guiCommunicator_, SIGNAL(dataFor2DSolutionPlotComputed(int,const double*,const double*,const double*,int,int,bool)),
               plotManager_, SLOT(setup2DsolutionPlotPatch(int,const double*,const double*,const double*,int,int,bool)))

    && connect(ui_->pushButton_addProblem, &QPushButton::clicked,
               problemDefinitionManager_, &ProblemDefinitionManager::handleAddNewProblemDefinitionClicked)

    && connect(ui_->pushButton_deleteProblem, &QPushButton::clicked,
               problemDefinitionManager_, &ProblemDefinitionManager::handleDeleteProblemDefinitionClicked)

    && connect(ui_->pushButton_editProblem, &QPushButton::clicked,
               problemDefinitionManager_, &ProblemDefinitionManager::handleEditProblemDefinitionClicked)

    && connect(problemDefinitionManager_, &ProblemDefinitionManager::selectedProblemChanged,
               this, &MainWindow::handleSelectedProblemChanged)

    && connect(problemDefinitionManager_, &ProblemDefinitionManager::noProblemSelected,
               this, &MainWindow::handleNoProblemSelected)

    && connect(ui_->actionHide_all_convergence_plots, &QAction::triggered,
               ui_->computationTable, &ComputationTable::uncheckAllCheckboxes)

    && connect(ui_->actionShow_all_convergence_plots, &QAction::triggered,
               ui_->computationTable, &ComputationTable::checkAllCheckboxes);


    if (!success)
    {
        showErrorMessage("Some signal/slot connections failed! Essential functions may not work properly.");
    }

    guiCommunicator_->moveToThread(&workerThread_);
    workerThread_.start();
}



void MainWindow::setupDockWidgets()
{
    tabifyDockWidget(ui_->dockWidget_convergence, ui_->dockWidget_tempConvergence);

    ui_->dockWidget_convergence->raise();

    connect(ui_->dockWidget_solutionPlot, SIGNAL(topLevelChanged(bool)), this, SLOT(showDockWidgetInWindow(bool)));
    connect(ui_->dockWidget_convergence, SIGNAL(topLevelChanged(bool)), this, SLOT(showDockWidgetInWindow(bool)));
    connect(ui_->dockWidget_tempConvergence, SIGNAL(topLevelChanged(bool)), this, SLOT(showDockWidgetInWindow(bool)));
}



QString MainWindow::inputSummary()
{
    QString jmaxString = QStringLiteral("j%1%2%3 = %4").arg(QChar(0x2098), QChar(0x2090), QChar(0x2093))
                   .arg(ui_->spinBox_jmax->value());
    if (ui_->spinBox_pmax->isEnabled())
        jmaxString += QStringLiteral(", p%1%2%3 = %4").arg(QChar(0x2098), QChar(0x2090), QChar(0x2093))
                .arg(ui_->spinBox_pmax->value());

    if (ui_->comboBox_discretizationType->currentText().contains(QStringLiteral("aggregated"), Qt::CaseInsensitive))
    {
        if (ui_->comboBox_domain->currentText().contains(QStringLiteral("ring"), Qt::CaseInsensitive))
            jmaxString.append(QStringLiteral(", patch overlap: 1.0"));
        else
            jmaxString.append(QStringLiteral(", patch overlap: %1").arg(ui_->doubleSpinBox_overlap->value()));
    }

    QString parameterString = parameterWidgetManager_->getParamStringForMethod(ui_->comboBox_method->currentText());
    if (!parameterString.isEmpty())
        parameterString.append("<br>");

    QString epsilonString;
    int factor = ui_->spinBox_epsilonFactor->value();
    QString text = ui_->spinBox_epsilonExponent->text();

    if (factor == 1)
    {
        epsilonString = QString("ε = %1").arg(text);
    }
    else
    {
        epsilonString = QString("ε = %1%2%3").arg(factor).arg(QChar(0x22C5), text);
    }

    QString summary = QStringLiteral("%1 on %2<br>"
                                     "<span style=\" font-style:italic;\">%3</span>"
                                     "<span style=\" font-family:'STIX';\">%4</span>"
                                     "<br>"
                                     "%5<br>"
                                     "%6<br>"
                                     "%7<br>"
                                     "<br>"
                                     "%8<br>"
                                     "%9"
                                     "%10");

    summary = summary.arg(ui_->comboBox_problemClass->currentText(),
                          ui_->comboBox_domain->currentText(),
                          ui_->comboBox_problem->currentText(),
                          ui_->textEdit_definitions->toHtml(),
                          ui_->comboBox_discretizationType->currentText(),
                          ui_->comboBox_Basis1D->currentText(),
                          jmaxString,
                          ui_->comboBox_method->currentText(),
                          parameterString)
                     .arg(epsilonString);

    return summary;
}



QString MainWindow::inputSummaryMatlab()
{
    QString jmaxString = QStringLiteral("j_{max} = %1").arg(ui_->spinBox_jmax->value());
    if (ui_->spinBox_pmax->isEnabled())
        jmaxString += QStringLiteral(", p_{max} = %1").arg(ui_->spinBox_pmax->value());

    if (ui_->comboBox_discretizationType->currentText().contains(QStringLiteral("aggregated"), Qt::CaseInsensitive))
    {
        if (ui_->comboBox_domain->currentText().contains(QStringLiteral("ring"), Qt::CaseInsensitive))
            jmaxString.append(QStringLiteral(", patch overlap: 1.0"));
        else
            jmaxString.append(QStringLiteral(", patch overlap: %1").arg(ui_->doubleSpinBox_overlap->value()));
    }

    QString epsilonString;
    int factor = ui_->spinBox_epsilonFactor->value();

    if (factor == 1)
    {
        epsilonString = QString("10^{-%1}").arg(ui_->spinBox_epsilonExponent->value());
    }
    else
    {
        epsilonString = QString("%1*10^{-%2}").arg(factor).arg(ui_->spinBox_epsilonExponent->value());
    }

    QString parameterString = parameterWidgetManager_->getParamStringForMethod(ui_->comboBox_method->currentText());

    QString summary = QStringLiteral("{'%1 on %2 (%3)', '%4', '%5', '%6, %7', '%8, \\epsilon = %9', '%10'}");

    summary = summary.arg(ui_->comboBox_problemClass->currentText(),
                          ui_->comboBox_domain->currentText(),
                          ui_->comboBox_problem->currentText(),
                          ui_->textEdit_definitions->toPlainText(),
                          ui_->comboBox_discretizationType->currentText(),
                          ui_->comboBox_Basis1D->currentText(),
                          jmaxString,
                          ui_->comboBox_method->currentText(),
                          epsilonString)
                     .arg(parameterString);

    return summary;
}



void MainWindow::setGuiInput(GuiInputData& input)
{
    input.computationNumber = computationNumber_;

    input.domain = ui_->comboBox_domain->currentText();
    input.problemType = ui_->comboBox_problemClass->currentText();
    input.problemName = ui_->comboBox_problem->currentText();
    input.discretizationType = ui_->comboBox_discretizationType->currentText();
    input.basis1D = ui_->comboBox_Basis1D->currentText();
    input.method = ui_->comboBox_method->currentText();

    input.exampleIndex = problemDefinitionManager_->currentExampleIndex();
    input.discretizationTypeIndex = ui_->comboBox_discretizationType->currentIndex();

    input.functionDefinitions = problemDefinitionManager_->getFuncDefStringsForSelectedProblem();

    input.jmax = ui_->spinBox_jmax->value();
    if (ui_->spinBox_pmax->isEnabled())
        input.pmax = ui_->spinBox_pmax->value();
    else
        input.pmax = -1;

    if (ui_->doubleSpinBox_overlap->isEnabled())
        input.overlap = ui_->doubleSpinBox_overlap->value();
    else
        input.overlap = 1.0;

    input.epsilon = double (ui_->spinBox_epsilonFactor->value()) * pow(10.0, -(ui_->spinBox_epsilonExponent->value()));

    input.plotSolutionRequested = ui_->checkBox_plotSolution->isChecked();
    input.resolution = ui_->spinBox_resolution->value();

    input.normEstimatesProvided = ui_->checkBox_normEstimates->isChecked();
    input.reuse_if_possible = ui_->checkBox_reuse->isChecked();

    input.matlabSummary = inputSummaryMatlab();
}



void MainWindow::handleSelectedProblemChanged(const QStringList& fullFunctionStrings, bool isExampleProblem)
{
    if (isExampleProblem)
    {
        ui_->pushButton_editProblem->setEnabled(false);
        ui_->pushButton_deleteProblem->setEnabled(false);
    }
    else
    {
        ui_->pushButton_editProblem->setEnabled(true);
        ui_->pushButton_deleteProblem->setEnabled(true);
    }

    ui_->textEdit_definitions->clear();
    for (const QString& funcDef : fullFunctionStrings)
        ui_->textEdit_definitions->append(funcDef);
}



void MainWindow::handleNoProblemSelected()
{
    ui_->pushButton_editProblem->setEnabled(false);
    ui_->pushButton_deleteProblem->setEnabled(false);
    ui_->textEdit_definitions->clear();
}



void MainWindow::handleEndOfComputation()
{
    ui_->pushButton_start->setEnabled(true);
    computationIsRunning_ = false;
    adaptToSelectedComputation();
}



void MainWindow::handleComputationSelected(const ComputationLogData* log)
{
    selectedComputation_ = log;
    adaptToSelectedComputation();
}



void MainWindow::adaptToSelectedComputation()
{
    if (selectedComputation_)
    {
        if (selectedComputation_->coeffState == CoeffComputationState::RUNNING ||
            selectedComputation_->plotState == PlotComputationState::RUNNING)
        {
            ui_->pushButton_deleteComputationEntry->setEnabled(false);
        }
        else
        {
            ui_->pushButton_deleteComputationEntry->setEnabled(true);
        }

        if (computationIsRunning_ || selectedComputation_->coeffState == CoeffComputationState::ERROR_PROBLEM_CREATION
            || selectedComputation_->coeffState == CoeffComputationState::ERROR_NO_SOLUTION_CREATED)
        {
            ui_->pushButton_computeSolutionPlotSamples->setEnabled(false);
        }
        else
        {
            ui_->pushButton_computeSolutionPlotSamples->setEnabled(true);
        }

        ui_->comboBox_fileExport->setItems(selectedComputation_->coeffState, selectedComputation_->plotState,
                                           selectedComputation_->optionalLogs);
        if (ui_->comboBox_fileExport->count() == 0)
        {
            ui_->pushButton_export->setEnabled(false);
        }
        else
        {
            ui_->pushButton_export->setEnabled(true);
        }

        ui_->textEdit_details->setHtml(selectedComputation_->inputSummary);
        ui_->label_textLog->setText(QString("Text Log (Computation %1)").arg(selectedComputation_->input.computationNumber));

        if (selectedComputation_->coeffState == CoeffComputationState::RUNNING)
        {
            ui_->stackedWidget_textLog->setCurrentWidget(ui_->page_runningComputationTextLog);
        }
        else
        {
            ui_->textEdit_textLog->setHtml(selectedComputation_->textLog);
            ui_->stackedWidget_textLog->setCurrentWidget(ui_->page_textLog);
        }


        plotManager_->showSolutionPlot(selectedComputation_->input.computationNumber);

        if (selectedComputation_->plotState == PlotComputationState::COMPLETE)
            ui_->dockWidget_solutionPlot->setWindowTitle(QStringLiteral("Solution plot (Computation %1)").arg(selectedComputation_->input.computationNumber));
        else
            ui_->dockWidget_solutionPlot->setWindowTitle(QStringLiteral("Marburg Software Library"));
    }
    else
    {
        ui_->textEdit_details->clear();
        ui_->label_textLog->setText("Text Log");
        ui_->textEdit_textLog->clear();
        ui_->stackedWidget_textLog->setCurrentWidget(ui_->page_textLog);
        ui_->comboBox_fileExport->clear();
        ui_->pushButton_export->setEnabled(false);

        ui_->pushButton_computeSolutionPlotSamples->setEnabled(false);
        ui_->pushButton_deleteComputationEntry->setEnabled(false);

        ui_->stackedWidget_solutionPlot->setCurrentIndex(0);
        ui_->dockWidget_solutionPlot->setWindowTitle(QStringLiteral("Marburg Software Library"));
        ui_->label_samplePointsCount->setText("");
    }


}



void MainWindow::showErrorMessage(const QString& message)
{
    QMessageBox::critical(this, "Error", message);
}



void MainWindow::showStatusMessage(const QString& message)
{
    ui_->statusBar->showMessage(message, 8000);
}



void MainWindow::showDockWidgetInWindow(bool topLevel)
{
    if (topLevel)
    {
        QWidget* widget = qobject_cast<QWidget*>(sender());
        widget->setWindowFlags(Qt::Window | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint);
        widget->show();
    }

}



void MainWindow::on_pushButton_start_clicked()
{
    if (ui_->comboBox_problem->currentIndex() < 0)
    {
        QMessageBox::warning(this, QStringLiteral("No problem specified!"),
                             QStringLiteral("There is no problem definition specified.\n"
                                            "Click on the \'+\'-button to add a new custom problem definition."),
                             QMessageBox::Ok);
        return;
    }

    GuiInputData input;
    setGuiInput(input);

    QString summary = inputSummary();
    QTextDocument doc;
    doc.setHtml(summary);

    if (input.normEstimatesProvided)
    {
        matrixNormDialog_->resetValuesToOne();
        if (matrixNormDialog_->exec() == QDialog::Accepted)
        {
            input.normA = matrixNormDialog_->getNormA();
            input.normAinv = matrixNormDialog_->getNormAinv();
            summary += QString("<br><br>Matrix norms (user provided):<br>||A|| = %3, ||A%1%2|| = %4").arg(QChar(0x207B), QChar(0x00B9)).arg(input.normA).arg(input.normAinv);
        }
        else
        {
            return;
        }
    }

    ui_->pushButton_start->setEnabled(false);
    ui_->pushButton_computeSolutionPlotSamples->setEnabled(false);
    computationIsRunning_ = true;
    ui_->textEdit_runningComputationTextLog->clear();
    ui_->textEdit_runningComputationTextLog->append(QString("<span style=\" font-style:italic; text-decoration: underline;\">Computation %1</span><br>").arg(computationNumber_));


    plotManager_->addNewConvergencePlots(computationNumber_, doc.toPlainText());
    computationManager_->startNewCoeffComputation(selectedProblemType_, input, summary);

    ui_->tabWidget->setCurrentIndex(1);
    computationNumber_++;
}



void MainWindow::on_comboBox_domain_currentIndexChanged(const QString& domainName)
{
    ui_->comboBox_problemClass->clear();
    ui_->comboBox_problemClass->addItems(plugin_->getProblemTypeListForDomain(domainName));
    ui_->comboBox_problemClass->setCurrentIndex(0);
}



void MainWindow::on_comboBox_problemClass_currentIndexChanged(int index)
{
    if (index >= 0)
    {
        selectedProblemType_ = plugin_->getProblemTypeModule(ui_->comboBox_domain->currentText(), index);
        problemDefinitionManager_->handleProblemTypeSelection(selectedProblemType_);
        ui_->textEdit_descriptionProblemtype->setPlainText(selectedProblemType_->getProblemTypeDescription());

        ui_->comboBox_discretizationType->clear();
        ui_->comboBox_discretizationType->addItems(selectedProblemType_->getDiscretizationTypeList());
        ui_->comboBox_discretizationType->setCurrentIndex(0);
    }
}



void MainWindow::on_comboBox_discretizationType_currentIndexChanged(int index)
{
    if (index >= 0)
    {
        ui_->comboBox_Basis1D->clear();
        ui_->comboBox_Basis1D->addItems(selectedProblemType_->get1DBasisList(index));
        ui_->comboBox_Basis1D->setCurrentIndex(0);

        ui_->comboBox_method->clear();
        ui_->comboBox_method->addItems(selectedProblemType_->getMethodList(index));
        ui_->comboBox_method->setCurrentIndex(0);

        int jmax = int (selectedProblemType_->getJmaxStandard(index));
        ui_->spinBox_jmax->setValue(jmax);
        ui_->spinBox_resolution->setValue(jmax);

        if (ui_->comboBox_discretizationType->currentText().contains(QStringLiteral("quark"), Qt::CaseInsensitive))
        {
            ui_->stackedWidget_pmax_overlap->setCurrentWidget(ui_->page_pmax);
            ui_->spinBox_pmax->setEnabled(true);
        }
        else
        {
            ui_->spinBox_pmax->setEnabled(false);

            if (ui_->comboBox_discretizationType->currentText().contains(QStringLiteral("aggregated"), Qt::CaseInsensitive))
            {
                ui_->stackedWidget_pmax_overlap->setCurrentWidget(ui_->page_overlap);
                if (ui_->comboBox_domain->currentText().contains(QStringLiteral("ring"), Qt::CaseInsensitive))
                {
                    ui_->doubleSpinBox_overlap->setValue(1.0);
                    ui_->doubleSpinBox_overlap->setEnabled(false);
                }
                else
                {
                    ui_->doubleSpinBox_overlap->setEnabled(true);
                }
            }
            else
            {
                ui_->doubleSpinBox_overlap->setEnabled(false);
                ui_->stackedWidget_pmax_overlap->setCurrentWidget(ui_->page_void);
            }
        }
    }
}



void MainWindow::on_comboBox_method_currentIndexChanged(int index)
{
    if (index >= 0)
    {
        int discrIndex = ui_->comboBox_discretizationType->currentIndex();
        ui_->spinBox_epsilonExponent->setValue(selectedProblemType_->getEpsilonStandardExponent(discrIndex, index));
        ui_->spinBox_epsilonFactor->setValue(1);

        QString methodName = ui_->comboBox_method->currentText();
        parameterWidgetManager_->displayWidgetForMethod(methodName);
    }
}



void MainWindow::on_pushButton_computeSolutionPlotSamples_clicked()
{
    int defaultResolution = selectedComputation_->input.jmax;

    if (selectedComputation_->input.basis1D.contains(QStringLiteral("DS")))
        defaultResolution++;

    bool ok;
    int resolution = QInputDialog::getInt(this, "Enter resolution",
                                          "Plot sampling resolution:", defaultResolution, 1, 21, 1, &ok);
    if (ok)
    {
        ui_->pushButton_start->setEnabled(false);
        computationIsRunning_ = true;
        computationManager_->startNewPlotComputation(selectedComputation_->input.computationNumber,
                                                     resolution, false);
    }
}



void MainWindow::on_pushButton_export_clicked()
{
    int computationNo = selectedComputation_->input.computationNumber;
    QString plotName = ui_->comboBox_fileExport->currentText();

    QString plotNameWithoutSpaces = plotName;
    plotNameWithoutSpaces.replace(QChar(' '), QChar('_'));
    QString suggestedFilename = QString("Computation_%1_%2.m").arg(computationNo).arg(plotNameWithoutSpaces);
    QDir directory(qApp->applicationDirPath());

    QString filename = QFileDialog::getSaveFileName(this, "Choose path and file name",
                                                    directory.absoluteFilePath(suggestedFilename),
                                                    "Matlab files (*.m)");
    // TODO: check Matlab file name conformity!

    if (filename.isEmpty())
        return;

    if (plotName == "Coefficients plot")
    {
        guiCommunicator_->exportIndexPlotToMatlabFile(computationNo, filename);
    }
    else
    {
        if (plotName == "Convergence logs")
        {
            guiCommunicator_->exportConvergenceLogsToMatlabFile(computationNo, filename);
        }
        else
        {
            if (plotName == "Solution plot")
            {
                guiCommunicator_->exportSolutionPlotToMatlabFile(computationNo, filename);
            }
            else
            {
                int optionalLogIndex = ui_->comboBox_fileExport->currentIndex() - 4;
                guiCommunicator_->exportOptionalLogToMatlabFile(computationNo, optionalLogIndex, filename);
            }
        }
    }
}



void MainWindow::on_pushButton_deleteComputationEntry_clicked()
{
    int computationNo = selectedComputation_->input.computationNumber;

    if (selectedComputation_->coeffState != CoeffComputationState::ERROR_PROBLEM_CREATION)
    {
        int answer = QMessageBox::question(this, "Delete computation entry?",
                                           QStringLiteral("Do you really want to delete the computation "
                                                          "entry for Computation %1?").arg(computationNo),
                                           QMessageBox::Yes, QMessageBox::No);
        if (answer == QMessageBox::No)
            return;
    }

    plotManager_->removeConvergencePlots(computationNo);
    plotManager_->removeSolutionPlot(computationNo);
    computationManager_->deleteComputation(computationNo);
    guiCommunicator_->deleteSolution(computationNo);
}



void MainWindow::on_checkBox_plotSolution_toggled(bool checked)
{
    ui_->widget_plotResolution->setEnabled(checked);
}



void MainWindow::on_actionClear_computation_history_triggered()
{
    int runningComputationNo = -1;
    QList<int> computationNoList = computationManager_->getComputationNumbers(runningComputationNo);

    if (computationNoList.isEmpty())
        return;

    ui_->tabWidget->setCurrentWidget(ui_->tabWidgetPage_ComputationHistory);

    QString question;

    if (computationIsRunning_)
        question = QStringLiteral("Do you really want to delete all computation "
                                  "entries (except for the running computation)?");
    else
        question = QStringLiteral("Do you really want to delete all computation "
                                  "entries?");

    int answer = QMessageBox::question(this, "Clear computation history?", question,
                                       QMessageBox::Yes, QMessageBox::No);
    if (answer == QMessageBox::No)
        return;

    for (int computationNo : computationNoList)
    {
        if (computationNo != runningComputationNo)
        {
            plotManager_->removeConvergencePlots(computationNo);
            plotManager_->removeSolutionPlot(computationNo);
            computationManager_->deleteComputation(computationNo);
            guiCommunicator_->deleteSolution(computationNo);
        }
    }
}



void MainWindow::on_actionDefining_custom_problems_triggered()
{
    if (!infoDefiningDialog_)
        infoDefiningDialog_ = new InfoDefiningDialog(this);

    infoDefiningDialog_->show();
}



void MainWindow::on_actionControls_triggered()
{
    if (!infoControlsDialog_)
        infoControlsDialog_ = new InfoControlsDialog(this);

    infoControlsDialog_->show();
}



void MainWindow::on_actionAbout_MSL_GUI_triggered()
{
    if (!infoAboutDialog_)
        infoAboutDialog_ = new InfoAboutDialog(this);

    infoAboutDialog_->show();
}
