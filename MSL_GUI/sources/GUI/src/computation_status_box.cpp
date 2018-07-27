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


#include <QLabel>
#include <QPushButton>
#include <QMessageBox>
#include <QTimer>

#include "computation_status_box.h"



ComputationStatusBox::ComputationStatusBox(QWidget* parent) :
    QGroupBox(parent),
    timer_(new QTimer(this)),
    zeroTime_(0, 0, 0),
    seconds_(0),
    coeffComputationIsRunning_(false),
    currentComputationNo_(-1)
{
    connect(timer_, &QTimer::timeout, this, &ComputationStatusBox::updateElapsedTime);
}



void ComputationStatusBox::setWidgets(QWidget* coeffStatusWidget, QWidget* plotStatusWidget)
{
    coeffStatusWidget_ = coeffStatusWidget;
    plotStatusWidget_ = plotStatusWidget;
}



void ComputationStatusBox::setAbortButton(QPushButton* abortButton)
{
    abortButton_ = abortButton;
    connect(abortButton_, &QPushButton::clicked, this, &ComputationStatusBox::handleAbortClicked);
    abortButton_->setVisible(false);
}



void ComputationStatusBox::setLabels(QLabel* coeffStatusLabel, QLabel* coeffTimeLabel,
                                     QLabel* plotStatusLabel, QLabel* plotTimeLabel)
{
    coeffStatusLabel_ = coeffStatusLabel;
    coeffTimeLabel_ = coeffTimeLabel;
    plotStatusLabel_ = plotStatusLabel;
    plotTimeLabel_ = plotTimeLabel;

    coeffStatusLabel_->setText("No computation started.");
    coeffTimeLabel_->setVisible(false);

    plotStatusLabel_->setVisible(false);
    plotTimeLabel_->setVisible(false);
}



void ComputationStatusBox::showCoefficientsComputing(int computationNo)
{
    setTitle(QString("Status (Computation %1)").arg(computationNo));
    plotStatusLabel_->setVisible(false);
    plotTimeLabel_->setVisible(false);

    coeffStatusLabel_->setText("▶  Computing solution coefficients...");
    abortButton_->setEnabled(true);
    abortButton_->setVisible(true);
    coeffTimeLabel_->setText("00:00");
    coeffTimeLabel_->setVisible(true);

    coeffComputationIsRunning_ = true;
    currentComputationNo_ = computationNo;
    restartTimer();
}



void ComputationStatusBox::showCoefficientsComputingEnded(const QString& endStatus)
{
    stopTimer();
    coeffStatusLabel_->setText(endStatus);
    abortButton_->setVisible(false); 
    coeffComputationIsRunning_ = false;
}



void ComputationStatusBox::showPlotComputing(int computationNo, int resolution,
                                             bool subsequent, const QString& coeffEndState)
{
    if (!subsequent)
    {
        setTitle(QString("Status (Computation %1)").arg(computationNo));
        coeffStatusLabel_->setText(coeffEndState);
        coeffTimeLabel_->setVisible(false);
    }

    coeffStatusWidget_->setEnabled(false);

    plotStatusLabel_->setText(QString("▶  Computing plot samples... (resolution: %1)").arg(resolution));
    plotTimeLabel_->setText("00:00");
    plotStatusLabel_->setVisible(true);
    plotTimeLabel_->setVisible(true);

    restartTimer();
}



void ComputationStatusBox::showPlotComputingEnded(const QString& endStatus)
{
    stopTimer();
    plotStatusLabel_->setText(endStatus);
    coeffStatusWidget_->setEnabled(true);
}



QString ComputationStatusBox::getCoeffTimeString() const
{
    return coeffTimeLabel_->text();
}



void ComputationStatusBox::updateElapsedTime()
{
    seconds_++;
    QTime currentTime(zeroTime_.addSecs(seconds_));

    QString time = (seconds_ < 3600) ? currentTime.toString("mm:ss") : currentTime.toString("H:mm:ss");

    if (coeffComputationIsRunning_)
    {
        coeffTimeLabel_->setText(time);
        emit coeffComputationTimeChanged(currentComputationNo_, time);
    }
    else
    {
        plotTimeLabel_->setText(time);
    }
}



void ComputationStatusBox::handleAbortClicked()
{
    QMessageBox::StandardButton button = QMessageBox::question(this, QStringLiteral("Abort computation?"),
                                                               QStringLiteral("Do you really want to abort the current computation?"));
    if (button == QMessageBox::Yes)
    {
        abortButton_->setEnabled(false);
        coeffStatusLabel_->setText("▶  Abort requested. Waiting for computation to stop...");
        emit abortRequested();
    }
}



void ComputationStatusBox::stopTimer()
{
    timer_->stop();
}



void ComputationStatusBox::restartTimer()
{
    seconds_ = 0;
    timer_->start(1000);
}
