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


#ifndef COMPUTATION_STATUS_BOX_H
#define COMPUTATION_STATUS_BOX_H

#include <QGroupBox>
#include <QTime>

class QLabel;
class QPushButton;


class ComputationStatusBox : public QGroupBox
{
    Q_OBJECT
public:
    explicit ComputationStatusBox(QWidget* parent = nullptr);

    void setWidgets(QWidget* coeffStatusWidget, QWidget* plotStatusWidget);
    void setAbortButton(QPushButton* abortButton);
    void setLabels(QLabel* coeffStatusLabel, QLabel* coeffTimeLabel,
                   QLabel* plotStatusLabel, QLabel* plotTimeLabel);

    void showCoefficientsComputing(int computationNo);
    void showCoefficientsComputingEnded(const QString& endStatus);
    void showPlotComputing(int computationNo, int resolution, bool subsequent, const QString& coeffEndState);
    void showPlotComputingEnded(const QString& endStatus);

    QString getCoeffTimeString() const;

signals:
    void abortRequested();
    void coeffComputationTimeChanged(int computationNo, const QString& coeffTime);

private slots:
    void updateElapsedTime();
    void handleAbortClicked();


private:
    void stopTimer();
    void restartTimer();

    QTimer* timer_;
    const QTime zeroTime_;
    int seconds_;

    bool coeffComputationIsRunning_;
    int currentComputationNo_;

    QWidget* coeffStatusWidget_;
    QWidget* plotStatusWidget_;

    QLabel* coeffStatusLabel_;
    QLabel* coeffTimeLabel_;
    QLabel* plotStatusLabel_;
    QLabel* plotTimeLabel_;

    QPushButton* abortButton_;
};

#endif // COMPUTATION_STATUS_BOX_H
