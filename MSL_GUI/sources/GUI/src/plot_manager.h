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


#ifndef PLOT_MANAGER_H
#define PLOT_MANAGER_H

#include <QObject>

#include "convergence_chart.h"

namespace QtDataVisualization {
class Q3DSurface;
class QSurface3DSeries;
}

using std::map;
using std::vector;

class PlotManager : public QObject
{
    Q_OBJECT
public:
    explicit PlotManager(QObject* parent, QChartView* convergenceView, QChartView* convergenceTimeView,
                         QStackedWidget* solutionStackedWidget, QLabel* samplePointLabel);

    void addNewConvergencePlots(int computationNo, const QString& calloutText);
    void removeConvergencePlots(int computationNo);
    void removeSolutionPlot(int computationNo);

signals:

public slots:
    void addToCurrentConvPlot(double x, double y);
    void addToCurrentConvTimePlot(double x, double y);

    void toggleConvPlotsVisibility(int computationNo);

    void setup1DsolutionPlotPatch(int computationNo, const double* xValues, const double* yValues,
                                  int samplingPointCount, bool newPlot);

    void setup2DsolutionPlotPatch(int computationNo, const double* xValues, const double* yValues,
                                  const double* zValues, int sampleCountX, int sampleCountY, bool newPlot);

    void showSolutionPlot(int computationNo);

private:
    ConvergenceChart* convergenceChart_;
    ConvergenceChart* convergenceTimeChart_;

    QStackedWidget* solutionStackedWidget_;
    QLabel* samplePointLabel_;

    QChart* solution1Dchart_;
    // maps computation numbers to (as the case may be multi-patch and hence multi-series) 1D solution plots:
    map<int, vector<QLineSeries*> > solutionPlots1D_;
    int current1DsolutionPlotNumber_;

    QtDataVisualization::Q3DSurface* surfaceGraph_;
    // maps computation numbers to (as the case may be multi-patch and hence multi-series) 2D solution plots:
    map<int, vector< QtDataVisualization::QSurface3DSeries* > > solutionPlots2D_;
    int current2DsolutionPlotNumber_;
    QLinearGradient colorGradient_;

    map<int, int> samplePointCounts_;   // maps computation numbers to sample point counts

    void setupSurfaceGraph();
    void switchTo1DsolutionPlot(int computationNo);
    void switchTo2DsolutionPlot(int computationNo);
};

#endif // PLOT_MANAGER_H
