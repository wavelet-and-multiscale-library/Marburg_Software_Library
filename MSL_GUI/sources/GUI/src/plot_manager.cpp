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


#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurface3DSeries>
using namespace QtDataVisualization;

#include "computation_table.h"
#include "plot_manager.h"


PlotManager::PlotManager(QObject* parent, QChartView* convergenceView, QChartView* convergenceTimeView,
                         QStackedWidget* solutionStackedWidget, QLabel* samplePointLabel)
    : QObject(parent), solutionStackedWidget_(solutionStackedWidget), samplePointLabel_(samplePointLabel),
      current1DsolutionPlotNumber_(-1), current2DsolutionPlotNumber_(-1)
{
    convergenceChart_ = new ConvergenceChart("Convergence", "# degrees of freedom",
                                             "(residual) error", 1.0, 1000.0, 1e-2, 1.0);
    convergenceView->setChart(convergenceChart_);

    convergenceTimeChart_ = new ConvergenceChart("Temporal convergence", "CPU time (in seconds)",
                                                 "(residual) error", 1.0, 100.0, 1e-2, 1.0);
    convergenceTimeView->setChart(convergenceTimeChart_);


    solution1Dchart_ = new QChart();
    solution1Dchart_->setTitle("Solution");
    solution1Dchart_->legend()->hide();
    QChartView* solution1DchartView = new QChartView(solution1Dchart_);
    solutionStackedWidget_->insertWidget(1, solution1DchartView);

    setupSurfaceGraph();
}



void PlotManager::addNewConvergencePlots(int computationNo, const QString& calloutText)
{
    QString calloutHeader;
    if (computationNo < 10)
    {
        calloutHeader = QString("&C&o&m&p&u&t&a&t&i&o&n& &%1\n\n").arg(computationNo);
    }
    else
    {
        calloutHeader = QString("&C&o&m&p&u&t&a&t&i&o&n& &%1&%2\n\n").arg(computationNo/10).arg(computationNo%10);
    }

    QString callout = calloutHeader + calloutText;

    QColor plotColor = ComputationTable::getColorForComputation(computationNo);
    convergenceChart_->addNewSeries(computationNo, callout, plotColor);
    convergenceTimeChart_->addNewSeries(computationNo, callout, plotColor);
}



void PlotManager::removeConvergencePlots(int computationNo)
{
    convergenceChart_->deleteSeries(computationNo);
    convergenceTimeChart_->deleteSeries(computationNo);
}



void PlotManager::removeSolutionPlot(int computationNo)
{
    if (computationNo == current1DsolutionPlotNumber_)
    {
        for (auto patch : solutionPlots1D_.at(computationNo))
            solution1Dchart_->removeSeries(patch);
        current1DsolutionPlotNumber_ = -1;
    }
    if (computationNo == current2DsolutionPlotNumber_)
    {
        for (auto patch : solutionPlots2D_.at(computationNo))
            surfaceGraph_->removeSeries(patch);
        current2DsolutionPlotNumber_ = -1;
    }

    if (solutionPlots1D_.count(computationNo))
    {
        for (auto patch : solutionPlots1D_.at(computationNo))
            delete patch;
        solutionPlots1D_.erase(computationNo);
    }
    if (solutionPlots2D_.count(computationNo))
    {
        for (auto patch : solutionPlots2D_.at(computationNo))
            delete patch;
        solutionPlots2D_.erase(computationNo);
    }
}



void PlotManager::addToCurrentConvPlot(double x, double y)
{
    convergenceChart_->addDataToCurrentSeries(x, y);
}



void PlotManager::addToCurrentConvTimePlot(double x, double y)
{
    convergenceTimeChart_->addDataToCurrentSeries(x, y);
}



void PlotManager::toggleConvPlotsVisibility(int computationNo)
{
    convergenceChart_->togglePlotVisibility(computationNo);
    convergenceTimeChart_->togglePlotVisibility(computationNo);
}



void PlotManager::setup1DsolutionPlotPatch(int computationNo, const double* xValues, const double* yValues,
                                           int samplingPointCount, bool newPlot)
{
    QLineSeries* solPlotPatch = new QLineSeries();
    solPlotPatch->setColor(ComputationTable::getColorForComputation(computationNo));
    QPen pen = solPlotPatch->pen();
    pen.setWidth(3);
    solPlotPatch->setPen(pen);

    //  Limit the number of sampling points:
    //    int d;
    //    if (samplingPointCount <= 2049)
    //    {
    //        d = 1;
    //    }
    //    else
    //    {
    //        d = (samplingPointCount - 1)/2048;
    //    }
    //    int i;
    //    for (i = 0; i < samplingPointCount; i += d)
    //    {
    //        solPlotPatch->append(xValues[i], yValues[i]);
    //    }
    //    solPlotPatch->setName(QString("Computation %1").arg(computationNo);


    for (int i = 0; i < samplingPointCount; i++)
    {
        solPlotPatch->append(xValues[i], yValues[i]);
    }

    solPlotPatch->setName(QString("Computation %1").arg(computationNo));

    if (solutionPlots1D_.count(computationNo) > 0)
    {
        if (newPlot)
        {
            if (computationNo == current1DsolutionPlotNumber_)
            {
                for (auto patch : solutionPlots1D_.at(computationNo))
                {
                    solution1Dchart_->removeSeries(patch);
                }

                current1DsolutionPlotNumber_ = -1;
            }

            for (auto patch : solutionPlots1D_.at(computationNo))
            {
                delete patch;
            }
            solutionPlots1D_.erase(computationNo);
            solutionPlots1D_[computationNo] = vector<QLineSeries*>(1, solPlotPatch);

            samplePointCounts_[computationNo] = 0;
        }
        else
        {
            solutionPlots1D_.at(computationNo).push_back(solPlotPatch);
        }
    }
    else
    {
        solutionPlots1D_[computationNo] = vector<QLineSeries*>(1, solPlotPatch);
        samplePointCounts_[computationNo] = 0;
    }

    samplePointCounts_[computationNo] += samplingPointCount;
}



void PlotManager::setup2DsolutionPlotPatch(int computationNo, const double* xValues, const double* yValues,
                                           const double* zValues, int sampleCountX, int sampleCountY, bool newPlot)
{
    QSurface3DSeries* solPlotPatch = new QSurface3DSeries();
    solPlotPatch->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
    solPlotPatch->setFlatShadingEnabled(true);

    solPlotPatch->setBaseGradient(colorGradient_);
    solPlotPatch->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

    QSurfaceDataArray* dataArray = new QSurfaceDataArray;
    dataArray->reserve(sampleCountY);

    for (int row = 0; row < sampleCountY; row++)
    {
        QSurfaceDataRow* newRow = new QSurfaceDataRow(sampleCountX);
        float y = float (yValues[row]);
        for (int col = 0; col < sampleCountX; col++)
        {
            float x = float (xValues[col*sampleCountY]);
            float z = float (zValues[col*sampleCountY + row]);
            (*newRow)[col].setPosition(QVector3D(x, z, y));
        }
        *dataArray << newRow;
    }

    solPlotPatch->dataProxy()->resetArray(dataArray);

    if (solutionPlots2D_.count(computationNo) > 0)
    {
        if (newPlot)
        {
            if (computationNo == current2DsolutionPlotNumber_)
            {
                for (auto patch : solutionPlots2D_.at(computationNo))
                {
                    surfaceGraph_->removeSeries(patch);
                }

                current2DsolutionPlotNumber_ = -1;
            }

            for (auto patch : solutionPlots2D_.at(computationNo))
            {
                delete patch;
            }
            solutionPlots2D_.erase(computationNo);
            solutionPlots2D_[computationNo] = vector<QSurface3DSeries*>(1, solPlotPatch);

            samplePointCounts_[computationNo] = 0;
        }
        else
        {
            solutionPlots2D_.at(computationNo).push_back(solPlotPatch);
        }
    }
    else
    {
        solutionPlots2D_[computationNo] = vector<QSurface3DSeries*>(1, solPlotPatch);
        samplePointCounts_[computationNo] = 0;
    }

    samplePointCounts_[computationNo] += sampleCountX*sampleCountY;
}




void PlotManager::showSolutionPlot(int computationNo)
{
    if (solutionPlots1D_.count(computationNo) > 0)
    {
        switchTo1DsolutionPlot(computationNo);
        samplePointLabel_->setText(QStringLiteral("Computed sample points: %1").arg(samplePointCounts_[computationNo]));
    }
    else
    {
        if (solutionPlots2D_.count(computationNo) > 0)
        {
            switchTo2DsolutionPlot(computationNo);
            samplePointLabel_->setText(QStringLiteral("Computed sample points: %1").arg(samplePointCounts_[computationNo]));
        }
        else
        {
            solutionStackedWidget_->setCurrentIndex(0);
            samplePointLabel_->setText("");
        }
    }
}



void PlotManager::setupSurfaceGraph()
{
    surfaceGraph_ = new Q3DSurface();
//    surfaceGraph_->activeTheme()->setType(Q3DTheme::ThemeArmyBlue);
    surfaceGraph_->activeTheme()->setType(Q3DTheme::ThemeDigia);

    surfaceGraph_->setAxisX(new QValue3DAxis);
    surfaceGraph_->setAxisY(new QValue3DAxis);
    surfaceGraph_->setAxisZ(new QValue3DAxis);

    surfaceGraph_->axisX()->setLabelFormat("%.2f");
    surfaceGraph_->axisZ()->setLabelFormat("%.2f");
    surfaceGraph_->axisX()->setLabelAutoRotation(30);
    surfaceGraph_->axisY()->setLabelAutoRotation(90);
    surfaceGraph_->axisZ()->setLabelAutoRotation(30);

    colorGradient_.setColorAt(0.0, Qt::black);
    colorGradient_.setColorAt(0.33, Qt::blue);
    colorGradient_.setColorAt(0.67, Qt::red);
    colorGradient_.setColorAt(1.0, Qt::yellow);

    QWidget* container = QWidget::createWindowContainer(surfaceGraph_);
    solutionStackedWidget_->insertWidget(2, container);
}



void PlotManager::switchTo1DsolutionPlot(int computationNo)
{
    if (computationNo != current1DsolutionPlotNumber_)
    {
        if (current1DsolutionPlotNumber_ > 0)
        {    
            for (auto patch : solutionPlots1D_.at(current1DsolutionPlotNumber_))
                solution1Dchart_->removeSeries(patch);
        }

        for (auto patch : solutionPlots1D_.at(computationNo))
            solution1Dchart_->addSeries(patch);
        current1DsolutionPlotNumber_ = computationNo;

        solution1Dchart_->createDefaultAxes();
    }
    solutionStackedWidget_->setCurrentIndex(1);
}



void PlotManager::switchTo2DsolutionPlot(int computationNo)
{
    if (computationNo != current2DsolutionPlotNumber_)
    {
        if (current2DsolutionPlotNumber_ > 0)
        {
            for (auto patch : solutionPlots2D_.at(current2DsolutionPlotNumber_))
                surfaceGraph_->removeSeries(patch);
        }

        for (auto patch : solutionPlots2D_.at(computationNo))
            surfaceGraph_->addSeries(patch);
        current2DsolutionPlotNumber_ = computationNo;
    }
    solutionStackedWidget_->setCurrentIndex(2);
}
