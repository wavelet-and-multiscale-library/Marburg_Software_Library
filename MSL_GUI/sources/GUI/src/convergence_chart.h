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


#ifndef CONVERGENCE_CHART_H
#define CONVERGENCE_CHART_H


#include <QtCharts>

using namespace QtCharts;



class ConvergenceChart : public QChart
{
    Q_OBJECT
public:
    ConvergenceChart(const QString& title, const QString& xAxisTitle, const QString& yAxisTitle,
                     qreal xMin, qreal xMax, qreal yMin, qreal yMax);

    void addNewSeries(int computationNo, const QString& calloutText, const QColor& color);
    void deleteSeries(int computationNo);
    void togglePlotVisibility(int computationNo);
    void adjustAxesRangesToVisibleSeries();

public slots:
    void addDataToCurrentSeries(double x, double y);

protected:
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event) override;

private:
    QLogValueAxis* axisX_;
    QLogValueAxis* axisY_;

    struct ConvergencePlot
    {
        QLineSeries* series;
        qreal xMin, xMax, yMin, yMax;
        //QString calloutText;
    };

    std::map<int, ConvergencePlot> convPlots_;    // maps computation numbers to convergence plots
    int currentComputationNo_;
    //const qreal xMinDefault_, xMaxDefault_, yMinDefault_, yMaxDefault_;


    class Callout : public QGraphicsItem
    {
    public:
        Callout(ConvergenceChart* parent, int pointsize);

        void setText(const QString &text);
        void setAnchor(QPointF point);
        void updateGeometry();

        QRectF boundingRect() const override;
        void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,QWidget *widget) override;

    protected:
        //void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event);
        void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
    private:
        QString m_text;
        QRectF m_textRect;
        QRectF m_rect;
        QPointF m_anchor;
        QFont m_font;
        ConvergenceChart* m_chart;
        int m_pointsize;
    };


    Callout* tooltip_;
    std::map<QObject*, QString> calloutTexts_;  // maps LineSeries to callout texts

private slots:
    void tooltip(QPointF point, bool show);
    void keepCallout();
    void removeCalloutFromScene(Callout* callout);
};


#endif // CONVERGENCE_CHART_H
