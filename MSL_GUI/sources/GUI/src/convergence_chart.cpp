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


#include "convergence_chart.h"

#include <cfloat>

ConvergenceChart::ConvergenceChart(const QString& title, const QString& xAxisTitle, const QString& yAxisTitle,
                                   qreal xMin, qreal xMax, qreal yMin, qreal yMax)
    : QChart(), tooltip_(nullptr)
{
    setTitle(title);

    axisX_ = new QLogValueAxis();
    axisX_->setTitleText(xAxisTitle);
    axisX_->setLabelFormat("%.0e");
    axisX_->setBase(10.0);
    axisX_->setMinorTickCount(-1);
    addAxis(axisX_, Qt::AlignBottom);

    axisY_ = new QLogValueAxis();
    axisY_->setTitleText(yAxisTitle);
    axisY_->setLabelFormat("%.0e");
    axisY_->setBase(10.0);
    axisY_->setMinorTickCount(-1);
    addAxis(axisY_, Qt::AlignLeft);

    axisX_->setRange(xMin, xMax);
    axisY_->setRange(yMin, yMax);

    setAcceptHoverEvents(true);
}



void ConvergenceChart::addNewSeries(int computationNo, const QString& calloutText, const QColor& color)
{
    currentComputationNo_ = computationNo;
    QLineSeries* newSeries = new QLineSeries();
    newSeries->setName(QString("Computation %1").arg(computationNo));
    newSeries->setColor(color);

    addSeries(newSeries);
    newSeries->attachAxis(axisX_);
    newSeries->attachAxis(axisY_);

    newSeries->setPointsVisible(true);
    QPen pen = newSeries->pen();
    pen.setWidth(2);
    newSeries->setPen(pen);

    ConvergencePlot newPlot;
    newPlot.series = newSeries;
    newPlot.xMin = DBL_MAX;
    newPlot.xMax = -DBL_MAX;
    newPlot.yMin = DBL_MAX;
    newPlot.yMax = -DBL_MAX;
    //newPlot.calloutText = calloutText;
    convPlots_[computationNo] = newPlot;

    calloutTexts_[newSeries] = calloutText;

    connect(newSeries, &QLineSeries::clicked, this, &ConvergenceChart::keepCallout);
    connect(newSeries, &QLineSeries::hovered, this, &ConvergenceChart::tooltip);
}



void ConvergenceChart::deleteSeries(int computationNo)
{
    removeSeries(convPlots_.at(computationNo).series);
    delete convPlots_.at(computationNo).series;
    convPlots_.erase(computationNo);
}



void ConvergenceChart::togglePlotVisibility(int computationNo)
{
    QLineSeries* series = convPlots_.at(computationNo).series;
    series->setVisible(!series->isVisible());
}



void ConvergenceChart::adjustAxesRangesToVisibleSeries()
{
    qreal xMin = DBL_MAX;
    qreal xMax = -DBL_MAX;
    qreal yMin = DBL_MAX;
    qreal yMax = -DBL_MAX;

    for(const auto& plotPair : convPlots_)
    {
        const ConvergencePlot& plot = plotPair.second;
        if (plot.series->isVisible())
        {
            if (plot.xMin < xMin)
                xMin = plot.xMin;
            if (plot.xMax > xMax)
                xMax = plot.xMax;
            if (plot.yMin < yMin)
                yMin = plot.yMin;
            if (plot.yMax > yMax)
                yMax = plot.yMax;
        }
    }

    if (xMax > -DBL_MAX)
    {
        axisX_->setRange(xMin, xMax);
        axisY_->setRange(yMin, yMax);
    }
}



void ConvergenceChart::addDataToCurrentSeries(double x, double y)
{
    if (x <= 0.0 || y <= 0.0)
        return;

    ConvergencePlot& currentPlot = convPlots_.at(currentComputationNo_);
    currentPlot.series->append(x, y);

    if (x < currentPlot.xMin)
        currentPlot.xMin = x;
    if (x > currentPlot.xMax)
        currentPlot.xMax = x;
    if (y < currentPlot.yMin)
        currentPlot.yMin = y;
    if (y > currentPlot.yMax)
        currentPlot.yMax = y;

    if (currentPlot.series->isVisible())
    {
        if (x < axisX_->min())
            axisX_->setMin(x);
        if (x > axisX_->max())
            axisX_->setMax(x);
        if (y < axisY_->min())
            axisY_->setMin(y);
        if (y > axisY_->max())
            axisY_->setMax(y);
    }

}




void ConvergenceChart::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event)
{
    adjustAxesRangesToVisibleSeries();
    event->setAccepted(true);
}




void ConvergenceChart::tooltip(QPointF point, bool show)
{
    if (tooltip_ == nullptr)
        tooltip_ = new Callout(this, 8);

    if (show)
    {
        tooltip_->setText(calloutTexts_.at(sender()));
        tooltip_->setAnchor(point);
        tooltip_->setZValue(11);
        tooltip_->updateGeometry();
        tooltip_->show();
    }
    else
    {
        tooltip_->hide();
    }
}



void ConvergenceChart::keepCallout()
{
    tooltip_ = new Callout(this, 8);
}



void ConvergenceChart::removeCalloutFromScene(Callout* callout)
{
    scene()->removeItem(callout);
}





// Implementation for class ConvergenceChart::Callout


ConvergenceChart::Callout::Callout(ConvergenceChart* parent, int pointsize)
    : QGraphicsItem(parent),
      m_anchor(1.0, 1.0),
      m_chart(parent),
      m_pointsize(pointsize)
{
    m_font.setPointSize(pointsize);
    setFlag(QGraphicsItem::ItemIsMovable);
}

QRectF ConvergenceChart::Callout::boundingRect() const
{
    QPointF anchor = mapFromParent(m_chart->mapToPosition(m_anchor));
    QRectF rect;
    rect.setLeft(qMin(m_rect.left(), anchor.x()));
    rect.setRight(qMax(m_rect.right(), anchor.x()));
    rect.setTop(qMin(m_rect.top(), anchor.y()));
    rect.setBottom(qMax(m_rect.bottom(), anchor.y()));
    return rect;
}

void ConvergenceChart::Callout::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    Q_UNUSED(option)
    Q_UNUSED(widget)
    QPainterPath path;
    path.addRoundedRect(m_rect, 5, 5);

    QPointF anchor = mapFromParent(m_chart->mapToPosition(m_anchor));
    if (!m_rect.contains(anchor)) {
        QPointF point1, point2;

        // establish the position of the anchor point in relation to m_rect
        bool above = anchor.y() <= m_rect.top();
        bool aboveCenter = anchor.y() > m_rect.top() && anchor.y() <= m_rect.center().y();
        bool belowCenter = anchor.y() > m_rect.center().y() && anchor.y() <= m_rect.bottom();
        bool below = anchor.y() > m_rect.bottom();

        bool onLeft = anchor.x() <= m_rect.left();
        bool leftOfCenter = anchor.x() > m_rect.left() && anchor.x() <= m_rect.center().x();
        bool rightOfCenter = anchor.x() > m_rect.center().x() && anchor.x() <= m_rect.right();
        bool onRight = anchor.x() > m_rect.right();

        // get the nearest m_rect corner.
        qreal x = (onRight + rightOfCenter) * m_rect.width();
        qreal y = (below + belowCenter) * m_rect.height();
        bool cornerCase = (above && onLeft) || (above && onRight) || (below && onLeft) || (below && onRight);
        bool vertical = qAbs(anchor.x() - x) > qAbs(anchor.y() - y);

        qreal x1 = x + leftOfCenter * 10 - rightOfCenter * 20 + cornerCase * !vertical * (onLeft * 10 - onRight * 20);
        qreal y1 = y + aboveCenter * 10 - belowCenter * 20 + cornerCase * vertical * (above * 10 - below * 20);;
        point1.setX(x1);
        point1.setY(y1);

        qreal x2 = x + leftOfCenter * 20 - rightOfCenter * 10 + cornerCase * !vertical * (onLeft * 20 - onRight * 10);;
        qreal y2 = y + aboveCenter * 20 - belowCenter * 10 + cornerCase * vertical * (above * 20 - below * 10);;
        point2.setX(x2);
        point2.setY(y2);

        path.moveTo(point1);
        path.lineTo(anchor);
        path.lineTo(point2);
        path = path.simplified();
    }
    painter->setBrush(QColor(255, 255, 255));
    painter->drawPath(path);

    QFont font = painter->font();
    font.setPointSize(m_pointsize);
    painter->setFont(font);

    painter->drawText(m_textRect, Qt::TextShowMnemonic, m_text);
}

//void ConvergenceChart::Callout::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event)
//{
//    m_chart->removeCalloutFromScene(this);
//    event->setAccepted(true);
//}

void ConvergenceChart::Callout::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    if ((event->buttons() & Qt::RightButton)
         && event->modifiers().testFlag(Qt::ControlModifier))
    {
        m_chart->removeCalloutFromScene(this);
    }
    event->setAccepted(true);
}


void ConvergenceChart::Callout::setText(const QString &text)
{
    QFontMetrics metrics(m_font);
    QStringList textLines = text.split('\n');
    for (int i = 0; i < textLines.length(); ++i)
        textLines[i] = metrics.elidedText(textLines.at(i), Qt::ElideRight, 230, Qt::TextShowMnemonic);
    m_text = textLines.join('\n');

    m_textRect = metrics.boundingRect(QRect(0, 0, 150, 150), Qt::AlignLeft | Qt::TextShowMnemonic, m_text);
    m_textRect.translate(5, 5);
    prepareGeometryChange();
    m_rect = m_textRect.adjusted(-5, -5, 5, 5);
}

void ConvergenceChart::Callout::setAnchor(QPointF point)
{
    m_anchor = point;
}

void ConvergenceChart::Callout::updateGeometry()
{
    prepareGeometryChange();
    setPos(m_chart->mapToPosition(m_anchor) + QPoint(10, -50));
}
