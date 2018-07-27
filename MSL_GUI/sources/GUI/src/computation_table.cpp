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


#include "computation_table.h"

#include <QHeaderView>
#include <QHBoxLayout>
#include <QCheckBox>


ComputationTable::ComputationTable(QWidget* parent) :
    QTableWidget(0, 5, parent), checkboxMapper_(new QSignalMapper(this))
{
    setSelectionMode(QAbstractItemView::SingleSelection);
    setSelectionBehavior(QAbstractItemView::SelectRows);

    setStyleSheet("QTableWidget::item{ selection-color: black; selection-background-color: rgb(230,240,255)}");

    setHorizontalHeaderLabels({"Computation", "Show\nconvergence", "Status", "Time", "Plot\n(resolution)"});
    horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
    horizontalHeader()->setSectionResizeMode(3, QHeaderView::ResizeToContents);
    horizontalHeader()->setSectionResizeMode(4, QHeaderView::Stretch);
    horizontalHeader()->setHighlightSections(false);

    connect(checkboxMapper_, SIGNAL(mapped(int)), this, SIGNAL(checkboxToogled(int)));
    connect(this, &QTableWidget::currentCellChanged, this, &ComputationTable::handleCurrentCellChanged);
}



QColor ComputationTable::getColorForComputation(int computationNo)
{
    switch (computationNo%10)
    {
    case 1:
        return QColor(Qt::red);
    case 2:
        return QColor(Qt::blue);
    case 3:
        return QColor(255, 175, 0);
    case 4:
        return QColor(Qt::green);
    case 5:
        return QColor(Qt::magenta);
    case 6:
        return QColor(Qt::cyan);
    case 7:
        return QColor(Qt::darkGray);
    case 8:
        return QColor(Qt::darkRed);
    case 9:
        return QColor(Qt::darkBlue);
    case 0:
        return QColor(Qt::darkGreen);
    default:
        return QColor(Qt::black);
    }
}



void ComputationTable::addEntry(int computationNo)
{
    setRowCount(rowCount()+1);
    int rowIndex = rowCount() - 1;

    computationNoList_.append(computationNo);

    // First column:
    QPixmap pixmap(32,32);
    QColor iconColor = getColorForComputation(computationNo);
    pixmap.fill(iconColor);
    QIcon icon(pixmap);
    QTableWidgetItem* compItem = new QTableWidgetItem(icon, QString("Computation %1").arg(computationNo));
    compItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    setItem(rowIndex, 0, compItem);

    // Second column:
    QCheckBox* checkbox = new QCheckBox();
    checkbox->setChecked(true);
    checkboxMapper_->setMapping(checkbox, computationNo);
    connect(checkbox, SIGNAL(stateChanged(int)), checkboxMapper_, SLOT(map()));
    checkBoxList_.append(checkbox);

    QWidget* widget = new QWidget();
    QHBoxLayout* layout = new QHBoxLayout();
    layout->addWidget(checkbox);
    layout->setAlignment(Qt::AlignCenter);
    layout->setContentsMargins(0, 0, 0, 0);
    widget->setLayout(layout);
    setCellWidget(rowIndex, 1, widget);

    // Third column:
    QTableWidgetItem* statusItem = new QTableWidgetItem("Computing...");
    statusItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    statusItem->setTextAlignment(Qt::AlignCenter);
    setItem(rowIndex, 2, statusItem);

    // Fourth column:
    QTableWidgetItem* timeItem = new QTableWidgetItem("00:00");
    timeItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    timeItem->setTextAlignment(Qt::AlignCenter);
    setItem(rowIndex, 3, timeItem);

    // Fifth column:
    QTableWidgetItem* resolutionItem = new QTableWidgetItem("n/a");
    resolutionItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    resolutionItem->setTextAlignment(Qt::AlignCenter);
    setItem(rowIndex, 4, resolutionItem);

    setCurrentCell(rowIndex, 0);
}



void ComputationTable::removeEntry(int computationNo)
{
    int rowIndex = getRowIndex(computationNo);

    if (rowIndex >= 0 && rowIndex < rowCount())
    {
        removeRow(rowIndex);
        computationNoList_.removeAt(rowIndex);
        checkBoxList_.removeAt(rowIndex);
    }
}



int ComputationTable::getSelectedComputationNo() const
{
    int selectedRow = currentRow();

    if (selectedRow >= 0 && selectedRow < rowCount())
    {
        return computationNoList_.at(selectedRow);
    }
    else
    {
        return -1;
    }
}



void ComputationTable::setStatusEntry(int computationNo, const QString& status)
{
    item(getRowIndex(computationNo), 2)->setText(status);
}



void ComputationTable::setTimeEntry(int computationNo, const QString& time)
{
    item(getRowIndex(computationNo), 3)->setText(time);
}



void ComputationTable::setPlotEntry(int computationNo, const QString& plotStatus)
{
    item(getRowIndex(computationNo), 4)->setText(plotStatus);
}



void ComputationTable::checkAllCheckboxes()
{
    for (auto checkbox : checkBoxList_)
        checkbox->setChecked(true);
}



void ComputationTable::uncheckAllCheckboxes()
{
    for (auto checkbox : checkBoxList_)
        checkbox->setChecked(false);
}



void ComputationTable::handleCurrentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn)
{
    (void) currentColumn;
    (void) previousColumn;

    if (currentRow != previousRow)
    {
        emit selectedComputationChanged(getSelectedComputationNo());
    }
}



int ComputationTable::getRowIndex(int computationNo) const
{
    return computationNoList_.indexOf(computationNo);
}
