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


#include "export_combobox.h"

#include <QStandardItemModel>

void ExportItemDelegate::paint(QPainter* painter, const QStyleOptionViewItem& option,
                               const QModelIndex& index) const
{
    QString type = index.data(Qt::AccessibleDescriptionRole).toString();
    if (type == QLatin1String("parent"))
    {
        QStyleOptionViewItem parentOption = option;
        parentOption.state |= QStyle::State_Enabled;
        QItemDelegate::paint(painter, parentOption, index);
    }
    else
    {
        if (type == QLatin1String("child"))
        {
            QStyleOptionViewItem childOption = option;
            int indent = option.fontMetrics.width(QString(4, QChar(' ')));
            childOption.rect.adjust(indent, 0, 0, 0);
            childOption.textElideMode = Qt::ElideNone;
            QItemDelegate::paint(painter, childOption, index);
        }
        else
        {
            QItemDelegate::paint(painter, option, index);
        }
    }
}




ExportComboBox::ExportComboBox(QWidget* parent)
    : QComboBox(parent)
{
    itemModel_ = qobject_cast<QStandardItemModel*>(model());

    ExportItemDelegate* itemDelegate = new ExportItemDelegate();
    setItemDelegate(itemDelegate);
//    addItems({"Coefficients plot", "Convergence logs", "Solution plot"});
}



void ExportComboBox::addParentItem(const QString& text)
{
    QStandardItem* item = new QStandardItem(text);
    item->setFlags(item->flags() & ~(Qt::ItemIsEnabled | Qt::ItemIsSelectable));
    item->setData("parent", Qt::AccessibleDescriptionRole);

    QFont font = item->font();
    font.setBold(true);
    item->setFont(font);

    itemModel_->appendRow(item);
}



void ExportComboBox::addChildItem(const QString& text)
{
    //QStandardItem* item = new QStandardItem(text + QString(4, QChar(' ')));
    QStandardItem* item = new QStandardItem(text);
    item->setData("child", Qt::AccessibleDescriptionRole);

    itemModel_->appendRow(item);
}



void ExportComboBox::setItems(CoeffComputationState::Enum coeffState, PlotComputationState::Enum plotState,
                              const QStringList& optionalLogs)
{
    clear();

    switch (coeffState)
    {
    case CoeffComputationState::RUNNING:
    case CoeffComputationState::ERROR_PROBLEM_CREATION:
    case CoeffComputationState::ERROR_NO_SOLUTION_CREATED:
        break;
    case CoeffComputationState::ERROR_COEFF_COMPUTATION:
    case CoeffComputationState::ABORTED:
    case CoeffComputationState::COMPLETE:
        addItems({"Coefficients plot", "Convergence logs", "Solution plot"});
        if (plotState != PlotComputationState::COMPLETE)
        {
            enableItem(2, false);
        }
        addOptionalLogs(optionalLogs);
        setCurrentIndex(0);
        break;
    }
}



void ExportComboBox::enableItem(int index, bool enable)
{
    QStandardItem* item = itemModel_->item(index);

    if (enable)
    {
       item->setFlags(item->flags() | Qt::ItemIsEnabled);
    }
    else
    {
       item->setFlags(item->flags() & ~Qt::ItemIsEnabled);
    }
}



void ExportComboBox::addOptionalLogs(const QStringList& optionalLogs)
{
    if (!optionalLogs.isEmpty())
    {
        addParentItem(QStringLiteral("Optional logs:"));
        for (const QString& log : optionalLogs)
            addChildItem(log);
    }
}
