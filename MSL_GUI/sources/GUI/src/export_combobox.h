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


#ifndef EXPORT_COMBOBOX_H
#define EXPORT_COMBOBOX_H


#include <QComboBox>
#include <QItemDelegate>
#include <QStandardItemModel>

#include "interfaces/computation_state_enums.h"


class ExportItemDelegate : public QItemDelegate
{
public:
    void paint(QPainter* painter, const QStyleOptionViewItem& option,
               const QModelIndex& index) const override;

};



class ExportComboBox : public QComboBox
{
public:
    ExportComboBox(QWidget* parent = nullptr);

    void addParentItem(const QString& text);
    void addChildItem(const QString& text);

    void setItems(CoeffComputationState::Enum coeffState, PlotComputationState::Enum plotState,
                  const QStringList& optionalLogs);

protected:
    void enableItem(int index, bool enable);
    void addOptionalLogs(const QStringList& optionalLogs);

    QStandardItemModel* itemModel_;
};

#endif // EXPORT_COMBOBOX_H
