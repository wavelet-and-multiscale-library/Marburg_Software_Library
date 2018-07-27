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


#ifndef COMPUTATION_TABLE_H
#define COMPUTATION_TABLE_H

#include <QTableWidget>
#include <QSignalMapper>

class QCheckBox;



class ComputationTable : public QTableWidget
{
    Q_OBJECT
public:
    ComputationTable(QWidget* parent = nullptr);

    static QColor getColorForComputation(int computationNo);

    void addEntry(int computationNo);
    void removeEntry(int computationNo);

    int getSelectedComputationNo() const;

public slots:
    void setStatusEntry(int computationNo, const QString& status);
    void setTimeEntry(int computationNo, const QString& time);
    void setPlotEntry(int computationNo, const QString& plotStatus);

    void checkAllCheckboxes();
    void uncheckAllCheckboxes();

signals:
    void selectedComputationChanged(int computationNo);
    void checkboxToogled(int computationNo);

private slots:
    void handleCurrentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

private:
    int getRowIndex(int computationNo) const;

    QList<int> computationNoList_;
    QList<QCheckBox*> checkBoxList_;
    QSignalMapper* checkboxMapper_;
};

#endif // COMPUTATION_TABLE_H
