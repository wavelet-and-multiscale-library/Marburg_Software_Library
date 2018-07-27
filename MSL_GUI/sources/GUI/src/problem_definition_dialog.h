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


#ifndef PROBLEM_DEFINITION_DIALOG_H
#define PROBLEM_DEFINITION_DIALOG_H

#include <QDialog>

class QLabel;
class QLineEdit;
class InfoDefiningDialog;

namespace Ui {
class ProblemDefinitionDialog;
}

class ProblemDefinitionDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ProblemDefinitionDialog(QWidget* parent, bool saveToFileAllowed);
    ~ProblemDefinitionDialog();

    void adaptTo(const QStringList& functionNames);
    void prepareEditing(const QString& currentProblemName,
                        const QStringList& currentProblemFuncDefs,
                        bool revertEnabled);
    void prepareNewDefinition();

signals:
    void saveDefinitionRequested(bool overwrite, bool saveToFile, const QString& oldName,
                                 const QString& newName, const QStringList& newFuncDefs) const;
    void revertToFileStateRequested() const;

private slots:
    void enableSaveButtons();
    void handleSaveButtonClicked() const;

    void on_pushButton_cancel_clicked();

    void on_pushButton_help_clicked();

private:
    Ui::ProblemDefinitionDialog* ui_;
    InfoDefiningDialog* infoDefiningDialog_;

    std::array<QLabel*, 4> functionLabels_;
    std::array<QLineEdit*, 4> functionLineEdits_;

    bool lastOpenedForEditing_;
    QString oldName_;
    const bool saveToFileAllowed_;
};

#endif // PROBLEM_DEFINITION_DIALOG_H
