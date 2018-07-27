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


#include "info_defining_dialog.h"
#include "problem_definition_dialog.h"
#include "ui_problem_definition_dialog.h"

ProblemDefinitionDialog::ProblemDefinitionDialog(QWidget* parent, bool saveToFileAllowed) :
    QDialog(parent),
    ui_(new Ui::ProblemDefinitionDialog),
    infoDefiningDialog_(nullptr),
    saveToFileAllowed_(saveToFileAllowed)
{
    ui_->setupUi(this);
    ui_->pushButton_saveToFile->setEnabled(saveToFileAllowed);

    functionLabels_.at(0) = ui_->label_f0;
    functionLabels_.at(1) = ui_->label_f1;
    functionLabels_.at(2) = ui_->label_f2;
    functionLabels_.at(3) = ui_->label_f3;

    functionLineEdits_.at(0) = ui_->lineEdit_f0;
    functionLineEdits_.at(1) = ui_->lineEdit_f1;
    functionLineEdits_.at(2) = ui_->lineEdit_f2;
    functionLineEdits_.at(3) = ui_->lineEdit_f3;

    connect(ui_->pushButton_revert, &QPushButton::clicked,
            this, &ProblemDefinitionDialog::revertToFileStateRequested);
    connect(ui_->pushButton_saveTemporarily, &QPushButton::clicked,
            this, &ProblemDefinitionDialog::handleSaveButtonClicked);
    connect(ui_->pushButton_saveToFile, &QPushButton::clicked,
            this, &ProblemDefinitionDialog::handleSaveButtonClicked);

    connect(ui_->lineEdit_problemName, &QLineEdit::textEdited,
            this, &ProblemDefinitionDialog::enableSaveButtons);
    connect(ui_->lineEdit_f0, &QLineEdit::textEdited,
            this, &ProblemDefinitionDialog::enableSaveButtons);
    connect(ui_->lineEdit_f1, &QLineEdit::textEdited,
            this, &ProblemDefinitionDialog::enableSaveButtons);
    connect(ui_->lineEdit_f2, &QLineEdit::textEdited,
            this, &ProblemDefinitionDialog::enableSaveButtons);
    connect(ui_->lineEdit_f3, &QLineEdit::textEdited,
            this, &ProblemDefinitionDialog::enableSaveButtons);

    setWindowModality(Qt::ApplicationModal);
}

ProblemDefinitionDialog::~ProblemDefinitionDialog()
{
    delete ui_;
}



void ProblemDefinitionDialog::adaptTo(const QStringList& functionNames)
{
    unsigned int i = 0;
    for (const QString& funcName : functionNames)
    {
        functionLabels_.at(i)->setVisible(true);
        functionLineEdits_.at(i)->setVisible(true);
        functionLabels_.at(i)->setText(funcName + " = ");
        functionLineEdits_.at(i)->clear();
        i++;
    }

    for (unsigned int j = i; j < functionLabels_.size(); j++)
    {
        functionLabels_.at(j)->setVisible(false);
        functionLineEdits_.at(j)->setVisible(false);
    }
}



void ProblemDefinitionDialog::prepareEditing(const QString& currentProblemName,
                                             const QStringList& currentProblemFuncDefs,
                                             bool revertEnabled)
{
    setWindowTitle("Edit problem definition");
    oldName_ = currentProblemName;

    for (int i = 0; i < currentProblemFuncDefs.size(); ++i)
        functionLineEdits_.at(i)->setText(currentProblemFuncDefs.at(i));

    ui_->lineEdit_problemName->setText(currentProblemName);

    ui_->pushButton_saveTemporarily->setEnabled(false);
    ui_->pushButton_saveToFile->setEnabled(false);
    ui_->pushButton_revert->setEnabled(revertEnabled);

    lastOpenedForEditing_ = true;
}



void ProblemDefinitionDialog::prepareNewDefinition()
{
    setWindowTitle("Define a new problem");
    oldName_ = QString();

    ui_->pushButton_saveTemporarily->setEnabled(false);
    ui_->pushButton_saveToFile->setEnabled(false);
    ui_->pushButton_revert->setEnabled(false);

    for (unsigned int i = 0; i < functionLineEdits_.size(); ++i)
    {
        functionLineEdits_.at(i)->clear();
    }
    ui_->lineEdit_problemName->clear();

    lastOpenedForEditing_ = false;
}



void ProblemDefinitionDialog::enableSaveButtons()
{
    ui_->pushButton_saveTemporarily->setEnabled(true);
    if (saveToFileAllowed_)
        ui_->pushButton_saveToFile->setEnabled(true);
}



void ProblemDefinitionDialog::handleSaveButtonClicked() const
{
    bool overwrite = lastOpenedForEditing_;
    bool saveToFile = false;
    if (sender() == ui_->pushButton_saveToFile)
        saveToFile = true;

    QStringList newDefinitions;
    for (unsigned int i = 0; i < functionLineEdits_.size(); ++i)
    {
        QLineEdit* lineEdit = functionLineEdits_.at(i);
        if (lineEdit->isVisible())
            newDefinitions.append(lineEdit->text());
    }

    emit saveDefinitionRequested(overwrite, saveToFile, oldName_,
                                 ui_->lineEdit_problemName->text(),
                                 newDefinitions);
}



void ProblemDefinitionDialog::on_pushButton_cancel_clicked()
{
    close();
}



void ProblemDefinitionDialog::on_pushButton_help_clicked()
{
    if (!infoDefiningDialog_)
        infoDefiningDialog_ = new InfoDefiningDialog(this);

    infoDefiningDialog_->show();
}
