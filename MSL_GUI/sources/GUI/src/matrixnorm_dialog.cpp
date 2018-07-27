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


#include "matrixnorm_dialog.h"
#include "ui_matrixnorm_dialog.h"

MatrixNormDialog::MatrixNormDialog(QWidget *parent) :
    QDialog(parent),
    ui_(new Ui::MatrixNormDialog)
{
    ui_->setupUi(this);
}



MatrixNormDialog::~MatrixNormDialog()
{
    delete ui_;
}



void MatrixNormDialog::resetValuesToOne()
{
    ui_->doubleSpinBox_normA->setValue(1.0);
    ui_->doubleSpinBox_normAinv->setValue(1.0);
}



double MatrixNormDialog::getNormA() const
{
    return ui_->doubleSpinBox_normA->value();
}



double MatrixNormDialog::getNormAinv() const
{
    return ui_->doubleSpinBox_normAinv->value();
}



void MatrixNormDialog::on_doubleSpinBox_normA_valueChanged(double value)
{
    ui_->doubleSpinBox_normAinv->setMinimum(1.0/value);
}



void MatrixNormDialog::on_doubleSpinBox_normAinv_valueChanged(double value)
{
    ui_->doubleSpinBox_normA->setMinimum(1.0/value);
}
