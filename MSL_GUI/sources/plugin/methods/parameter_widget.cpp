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


#include "parameter_widget.h"
#include <QLabel>

namespace methods
{

ParameterWidget::ParameterWidget(const QString& methodName, QWidget* parent)
    : QWidget(parent), layoutLeft_(new QFormLayout), layoutRight_(new QFormLayout), paramCount(0)
{
    layoutLeft_->setLabelAlignment(Qt::AlignRight);
    layoutRight_->setLabelAlignment(Qt::AlignRight);

    //QString widgetTitle = QString("Additional parameters for algorithm %1").arg(methodName);
    //setTitle(widgetTitle);
}



void ParameterWidget::addDoubleParameter(const QString& identifier, double startValue,
                                         double minValue, double maxValue, double step)
{
    addDoubleParameter(identifier, startValue, minValue, maxValue, identifier, step);
}



void ParameterWidget::addDoubleParameter(const QString& identifier, double startValue,
                                         double minValue, double maxValue, const QString& labelText, double step)
{
    QDoubleSpinBox* doubleSpinBox = new QDoubleSpinBox;
    doubleSpinBox->setRange(minValue, maxValue);
    doubleSpinBox->setValue(startValue);
    doubleSpinBox->setSingleStep(step);
    doubleSpinBox->setMaximumWidth(80);

    doubleBoxes_[identifier] = doubleSpinBox;

    if (paramCount%2 == 0)
        layoutLeft_->addRow(labelText + " =", doubleSpinBox);
    else
        layoutRight_->addRow(labelText + " =", doubleSpinBox);

    paramCount++;
}



void ParameterWidget::addIntParameter(const QString& identifier, int startValue,
                                      int minValue, int maxValue, int step)
{
    addIntParameter(identifier, startValue, minValue, maxValue, identifier, step);
}



void ParameterWidget::addIntParameter(const QString& identifier, int startValue,
                                      int minValue, int maxValue, const QString& labelText, int step)
{
    QSpinBox* spinBox = new QSpinBox;
    spinBox->setRange(minValue, maxValue);
    spinBox->setValue(startValue);
    spinBox->setSingleStep(step);
    spinBox->setMaximumWidth(80);

    intBoxes_[identifier] = spinBox;

    if (paramCount%2 == 0)
        layoutLeft_->addRow(labelText + " =", spinBox);
    else
        layoutRight_->addRow(labelText + " =", spinBox);

    paramCount++;
}



void ParameterWidget::addBoolParameter(const QString& identifier, bool startValue)
{
    addBoolParameter(identifier, startValue, identifier);
}



void ParameterWidget::addBoolParameter(const QString& identifier, bool startValue,
                                       const QString& labelText)
{
    QCheckBox* checkBox = new QCheckBox(labelText);
    checkBox->setChecked(startValue);

    boolBoxes_[identifier] = checkBox;

    if (paramCount%2 == 0)
        layoutLeft_->addRow(checkBox);
    else
        layoutRight_->addRow(checkBox);

    paramCount++;
}



void ParameterWidget::addStringParameter(const QString& identifier, const QStringList& stringValues)
{
    addStringParameter(identifier, stringValues, identifier);
}



void ParameterWidget::addStringParameter(const QString& identifier, const QStringList& stringValues,
                                         const QString& labelText)
{
    QComboBox* comboBox = new QComboBox;
    comboBox->addItems(stringValues);

    stringBoxes_[identifier] = comboBox;

    if (paramCount%2 == 0)
        layoutLeft_->addRow(labelText + ":", comboBox);
    else
        layoutRight_->addRow(labelText + ":", comboBox);

    paramCount++;
}



void ParameterWidget::setTheLayout()
{
    QHBoxLayout* horizontalLayout = new QHBoxLayout();
    horizontalLayout->addLayout(layoutLeft_);
    horizontalLayout->addSpacing(15);
    horizontalLayout->addLayout(layoutRight_);
    horizontalLayout->addStretch();
    setLayout(horizontalLayout);
}



double ParameterWidget::getDoubleParameter(const QString& identifier) const
{
    return doubleBoxes_.at(identifier)->value();
}



int ParameterWidget::getIntParameter(const QString& identifier) const
{
    return intBoxes_.at(identifier)->value();
}



bool ParameterWidget::getBoolParameter(const QString& identifier) const
{
    return boolBoxes_.at(identifier)->isChecked();
}



QString ParameterWidget::getStringParameter(const QString& identifier) const
{
    return stringBoxes_.at(identifier)->currentText();
}



QString ParameterWidget::getParamString() const
{
    QStringList parameters;

    for (int i = 0; i < layoutLeft_->rowCount(); i++)
    {
        {
            QWidget* field = layoutLeft_->itemAt(i, QFormLayout::FieldRole)->widget();

            QCheckBox* checkbox = qobject_cast<QCheckBox*>(field);
            if (checkbox)
            {
                parameters << checkbox->text() + ((checkbox->isChecked())? ": Yes" : ": No");
            }
            else
            {
                QLabel* label = qobject_cast<QLabel*>(layoutLeft_->itemAt(i, QFormLayout::LabelRole)->widget());
                QString param = label->text() + " %1";

                QAbstractSpinBox* spinbox = qobject_cast<QAbstractSpinBox*>(field);
                if (spinbox)
                {
                    parameters << param.arg(spinbox->text());
                }
                else
                {
                    QComboBox* combobox = qobject_cast<QComboBox*>(field);
                    if (combobox)
                    {
                        parameters << param.arg(combobox->currentText());
                    }
                }
            }
        }

        if (i < layoutRight_->rowCount())
        {
            QWidget* field = layoutRight_->itemAt(i, QFormLayout::FieldRole)->widget();

            QCheckBox* checkbox = qobject_cast<QCheckBox*>(field);
            if (checkbox)
            {
                parameters << checkbox->text() + ((checkbox->isChecked())? ": Yes" : ": No");
            }
            else
            {
                QLabel* label = qobject_cast<QLabel*>(layoutRight_->itemAt(i, QFormLayout::LabelRole)->widget());
                QString param = label->text() + " %1";

                QAbstractSpinBox* spinbox = qobject_cast<QAbstractSpinBox*>(field);
                if (spinbox)
                {
                    parameters << param.arg(spinbox->text());
                }
                else
                {
                    QComboBox* combobox = qobject_cast<QComboBox*>(field);
                    if (combobox)
                    {
                        parameters << param.arg(combobox->currentText());
                    }
                }
            }
        }
    }

    return parameters.join(" *** ");
}


}
