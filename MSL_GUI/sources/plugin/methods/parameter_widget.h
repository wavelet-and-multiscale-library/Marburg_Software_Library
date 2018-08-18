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


#ifndef METHODS_PARAMETER_WIDGET_H
#define METHODS_PARAMETER_WIDGET_H

#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QComboBox>
#include <QFormLayout>

namespace methods
{

class ParameterWidget : public QWidget
{
public:
    explicit ParameterWidget(const QString& methodName, QWidget* parent = nullptr);

    void addDoubleParameter(const QString& identifier, double startValue, double minValue, double maxValue, double step = 1.0);
    void addDoubleParameter(const QString& identifier, double startValue, double minValue, double maxValue,
                            const QString& labelText, double step = 1.0);

    void addIntParameter(const QString& identifier, int startValue, int minValue, int maxValue, int step = 1);
    void addIntParameter(const QString& identifier, int startValue, int minValue, int maxValue,
                         const QString& labelText, int step = 1);

    void addBoolParameter(const QString& identifier, bool startValue);
    void addBoolParameter(const QString& identifier, bool startValue, const QString& labelText);

    void addStringParameter(const QString& identifier, const QStringList& stringValues);
    void addStringParameter(const QString& identifier, const QStringList& stringValues,
                            const QString& labelText);

    void setTheLayout();

    double getDoubleParameter(const QString& identifier) const;
    int getIntParameter(const QString& identifier) const;
    bool getBoolParameter(const QString& identifier) const;
    QString getStringParameter(const QString& identifier) const;

    QString getParamString() const;

private:
    std::map<QString, QDoubleSpinBox*> doubleBoxes_;
    std::map<QString, QSpinBox*> intBoxes_;
    std::map<QString, QCheckBox*> boolBoxes_;
    std::map<QString, QComboBox*> stringBoxes_;

    QFormLayout* layoutLeft_;
    QFormLayout* layoutRight_;

    int paramCount;
};


}

#endif // METHODS_PARAMETER_WIDGET_H
