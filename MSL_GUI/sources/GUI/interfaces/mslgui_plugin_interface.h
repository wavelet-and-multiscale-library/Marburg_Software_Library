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


#ifndef MSLGUI_PLUGIN_INTERFACE_H
#define MSLGUI_PLUGIN_INTERFACE_H

#include <QtPlugin>
#include <QStringList>

class AbstractGuiCommunicator;
class AbstractProblemTypeModule;
class AbstractParameterWidgetManager;


class MslGuiPluginInterface
{
public:
    virtual ~MslGuiPluginInterface() {}

    virtual AbstractGuiCommunicator& getGuiCommunicator() = 0;
    virtual AbstractParameterWidgetManager& getParameterWidgetManager() = 0;

    virtual void initProblemTypeModules() = 0;

    virtual const QStringList& getDomainNameList() const = 0;
    virtual const QStringList& getProblemTypeListForDomain(const QString& domain) const = 0;

    virtual AbstractProblemTypeModule* getProblemTypeModule(const QString& domain,
                                                            int problemTypeIndex) const = 0;
};


#define MslGuiPluginInterface_iid "MslGui.PluginInterface/1.0"

Q_DECLARE_INTERFACE(MslGuiPluginInterface, MslGuiPluginInterface_iid)

#endif // MSLGUI_PLUGIN_INTERFACE_H
