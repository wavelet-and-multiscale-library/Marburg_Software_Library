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


#ifndef MSLGUI_PLUGIN_H
#define MSLGUI_PLUGIN_H

#include <QObject>
#include <QtPlugin>

#include "GUI/interfaces/mslgui_plugin_interface.h"
#include "gui_communicator.h"
#include "methods/parameter_widgetmanager.h"



class MslGuiPlugin : public QObject, public MslGuiPluginInterface
{
    Q_OBJECT
#if QT_VERSION >= 0x050000
    Q_PLUGIN_METADATA(IID MslGuiPluginInterface_iid FILE "MslGuiPlugin.json")
#endif // QT_VERSION >= 0x050000
    Q_INTERFACES(MslGuiPluginInterface)

public:
    MslGuiPlugin();
    ~MslGuiPlugin() override;

    AbstractGuiCommunicator& getGuiCommunicator() override;
    AbstractParameterWidgetManager& getParameterWidgetManager() override;

    void initProblemTypeModules() override;

    const QStringList& getDomainNameList() const override;
    const QStringList& getProblemTypeListForDomain(const QString& domain) const override;

    AbstractProblemTypeModule* getProblemTypeModule(const QString& domain, int problemTypeIndex)
    const override;

private:
    void addProblemType(const QString& domain, const QString& problemTypeName,
                        const QString& domainTag, const QString& problemTypeTag,
                        AbstractProblemTypeModule* problemType);

    GuiCommunicator communicator_;
    methods::ParameterWidgetManager widgetManager_;

    QStringList domainNames_;
    std::map<QString, QStringList> problemTypeNames_;
    std::map<QString, std::vector<AbstractProblemTypeModule*> > problemTypeData_;

};

#endif // MSLGUI_PLUGIN_H
