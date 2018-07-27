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


#include "mslgui_plugin.h"
#include "GUI/interfaces/abstract_problemtype_module.h"

#include "misc/string_conversion.h"
#include "muParserError.h"


#include "interval/sturm_bvp/problemtype.h"
#include "interval/poisson_bvp/problemtype.h"
#include "interval/biharmonic_bvp/problemtype.h"

#include "cube/elliptic_bvp/problemtype.h"
#include "cube/poisson_bvp/problemtype.h"

#include "ldomain/poisson_bvp/problemtype.h"
#include "ldomain/elliptic_bvp/problemtype.h"
#include "ldomain/biharmonic_bvp/problemtype.h"

#include "recring/poisson_bvp/problemtype.h"
#include "recring/biharmonic_bvp/problemtype.h"

#include "slitdomain/poisson_bvp/problemtype.h"


MslGuiPlugin::MslGuiPlugin() :
    communicator_(),
    widgetManager_()
{

}



MslGuiPlugin::~MslGuiPlugin()
{
    for (const auto& pair : problemTypeData_) {
        for (const auto& problemModule : pair.second) {
            delete problemModule;
        }
    }
}



AbstractGuiCommunicator& MslGuiPlugin::getGuiCommunicator()
{
    return communicator_;
}



AbstractParameterWidgetManager& MslGuiPlugin::getParameterWidgetManager()
{
    return widgetManager_;
}






void MslGuiPlugin::initProblemTypeModules()
{
#define ADD_PROBLEMTYPE(domainName, problemTypeName, ns_domain, ns_problemType) \
    addProblemType(domainName, problemTypeName, #ns_domain, #ns_problemType, new ns_domain::ns_problemType::ProblemType);


    ADD_PROBLEMTYPE("Interval (0,1)", "Poisson BVP", interval, poisson_bvp)
    ADD_PROBLEMTYPE("Interval (0,1)", "Sturm BVP", interval, sturm_bvp)
    ADD_PROBLEMTYPE("Interval (0,1)", "Biharmonic BVP", interval, biharmonic_bvp)
    ADD_PROBLEMTYPE("Cube (0,1)²", "Poisson BVP", cube, poisson_bvp)
    ADD_PROBLEMTYPE("Cube (0,1)²", "Elliptic BVP", cube, elliptic_bvp)
    ADD_PROBLEMTYPE("L-Domain", "Poisson BVP", ldomain, poisson_bvp)
    ADD_PROBLEMTYPE("L-Domain", "Elliptic BVP", ldomain, elliptic_bvp)
    ADD_PROBLEMTYPE("L-Domain", "Biharmonic BVP", ldomain, biharmonic_bvp)
    ADD_PROBLEMTYPE("Rectangular Ring", "Poisson BVP", recring, poisson_bvp)
    ADD_PROBLEMTYPE("Rectangular Ring", "Biharmonic BVP", recring, biharmonic_bvp)
    ADD_PROBLEMTYPE("Slit domain", "Poisson BVP", slitdomain, poisson_bvp)
}



const QStringList& MslGuiPlugin::getDomainNameList() const
{
    return domainNames_;
}



const QStringList& MslGuiPlugin::getProblemTypeListForDomain(const QString& domain) const
{
    return problemTypeNames_.at(domain);
}



AbstractProblemTypeModule* MslGuiPlugin::getProblemTypeModule(const QString& domain,
                                                              int problemTypeIndex) const
{
    const auto& vec = problemTypeData_.at(domain);
    return vec.at(problemTypeIndex);
}



void MslGuiPlugin::addProblemType(const QString& domain, const QString& problemTypeName,
                                  const QString& domainTag, const QString& problemTypeTag,
                                  AbstractProblemTypeModule* problemType)
{
    try
    {
        problemType->initialize(domainTag, problemTypeTag);
    }
    catch(const mu::ParserError& e)
    {
        QString error = QString("Parser related error during initialization of problem module!\n"
                                "Problem type: %1. Domain: %2.\n"
                                "Error details:\n\n%3\n\n"
                                "Problem module not loaded!").arg(problemTypeName, domain,
                                                                  qStringFromMuString(e.GetMsg()));

        emit communicator_.errorOccured(error);
        delete problemType;
        return;
    }
    catch(const std::exception& e)
    {
        QString error = QString("Error during initialization of problem module!\n"
                                "Problem type: %1. Domain: %2.\n"
                                "Error details:\n\n%3\n\n"
                                "Problem module not loaded!").arg(problemTypeName, domain,
                                                                  QString(e.what()));

        emit communicator_.errorOccured(error);
        delete problemType;
        return;
    }
    catch(...)
    {
        QString error = QString("Unknown error during initialization of problem module!\n"
                                "Problem type: %1. Domain: %2.\n\n"
                                "Problem module not loaded!").arg(problemTypeName, domain);

        emit communicator_.errorOccured(error);
        delete problemType;
        return;
    }

    if (!domainNames_.contains(domain)) {
        domainNames_.append(domain);
    }

    problemTypeNames_[domain].append(problemTypeName);
    problemTypeData_[domain].push_back(problemType);
}



#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(MslGuiPlugin, MslGuiPlugin)
#endif // QT_VERSION < 0x050000
