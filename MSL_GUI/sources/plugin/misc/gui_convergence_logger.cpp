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


#include "gui_convergence_logger.h"
#include "GUI/interfaces/abstract_gui_communicator.h"


GuiConvergenceLogger::GuiConvergenceLogger(const AbstractGuiCommunicator* communicator)
    : MathTL::ConvergenceLogger(true, true, true, true),
      communicator_(communicator)
{

}



void GuiConvergenceLogger::logConvergenceData(double degrees_of_freedom, double approx_error)
{
    double elapsedSeconds = getElapsedSeconds();
    convergenceLog_.addData(degrees_of_freedom, approx_error);
    convergenceTimeLog_.addData(elapsedSeconds, approx_error);

    lastApproxError_ = approx_error;
    iterations_++;

    emit communicator_->convergencePlotDataPairComputed(degrees_of_freedom, approx_error);
    emit communicator_->convergenceTimePlotDataPairComputed(elapsedSeconds, approx_error);
}



void GuiConvergenceLogger::checkAbortConditions()
{
    if (communicator_->gotAbortRequest())
    {
        logMessage("\nComputation aborted by request!\n");
        throw abort_request();
    }

    MathTL::ConvergenceLogger::checkAbortConditions();
}



const MathTL::DataLog& GuiConvergenceLogger::getOptionalDataLog(const QString& logName) const
{
    return MathTL::ConvergenceLogger::getOptionalDataLog(logName.toStdString());
}



const MathTL::TimeLog& GuiConvergenceLogger::getOptionalTimeLog(const QString& logName) const
{
    return MathTL::ConvergenceLogger::getOptionalTimeLog(logName.toStdString());
}
