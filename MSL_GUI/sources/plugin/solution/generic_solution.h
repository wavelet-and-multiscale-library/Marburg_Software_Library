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


#ifndef GENERIC_SOLUTION_H
#define GENERIC_SOLUTION_H

#include <memory>

#include "abstract_solution.h"
#include "GUI/interfaces/gui_inputdata.h"
#include "solution_tools.h"
#include "misc/gui_convergence_logger.h"
#include "main/gui_communicator.h"



template<class COEFF_VECTOR, class PLOT_DATA, class WSYSTEM>
class GenericSolution : public AbstractSolution
{
public:

    typedef COEFF_VECTOR CoeffVector;

    GenericSolution(const GuiInputData& input,
                    const CoeffVector& coeffs,
                    std::shared_ptr<WSYSTEM> wsystem,
                    GuiConvergenceLogger* logger)
        : input_(input), coeffs_(coeffs), wsystem_(wsystem), scaledCoeffs_(), plotData_(),
          logger_(logger), optionalDataLogNames_(logger->getOptionalDataLogNames()),
          optionalTimeLogNames_(logger->getOptionalTimeLogNames()), isIncompleteDueToAbort_(false),
          isIncompleteDueToError_(false)
    {

    }


    void setScaledCoefficients(const CoeffVector& scaledCoeffs)
    {
        scaledCoeffs_ = scaledCoeffs;
    }


    void computePlotData(int resolution) override
    {
        //plotData_ = WaveletTL::evaluate(*wsystem_, scaledCoeffs_, true, resolution);
        solution_tools::evaluate(*wsystem_, scaledCoeffs_, resolution, plotData_);
    }


    void sendPlotDataToGuiVia(const GuiCommunicator& communicator) const override
    {
        communicator.sendPlotDataToGui(plotData_, input_.computationNumber);
    }


    void writeSolutionPlotToMatlabFile(const char* filename) const override
    {
        std::ofstream outStream(filename);

        MathTL::matlab_output(outStream, plotData_);

        outStream << "title('Solution plot')" << std::endl;
        writeMatlabInputSummary(outStream);
        outStream.close();
    }


    void writeIndexPlotToMatlabFile(const char* filename) const override
    {
        std::ofstream outStream(filename);
        outStream << "figure\n" << std::endl;

        solution_tools::plotIndices(wsystem_.get(), coeffs_, input_.jmax, outStream);

        writeMatlabInputSummary(outStream);
        outStream.close();
    }


    void writeConvergenceLogsToMatlabFile(const char* filename) const override
    {
        std::ofstream outStream(filename);
        outStream << "figure\n" << std::endl;

        logger_->writeConvergencePlots(outStream);

        writeMatlabInputSummary(outStream);
        outStream.close();
    }


    void writeOptionalLogToMatlabFile(int logIndex, const char* filename) const override
    {
        std::ofstream outStream(filename);
        outStream << "figure\n" << std::endl;

        std::string logName;

        if (logIndex < optionalDataLogNames_.size())
        {
            logName = optionalDataLogNames_.at(logIndex);
            logger_->writeOptionalDataPlot(logName, outStream);
        }
        else
        {
            logName = optionalTimeLogNames_.at(logIndex - optionalDataLogNames_.size());
            logger_->writeOptionalTimePlot(logName, outStream);
        }

        outStream << "title('" << logName << "')" << std::endl;
        writeMatlabInputSummary(outStream);
        outStream.close();
    }


    bool isIncompleteDueToAbort() const override
    {
        return isIncompleteDueToAbort_;
    }


    bool isIncompleteDueToError() const override
    {
        return isIncompleteDueToError_;
    }


    void setIncompleteDueToAbort(bool incomplete)
    {
        isIncompleteDueToAbort_ = incomplete;
    }


    void setIncompleteDueToError(bool incomplete)
    {
        isIncompleteDueToError_ = incomplete;
    }


    QStringList getOptionalLogNames() const override
    {
        QStringList namesList;

        for (const auto& name : optionalDataLogNames_)
        {
            namesList.append(QString::fromStdString(name));
        }

        for (const auto& name : optionalTimeLogNames_)
        {
            namesList.append(QString::fromStdString(name));
        }

        return namesList;
    }


protected:
    const GuiInputData input_;
    const CoeffVector coeffs_;
    const std::shared_ptr<WSYSTEM> wsystem_;
    CoeffVector scaledCoeffs_;
    PLOT_DATA plotData_;

    const std::unique_ptr<GuiConvergenceLogger> logger_;
    const std::vector<std::string> optionalDataLogNames_;
    const std::vector<std::string> optionalTimeLogNames_;

    bool isIncompleteDueToAbort_;
    bool isIncompleteDueToError_;


    void writeMatlabInputSummary(std::ofstream& outStream) const
    {
        outStream << "summary = " << input_.matlabSummary.toStdString() << ";" << std::endl;
        outStream << "annotation('textbox', [0 0 1 0.1], 'String', summary, 'FontWeight', 'bold', 'FontSize', 6, 'EdgeColor', 'none', 'VerticalAlignment', 'bottom')" << std::endl;

        if (isIncompleteDueToAbort())
        {
            outStream << "str = 'Computation incomplete (aborted)';" << std::endl;
            outStream << "annotation('textbox', [0 0.9 1 0.1], 'String', str, 'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'right')" << std::endl;
        }
        if (isIncompleteDueToError())
        {
            outStream << "str = 'Computation incomplete (due to error)';" << std::endl;
            outStream << "annotation('textbox', [0 0.9 1 0.1], 'String', str, 'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'right')" << std::endl;
        }
    }
};




#endif // GENERIC_SOLUTION_H
