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


#ifndef GENERIC_DISCRETIZED_PROBLEM_H
#define GENERIC_DISCRETIZED_PROBLEM_H

#include <memory>

#include "misc/typelist.h"
#include "abstract_discretized_problem.h"
#include "solution/abstract_solution.h"
#include "solution/solution_tools.h"
#include "GUI/interfaces/gui_inputdata.h"
#include "GUI/interfaces/abstract_gui_communicator.h"
#include "muParserError.h"

#include "misc/string_conversion.h"
#include "misc/gui_convergence_logger.h"
#include "misc/basis_or_frame.h"



template< template<class P> class CACHED, class EQUATION,
          class RAW_PROBLEM, class SOLUTION, class METHOD_LIST >
class GenericDiscretizedProblem;



template< template<class P> class CACHED, class EQUATION,
          class RAW_PROBLEM, class SOLUTION, class ... METHODS >
class GenericDiscretizedProblem<CACHED, EQUATION, RAW_PROBLEM, SOLUTION, TypeList<METHODS...> >
        : public AbstractDiscretizedProblem
{

public:

    typedef RAW_PROBLEM RawProblemType;

    // Define inner type 'WaveletSystem' to be EQUATION::WaveletBasis if the latter exists,
    // otherwise define 'WaveletSystem' to be EQUATION::Frame
    typedef typename BasisOrFrame<EQUATION>::systemType WaveletSystem;

    typedef typename SOLUTION::CoeffVector CoeffVector;

    static const int space_dimension = EQUATION::space_dimension;


    AbstractSolution* computeSolution(const QString& method, double epsilon, int jmax,
                                      const AbstractGuiCommunicator* communicator) override;

    bool canBeReusedFor(const GuiInputData& newInput) const override;

    void updateTo(const GuiInputData& newInput) override;

    double norm_A() override
    {
        return cachedProblem_->norm_A();
    }

    double norm_Ainv() override
    {
        return cachedProblem_->norm_Ainv();
    }

    void setupProblem(const GuiInputData& input, const RAW_PROBLEM& rawProblem);


protected:

    virtual WaveletSystem* createWaveletSystem(int jmax) = 0;

    virtual EQUATION* createDiscretizedEquation(const RAW_PROBLEM& rawProblem,
                                                const WaveletSystem& basis) = 0;

    GuiInputData lastInput_;

    std::shared_ptr<WaveletSystem> wsystem_;
    std::unique_ptr<EQUATION> uncachedProblem_;
    std::unique_ptr< CACHED<EQUATION> > cachedProblem_;


private:

    template<class FIRST_METHOD, class ... OTHER_METHODS>
    void computeSolutionHelper(CoeffVector& coeffs,
                               MathTL::AbstractConvergenceLogger& logger,
                               const TypeList<FIRST_METHOD, OTHER_METHODS...>& methodList)
    {
        (void) methodList;

        if (FIRST_METHOD::name() == lastInput_.method) {
            FIRST_METHOD::solve(*cachedProblem_, lastInput_, coeffs, logger);
        }
        else {
            computeSolutionHelper(coeffs, logger, TypeList<OTHER_METHODS...>());
        }
    }

    // Base case for recursion:
    void computeSolutionHelper(CoeffVector& coeffs,
                               const MathTL::AbstractConvergenceLogger& logger,
                               const TypeList<>& methodList)
    {
        throw std::logic_error("Method \'" + lastInput_.method.toStdString()
                               + "\' was not found in the method list of the "
                                 "discretized problem class!");
        (void) coeffs;
        (void) logger;
        (void) methodList;
    }
};



/*##################################################################################################
    Implementation
##################################################################################################*/


template< template<class P> class CACHED, class EQUATION, class RAW_PROBLEM, class SOLUTION, class ... METHODS >
AbstractSolution*
GenericDiscretizedProblem<CACHED, EQUATION, RAW_PROBLEM, SOLUTION, TypeList<METHODS...> >::
computeSolution(const QString& method, double epsilon, int jmax, const AbstractGuiCommunicator* communicator)
{
    lastInput_.method = method;
    lastInput_.epsilon = epsilon;
    lastInput_.jmax = jmax;

    CoeffVector coeffs;
    bool aborted = false;
    bool errorOccured = false;

    GuiConvergenceLogger* logger = new GuiConvergenceLogger(communicator);

    try
    {
        solution_tools::adjustCoeffVector(coeffs, *wsystem_);
        computeSolutionHelper(coeffs, *logger, TypeList<METHODS...>());
    }
    catch(const abort_request& e)
    {
        (void) e;
        aborted = true;
    }
    catch(const mu::ParserError& e)
    {
        QString errorMessage("Parser related error during computation of solution coefficients:\n\n");
        errorMessage.append(qStringFromMuString(e.GetMsg()));
        emit communicator->errorOccured(errorMessage);
        errorOccured = true;
    }
    catch(const std::exception& theException)
    {
        QString errorMessage("Error during computation of solution coefficients:\n\n");
        errorMessage.append(QString(theException.what()));
        emit communicator->errorOccured(errorMessage);
        errorOccured = true;
    }
    catch(...)
    {
        QString errorMessage("Unknown error during computation of solution coefficients.");
        emit communicator->errorOccured(errorMessage);
        errorOccured = true;
    }

    SOLUTION* solution = new SOLUTION(lastInput_, coeffs, wsystem_, logger);  // from here on logger is owned by solution

    try
    {
        solution_tools::scaleCoefficients(coeffs, *cachedProblem_, -1);
    }
    catch(...)
    {
        QString errorMessage("Error: scaling of solution coefficients failed!");
        emit communicator->errorOccured(errorMessage);
        errorOccured = true;
    }

    solution->setScaledCoefficients(coeffs);

    if (aborted)
        solution->setIncompleteDueToAbort(true);
    if (errorOccured)
        solution->setIncompleteDueToError(true);

    return solution;
}



template< template<class P> class CACHED, class EQUATION, class RAW_PROBLEM, class SOLUTION, class ... METHODS >
bool GenericDiscretizedProblem<CACHED, EQUATION, RAW_PROBLEM, SOLUTION, TypeList<METHODS...> >::
canBeReusedFor(const GuiInputData& newInput) const
{
    bool canReuse = (newInput.domain == lastInput_.domain)
                    && (newInput.problemType == lastInput_.problemType)
                    && (newInput.discretizationTypeIndex == lastInput_.discretizationTypeIndex)
                    && (newInput.basis1D == lastInput_.basis1D)
                    && (newInput.jmax == lastInput_.jmax)
                    && (newInput.pmax == lastInput_.pmax)
                    && (newInput.overlap == lastInput_.overlap)
                    && (((newInput.exampleIndex == lastInput_.exampleIndex) && (newInput.exampleIndex >= 0)) ||
                        ((newInput.exampleIndex < 0) && (lastInput_.exampleIndex < 0) && (newInput.functionDefinitions == lastInput_.functionDefinitions)));

    return canReuse;
}



template< template<class P> class CACHED, class EQUATION, class RAW_PROBLEM, class SOLUTION, class ... METHODS >
void GenericDiscretizedProblem<CACHED, EQUATION, RAW_PROBLEM, SOLUTION, TypeList<METHODS...> >::
updateTo(const GuiInputData& newInput)
{
    if (newInput.normEstimatesProvided)
    {
        cachedProblem_->set_normA(newInput.normA);
        cachedProblem_->set_normAinv(newInput.normAinv);
    }
    else
    {
        if (lastInput_.normEstimatesProvided)
        {
            cachedProblem_->set_normA(0.0);
            cachedProblem_->set_normAinv(0.0);
        }
    }

    lastInput_ = newInput;
}



template< template<class P> class CACHED, class EQUATION, class RAW_PROBLEM, class SOLUTION, class ... METHODS >
void GenericDiscretizedProblem<CACHED, EQUATION, RAW_PROBLEM, SOLUTION, TypeList<METHODS...> >::
setupProblem(const GuiInputData& input, const RAW_PROBLEM& rawProblem)
{
    lastInput_ = input;

    wsystem_.reset(createWaveletSystem(input.jmax));

    EQUATION* equation = createDiscretizedEquation(rawProblem, *wsystem_);

    uncachedProblem_.reset(equation);

    if (input.normEstimatesProvided)
    {
        cachedProblem_.reset(new CACHED<EQUATION>(equation, input.normA, input.normAinv));
    }
    else
    {
        cachedProblem_.reset(new CACHED<EQUATION>(equation));
    }
}


#endif // GENERIC_DISCRETIZED_PROBLEM_H
