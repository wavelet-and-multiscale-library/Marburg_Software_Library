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


#ifndef GENERIC_PROBLEMTYPE_MODULE_H
#define GENERIC_PROBLEMTYPE_MODULE_H


#include "GUI/interfaces/abstract_problemtype_module.h"
#include "GUI/interfaces/gui_inputdata.h"
#include "raw/raw_problem_module.h"
#include "parsed/parsed_expression_checker.h"
#include "discr/discretization_module_guidata.h"
#include "misc/typelist.h"
#include "misc/convert_basis_to_qstring.h"



template<class RP_MODULE, class ... DISCR_MODULES>
class GenericProblemTypeModule;


template<class ... RP_CONTAINERS, class ... DISCR_MODULES>
class GenericProblemTypeModule<RawProblemModule<RP_CONTAINERS...>, DISCR_MODULES...>
        : public AbstractProblemTypeModule
{

public:

    typedef RawProblemModule<RP_CONTAINERS...> RP_Module;

    static const unsigned int dimension = RP_Module::dimension;
    static const std::size_t n_functions = RP_Module::n_functions;


    GenericProblemTypeModule() : checker_(0), initialized_(false) { }


    ~GenericProblemTypeModule()
    {
        delete checker_;
    }


    void initialize(const QString& domainTag, const QString& problemTypeTag) override;


    const QString& getDomainTag() const
    {
        return domainTag_;
    }


    const QString& getProblemTypeTag() const
    {
        return problemTypeTag_;
    }


    const QString& getProblemTypeDescription() const override
    {
        return problemDescription_;
    }


    const QStringList& getFunctionNames() const override
    {
        return completeFunctionNames_;
    }


    const QStringList& getExampleNames() const override
    {
        return exampleNames_;
    }


    const std::vector<QStringList>& getExampleDefinitions() const
    {
        return exampleDefinitions_;
    }


    const QStringList& getDiscretizationTypeList() const override
    {
        return discretizationTypeNames_;
    }


    const QStringList& get1DBasisList(int discretizationIndex) const override
    {
        return discrModuleGuiDataSets_.at(discretizationIndex).basis1DList;
    }


    const QStringList& getMethodList(int discretizationIndex) const override
    {
        return discrModuleGuiDataSets_.at(discretizationIndex).methodList;
    }


    unsigned int getJmaxStandard(int discretizationIndex) const override
    {
        return discrModuleGuiDataSets_.at(discretizationIndex).jMaxStandard;
    }


    int getEpsilonStandardExponent(int discretizationIndex, int methodIndex) const override
    {
        const std::vector<int>& epsilonStandards =
                discrModuleGuiDataSets_.at(discretizationIndex).epsilonStandardExponents;
        return epsilonStandards.at(methodIndex);
    }


    bool detectsErrorsInProblemDefinition(const QString& problemName,
                                          const QStringList& functionDefs,
                                          QString& out_errors) const override;


    AbstractDiscretizedProblem* createDiscretizedProblem(const GuiInputData& input) override
    {
        return selectDiscrModuleAndCreateProblem<DISCR_MODULES...>(input);
    }



protected:

    QString problemDescription_;
    std::array<QString, n_functions> functionIdentifiers_;
    std::array<QString, dimension> functionVariables_;


    virtual void setupNamesAndDescription() = 0; // set up problemDescription_, functionIdentifiers_ and functionVariables_.
    virtual void setupExamples() = 0; // when implementing this method, use addExampleProblem(...)


    void addExampleProblem(const QString& exampleName,
                           const std::array<QString, n_functions>& funcDefinitions,
                           const typename RP_CONTAINERS::RawProblemType* ... exampleRepresentations)
    {
        exampleNames_.append(exampleName);

        QStringList definitions;
        for (int i = 0; i < n_functions; ++i)
            definitions.append(funcDefinitions.at(i));

        exampleDefinitions_.push_back(definitions);
        rawProblemModule_.addExampleProblem(exampleRepresentations...);
    }

private:

    void setupExpressionChecker()
    {
        try
        {
            checker_ = new ParsedExpressionChecker<dimension>(functionVariables_);
        }
        catch(const mu::ParserError& e)
        {
            QString varNames;
            for (const auto& var : functionVariables_)   {
                varNames.append(var);
                varNames.append(',');
            }
            varNames.chop(1);

            QString error = QString("Error during setup of ParsedExpressionChecker:\n%1\n"
                            "(Probably a problem with the given variable name(s) \'%2\')")
                            .arg(qStringFromMuString(e.GetMsg()), varNames);

            throw mu::ParserError(muStringFromQString(error));
        }
    }


    void setupDiscretizationModulesGuiData()
    {
        (void) std::initializer_list<int>{
               (discrModuleGuiDataSets_.push_back(DISCR_MODULES::getGuiData(dimension)), 0)...};

        (void) std::initializer_list<int>{
               (discretizationTypeNames_.append(DISCR_MODULES::getDiscretizationTypeName()), 0)...};
    }


    void setupCompleteFunctionNames()
    {
        QString varNames("(");
        for (const auto& var : functionVariables_)   {
            varNames.append(var);
            varNames.append(',');
        }
        varNames.chop(1);
        varNames.append(')');

        for (const auto& func : functionIdentifiers_)
            completeFunctionNames_ << (func + varNames);
    }


    template<class FIRST_DISCR_MODULE, class ... OTHER_DISCR_MODULES>
    AbstractDiscretizedProblem* selectDiscrModuleAndCreateProblem(const GuiInputData& input)
    {
        if (FIRST_DISCR_MODULE::getDiscretizationTypeName() == input.discretizationType)
        {
            typename FIRST_DISCR_MODULE::List_1D_Bases basisList;
            return select1DBasisAndCreateProblem<FIRST_DISCR_MODULE>(input, basisList);
        }
        else
        {
            return selectDiscrModuleAndCreateProblem<OTHER_DISCR_MODULES...>(input);
        }
    }


    // Base case for recursion:
    template<bool DUMMY = true>
    AbstractDiscretizedProblem* selectDiscrModuleAndCreateProblem(const GuiInputData& input)
    {
        throw std::logic_error("Discretization type \'" + input.discretizationType.toStdString()
                               + "\' was not found in the list of discretization modules for the "
                                 "chosen problemtype module!");
        return nullptr;
    }


    template<class DMODULE, class FIRST_1D_BASIS, class ... OTHER_1D_BASES>
    AbstractDiscretizedProblem*
    select1DBasisAndCreateProblem(const GuiInputData& input,
                                  const TypeList<FIRST_1D_BASIS, OTHER_1D_BASES...>& basisList)
    {
        (void) basisList;
        if (Convert<FIRST_1D_BASIS>::toQString() == input.basis1D)
        {
            typedef typename DMODULE::template Problem<FIRST_1D_BASIS> ProblemType;
            return createProblem<ProblemType>(input);
        }
        else
        {
            return select1DBasisAndCreateProblem<DMODULE>(input, TypeList<OTHER_1D_BASES...>());
        }
    }


    // Base case for recursion:
    template<class DMODULE>
    AbstractDiscretizedProblem* select1DBasisAndCreateProblem(const GuiInputData& input,
                                                              const TypeList<>& basisList)
    {
        (void) basisList;
        throw std::logic_error("1D basis \'" + input.basis1D.toStdString() + "\' was not found "
                               "in the list of 1D bases for the discretization module of type \'"
                               + input.discretizationType.toStdString() + "\'.");
        return nullptr;
    }


    template<class PROBLEMTYPE>
    AbstractDiscretizedProblem* createProblem(const GuiInputData& input)
    {
        typedef typename PROBLEMTYPE::RawProblemType RPType;
        PROBLEMTYPE* problem = new PROBLEMTYPE;

        if (input.exampleIndex >= 0)
        {
            const RPType& rawProblem = rawProblemModule_.template getExampleProblem<RPType>(input.exampleIndex);
            problem->setupProblem(input, rawProblem);
        }
        else
        {
            std::array<QString, n_functions> funcDefinitions;
            for (int i = 0; i < n_functions; ++i)
                funcDefinitions.at(i) = input.functionDefinitions.at(i);

            const RPType& rawProblem =
            rawProblemModule_.template getParsedProblem<RPType>(funcDefinitions);
            problem->setupProblem(input, rawProblem);
        }

        return problem;
    }


    RP_Module rawProblemModule_;
    ParsedExpressionChecker<dimension>* checker_;

    QString domainTag_, problemTypeTag_;
    QStringList exampleNames_;
    std::vector<QStringList> exampleDefinitions_;

    bool initialized_;
    QStringList completeFunctionNames_;

    QStringList discretizationTypeNames_;
    std::vector<DiscretizationModuleGuiData> discrModuleGuiDataSets_;
};



/*##################################################################################################
    Implementation
##################################################################################################*/



template<class ... RP_CONTAINERS, class ... DISCR_MODULES>
void GenericProblemTypeModule<RawProblemModule<RP_CONTAINERS...>, DISCR_MODULES...>::
initialize(const QString& domainTag, const QString& problemTypeTag)
{
    if(!initialized_)
    {
        domainTag_ = domainTag;
        problemTypeTag_ = problemTypeTag;
        setupNamesAndDescription();
        setupExpressionChecker();
        setupExamples();

        rawProblemModule_.setupParsedProblemAndFunctions(functionVariables_);

        setupDiscretizationModulesGuiData();
        setupCompleteFunctionNames();

        initialized_ = true;
    }
}



template<class ... RP_CONTAINERS, class ... DISCR_MODULES>
bool GenericProblemTypeModule<RawProblemModule<RP_CONTAINERS...>, DISCR_MODULES...>::
detectsErrorsInProblemDefinition(const QString& problemName,
                                      const QStringList& functionDefs,
                                      QString& out_errors) const
{
    out_errors = QString("Definition of problem \'%1\' contains "
                         "the following error(s):\n\n").arg(problemName);

    if (functionDefs.size() != n_functions) {
        QString error = QString("Wrong number of function definitions provided (expected %1, "
                                "provided %2).").arg(n_functions).arg(functionDefs.size());
        out_errors.append(error);

        return true;
    }

    bool errorsDetected = false;

    if (problemName.isEmpty())
    {
        out_errors.append("Problem name is empty.\n\n");
        errorsDetected = true;
    }

    for (int i = 0; i < n_functions; ++i) {
        QString syntaxErrors;
        if (checker_->detectsSyntaxErrors(functionDefs.at(i), syntaxErrors)) {
            QString error = QString("In definition of %1 := %2 :\n%3\n")
                        .arg(completeFunctionNames_.at(i), functionDefs.at(i), syntaxErrors);
            out_errors.append(error);
            errorsDetected = true;
        }
    }

    return errorsDetected;
}



#endif // GENERIC_PROBLEMTYPE_MODULE_H
