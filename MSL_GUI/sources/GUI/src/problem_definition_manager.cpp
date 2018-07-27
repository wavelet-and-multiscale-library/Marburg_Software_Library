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


#include <QMessageBox>
#include <QTextStream>
#include <QComboBox>

#include "interfaces/abstract_problemtype_module.h"
#include "problem_definition_dialog.h"

#include "problem_definition_manager.h"


const QString ProblemDefinitionManager::DEFINITIONS_FILE_NAME = QStringLiteral("custom_problem_definitions.xml");



ProblemDefinitionManager::ProblemDefinitionManager(QWidget* parent, QComboBox* problemComboBox)
    : QObject(parent),
      problemComboBox_(problemComboBox),
      definitionsDocument_("custom_problem_definitions"),
      definitionsFile_(DEFINITIONS_FILE_NAME),
      iconLocked_(QIcon(":/icons/locked.svg")),
      iconSaved_(QIcon(":/icons/document-save.svg")),
      iconUnsaved_(QIcon(":/icons/image-loading.svg"))
{
    bool documentOkay;

    if (QFile::exists(DEFINITIONS_FILE_NAME))
    {
        if (!definitionsFile_.open(QIODevice::ReadOnly))
        {
            documentOkay = false;
            QString errorMessage = QString("Could not read file \'%1\'. Custom problem definitions "
                                           "can not be saved!\n\n%2")
                                   .arg(DEFINITIONS_FILE_NAME, definitionsFile_.errorString());
            QMessageBox::critical(parent, "File reading error", errorMessage);
        }
        else
        {
            QString parseError;
            if (!definitionsDocument_.setContent(&definitionsFile_, &parseError))
            {
                documentOkay = false;
                definitionsDocument_ = QDomDocument();
                QString errorMessage = QString("Could not parse content of file \'%1\'. Custom problem "
                                               "definitions can not be saved!\n\n%2")
                        .arg(DEFINITIONS_FILE_NAME, parseError);
                QMessageBox::critical(parent, "Document parsing error", errorMessage);
            }
            else
            {
                documentOkay = true;
            }
            definitionsFile_.close();
        }
    }
    else
    {
        QDomElement root = definitionsDocument_.createElement("custom_problem_definitions");
        definitionsDocument_.appendChild(root);
        documentOkay = true;
    }

    definitionDialog_ = new ProblemDefinitionDialog(parent, documentOkay);

    connect(problemComboBox_, SIGNAL(currentIndexChanged(int)),
            this, SLOT(handleProblemSelection(int)));
    connect(definitionDialog_, &ProblemDefinitionDialog::saveDefinitionRequested,
            this, &ProblemDefinitionManager::handleSaveDefRequest);
    connect(definitionDialog_, &ProblemDefinitionDialog::revertToFileStateRequested,
            this, &ProblemDefinitionManager::handleRevertToFileStateRequest);
}



int ProblemDefinitionManager::currentExampleIndex()
{
    if (currentProblemDefinition_->status == LOCKED_EXAMPLE)
        return problemComboBox_->currentIndex();
    else
        return -1;
}



QStringList ProblemDefinitionManager::getFuncDefStringsForSelectedProblem()
{
    int problemIndex = problemComboBox_->currentIndex();
    if (problemIndex < 0)
        return QStringList();
    else
        return definitions_.at(selectedProblemType_).at(problemIndex).functionDefStrings;
}



void ProblemDefinitionManager::handleProblemTypeSelection(AbstractProblemTypeModule* selectedProblemType)
{
    selectedProblemType_ = selectedProblemType;
    functionNamesXml_ = selectedProblemType->getFunctionNames();
    functionNamesXml_.replaceInStrings(QRegExp("[(,)]"), "_");

    problemComboBox_->clear();

    if (definitions_.count(selectedProblemType) == 0)
    {
        addDefinitionsForSelectedProblemType();
    }

    for (const auto& definition : definitions_.at(selectedProblemType))
    {
        switch (definition.status)
        {
        case LOCKED_EXAMPLE:
            problemComboBox_->addItem(iconLocked_, definition.name);
            break;
        case IN_SYNC_WITH_FILE:
            problemComboBox_->addItem(iconSaved_, definition.name);
            break;
        case CHANGED_VS_FILE:
        case TEMPORARY:
            problemComboBox_->addItem(iconUnsaved_, definition.name);
            break;
        }
    }

    if (problemComboBox_->count() == 0)
        emit noProblemSelected();
    else
        problemComboBox_->setCurrentIndex(0);

    definitionDialog_->adaptTo(selectedProblemType->getFunctionNames());
}



void ProblemDefinitionManager::handleProblemSelection(int problemIndex)
{
    if (problemIndex >= 0)
    {
        currentProblemDefinition_ = &(definitions_[selectedProblemType_][problemIndex]);
        bool isExampleProblem = (currentProblemDefinition_->status == LOCKED_EXAMPLE);

        emit selectedProblemChanged(currentProblemDefinition_->fullFunctionStrings, isExampleProblem);
    }
    else
        currentProblemDefinition_ = nullptr;
}



void ProblemDefinitionManager::handleEditProblemDefinitionClicked()
{
    bool revertEnabled = (currentProblemDefinition_->status == CHANGED_VS_FILE);
    definitionDialog_->prepareEditing(currentProblemDefinition_->name,
                                      currentProblemDefinition_->functionDefStrings,
                                      revertEnabled);
    definitionDialog_->show();
}



void ProblemDefinitionManager::handleDeleteProblemDefinitionClicked()
{
    QString question;

    switch (currentProblemDefinition_->status)
    {
    case CHANGED_VS_FILE:
        handleRevertToFileStateRequest();
        return;
    case IN_SYNC_WITH_FILE:
        question = QStringLiteral("Do you really want to delete problem \'%1\' in file \'%2\'?")
                    .arg(currentProblemDefinition_->name, DEFINITIONS_FILE_NAME);
        break;
    case TEMPORARY:
        question = QStringLiteral("Do you really want to delete the temporary problem \'%1\'?")
                    .arg(currentProblemDefinition_->name);
        break;
    default:
        break;
    }

    int answer = QMessageBox::question(qobject_cast<QWidget*>(parent()), "Delete problem?",
                                       question, QMessageBox::Yes, QMessageBox::No);
    if (answer == QMessageBox::No)
        return;

    if (currentProblemDefinition_->status == IN_SYNC_WITH_FILE)
        deleteCurrentDefinitionFromFile();

    definitions_.at(selectedProblemType_).removeAt(problemComboBox_->currentIndex());
    currentProblemDefinition_ = nullptr;
    problemComboBox_->removeItem(problemComboBox_->currentIndex());
}



void ProblemDefinitionManager::handleAddNewProblemDefinitionClicked()
{
    definitionDialog_->prepareNewDefinition();
    definitionDialog_->show();
}



void ProblemDefinitionManager::handleSaveDefRequest(bool overwrite, bool saveToFile, const QString& oldName,
                                                    const QString& newName, const QStringList& newFuncDefs)
{
    QString errors;
    if (selectedProblemType_->detectsErrorsInProblemDefinition(newName, newFuncDefs, errors))
    {
        QMessageBox::warning(definitionDialog_, "Error in problem definition", errors);
        return;
    }

    if (overwrite && saveToFile && (currentProblemDefinition_->status != TEMPORARY))
    {
        QString question = QStringLiteral("Do you really want to overwrite definition of problem \'%1\' "
                                          "in file \'%2\'?").arg(oldName, DEFINITIONS_FILE_NAME);
        int answer = QMessageBox::question(definitionDialog_, "Overwrite definition in file?",
                                           question, QMessageBox::Yes, QMessageBox::No);
        if (answer == QMessageBox::No)
            return;
    }

    ProblemDefinition newDefinition;
    newDefinition.name = newName;
    newDefinition.functionDefStrings = newFuncDefs;
    const QStringList& funcNames = selectedProblemType_->getFunctionNames();
    for (int i = 0; i < funcNames.size(); ++i)
    {
        newDefinition.fullFunctionStrings.append(funcNames.at(i) + " = " + newFuncDefs.at(i));
    }

    if (saveToFile)
    {
        bool savingSuccessful;
        const QString& domainTag = selectedProblemType_->getDomainTag();
        const QString& problemTypeTag = selectedProblemType_->getProblemTypeTag();

        QDomElement domainElement;
        QDomNodeList domainList = definitionsDocument_.elementsByTagName(domainTag);
        if (domainList.isEmpty())
        {
            domainElement = definitionsDocument_.createElement(domainTag);
            definitionsDocument_.documentElement().appendChild(domainElement);
        }
        else
        {
            domainElement = domainList.at(0).toElement();
        }

        if (domainElement.isNull())
        {
            savingSuccessful = false;
        }
        else
        {
            QDomElement problemTypeElem;
            QDomNodeList problemTypeList = domainElement.elementsByTagName(problemTypeTag);
            if (problemTypeList.isEmpty())
            {
                problemTypeElem = definitionsDocument_.createElement(problemTypeTag);
                domainElement.appendChild(problemTypeElem);
            }
            else
            {
                problemTypeElem = problemTypeList.at(0).toElement();
            }

            if (problemTypeElem.isNull())
            {
                savingSuccessful = false;
            }
            else
            {
                QDomElement newDefElement = definitionsDocument_.createElement("problem_definition");
                newDefElement.setAttribute("name", newName);
                for (int i = 0; i < functionNamesXml_.size(); ++i)
                {
                    newDefElement.setAttribute(functionNamesXml_.at(i), newFuncDefs.at(i));
                }

                if (overwrite && !currentProblemDefinition_->fileRepresentative.isNull())
                        problemTypeElem.replaceChild(newDefElement, currentProblemDefinition_->fileRepresentative);
                else
                    problemTypeElem.appendChild(newDefElement);

                newDefinition.fileRepresentative = newDefElement;

                if (!definitionsFile_.open(QIODevice::WriteOnly))
                    savingSuccessful = false;
                else
                {
                    QTextStream stream(&definitionsFile_);
                    stream << definitionsDocument_.toString();
                    savingSuccessful = true;
                    definitionsFile_.close();
                }
            }
        }

        if (savingSuccessful)
            newDefinition.status = IN_SYNC_WITH_FILE;
        else
        {
            QString message = QStringLiteral("Could not save definition for problem \'%1\' in file \'%2\' "
                                             "(file is possibly corrupted).").arg(newName, DEFINITIONS_FILE_NAME);
            QMessageBox::warning(definitionDialog_, "Saving failed", message);

            if (overwrite)
            {
                if (currentProblemDefinition_->status == TEMPORARY)
                    newDefinition.status = TEMPORARY;
                else
                {
                    newDefinition.status = CHANGED_VS_FILE;
                    if (newDefinition.fileRepresentative.isNull())
                        newDefinition.fileRepresentative = currentProblemDefinition_->fileRepresentative;
                }
            }
            else
                newDefinition.status = TEMPORARY;
        }
    }
    else // if (!saveToFile)
    {
        if (overwrite)
        {
            if (currentProblemDefinition_->status == TEMPORARY)
                newDefinition.status = TEMPORARY;
            else
            {
                newDefinition.status = CHANGED_VS_FILE;
                newDefinition.fileRepresentative = currentProblemDefinition_->fileRepresentative;
            }
        }
        else
            newDefinition.status = TEMPORARY;
    }

    int index;
    if (overwrite)
    {
        *currentProblemDefinition_ = newDefinition;
        index = problemComboBox_->currentIndex();
    }
    else
    {
        definitions_.at(selectedProblemType_).append(newDefinition);
        index = definitions_.at(selectedProblemType_).size() - 1;
    }

    handleProblemTypeSelection(selectedProblemType_);
    problemComboBox_->setCurrentIndex(index);

    definitionDialog_->close();
}



void ProblemDefinitionManager::handleRevertToFileStateRequest()
{
    QString question = QStringLiteral("Do you really want to revert problem \'%1\' to file state from file \'%2\'?")
                .arg(currentProblemDefinition_->name, DEFINITIONS_FILE_NAME);

    int answer = QMessageBox::question(qobject_cast<QWidget*>(parent()), "Revert problem?",
                                       question, QMessageBox::Yes, QMessageBox::No);
    if (answer == QMessageBox::No)
        return;

    currentProblemDefinition_->name = currentProblemDefinition_->fileRepresentative.attribute("name");
    currentProblemDefinition_->fullFunctionStrings.clear();
    currentProblemDefinition_->functionDefStrings.clear();

    const QStringList& functionNames = selectedProblemType_->getFunctionNames();
    for (int i = 0; i < functionNames.size(); ++i)
    {
        QString funcDef = currentProblemDefinition_->fileRepresentative.attribute(functionNamesXml_.at(i));
        currentProblemDefinition_->functionDefStrings.append(funcDef);
        currentProblemDefinition_->fullFunctionStrings.append(functionNames.at(i) + " = " + funcDef);
    }
    currentProblemDefinition_->status = IN_SYNC_WITH_FILE;

    definitionDialog_->close();

    int index = problemComboBox_->currentIndex();
    handleProblemTypeSelection(selectedProblemType_);
    problemComboBox_->setCurrentIndex(index);
}



void ProblemDefinitionManager::addDefinitionsForSelectedProblemType()
{
    QList<ProblemDefinition>& definitionList = definitions_[selectedProblemType_];

    const QStringList& functionNames = selectedProblemType_->getFunctionNames();
    const QStringList& exampleNames = selectedProblemType_->getExampleNames();
    const std::vector<QStringList>& exampleDefs = selectedProblemType_->getExampleDefinitions();

    for (int i = 0; i < exampleNames.size(); ++i)
    {
        ProblemDefinition definition;
        definition.name = exampleNames.at(i);
        definition.status = LOCKED_EXAMPLE;
        definition.functionDefStrings = exampleDefs.at(i);
        for (int j = 0; j < functionNames.size(); ++j)
            definition.fullFunctionStrings.append(functionNames.at(j) + " = " + exampleDefs.at(i).at(j));
        definitionList.append(definition);
    }

    QDomElement problemTypeElement = getDomElementForSelectedProblemType();

    if (problemTypeElement.isNull())
        return;

    QDomNodeList nodeList = problemTypeElement.elementsByTagName("problem_definition");
    for (int i = 0; i < nodeList.size(); ++i)
    {
        QDomElement defElement = nodeList.at(i).toElement();
        ProblemDefinition definition;
        definition.name = defElement.attribute("name");
        definition.fileRepresentative = defElement;
        definition.status = IN_SYNC_WITH_FILE;
        for (int i = 0; i < functionNames.size(); ++i)
        {
            QString funcDef = defElement.attribute(functionNamesXml_.at(i));
            definition.functionDefStrings.append(funcDef);
            definition.fullFunctionStrings.append(functionNames.at(i) + " = " + funcDef);
        }
        definitionList.append(definition);
    }
}



void ProblemDefinitionManager::deleteCurrentDefinitionFromFile()
{
    QDomElement selectedProblemType = getDomElementForSelectedProblemType();
    if (!selectedProblemType.isNull())
    {
        selectedProblemType.removeChild(currentProblemDefinition_->fileRepresentative);

        definitionsFile_.open(QIODevice::WriteOnly);
        QTextStream stream(&definitionsFile_);
        stream << definitionsDocument_.toString();
        definitionsFile_.close();
    }
}



QDomElement ProblemDefinitionManager::getDomElementForSelectedProblemType()
{
    const QString& domainTag = selectedProblemType_->getDomainTag();
    const QString& problemTypeTag = selectedProblemType_->getProblemTypeTag();

    QDomNodeList domainList = definitionsDocument_.elementsByTagName(domainTag);
    if (domainList.isEmpty())
        return QDomElement();

    QDomElement domainElement = domainList.at(0).toElement();
    if (domainElement.isNull())
        return QDomElement();

    QDomNodeList problemTypeList = domainElement.elementsByTagName(problemTypeTag);
    if (problemTypeList.isEmpty())
        return QDomElement();

     return problemTypeList.at(0).toElement();
}
