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


#ifndef PROBLEM_DEFINITION_MANAGER_H
#define PROBLEM_DEFINITION_MANAGER_H

#include <QObject>
#include <QFile>
#include <QIcon>
#include <QtXml/QDomDocument>

class ProblemDefinitionDialog;
class QComboBox;
class AbstractProblemTypeModule;


enum DefinitionStatus
{
    LOCKED_EXAMPLE,
    IN_SYNC_WITH_FILE,
    CHANGED_VS_FILE,
    TEMPORARY
};


struct ProblemDefinition
{
   QString name;
   QStringList functionDefStrings;
   QStringList fullFunctionStrings;
   DefinitionStatus status;
   QDomElement fileRepresentative;
};



class ProblemDefinitionManager : public QObject
{
    Q_OBJECT
public:
    ProblemDefinitionManager(QWidget* parent, QComboBox* problemComboBox);

    int currentExampleIndex();
    QStringList getFuncDefStringsForSelectedProblem();

    static const QString DEFINITIONS_FILE_NAME;

public slots:
    void handleProblemTypeSelection(AbstractProblemTypeModule* selectedProblemType);
    void handleProblemSelection(int problemIndex);

    void handleEditProblemDefinitionClicked();
    void handleDeleteProblemDefinitionClicked();
    void handleAddNewProblemDefinitionClicked();

private slots:
    void handleSaveDefRequest(bool overwrite, bool saveToFile, const QString& oldName,
                              const QString& newName, const QStringList& newFuncDefs);
    void handleRevertToFileStateRequest();

signals:
    void selectedProblemChanged(const QStringList& fullFunctionStrings, bool isExampleProblem);
    void noProblemSelected();

private:
    void addDefinitionsForSelectedProblemType();
    void deleteCurrentDefinitionFromFile();
    QDomElement getDomElementForSelectedProblemType();


    std::map< AbstractProblemTypeModule*, QList<ProblemDefinition> > definitions_;
    AbstractProblemTypeModule* selectedProblemType_;
    ProblemDefinition* currentProblemDefinition_;

    QComboBox* problemComboBox_;
    QDomDocument definitionsDocument_;
    QFile definitionsFile_;
    ProblemDefinitionDialog* definitionDialog_;

    QStringList functionNamesXml_;

    const QIcon iconLocked_, iconSaved_, iconUnsaved_;
};

#endif // PROBLEM_DEFINITION_MANAGER_H
