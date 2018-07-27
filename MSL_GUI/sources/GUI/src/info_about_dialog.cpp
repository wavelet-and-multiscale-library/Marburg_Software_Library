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


#include <QDesktopServices>
#include "info_about_dialog.h"
#include "ui_info_about_dialog.h"

InfoAboutDialog::InfoAboutDialog(QWidget* parent) :
    QDialog(parent),
    ui_(new Ui::InfoAboutDialog)
{
    ui_->setupUi(this);
}

InfoAboutDialog::~InfoAboutDialog()
{
    delete ui_;
}

void InfoAboutDialog::on_label_credits_linkActivated(const QString &link)
{
    if (link == "GPLv3")
        ui_->tabWidget->setCurrentWidget(ui_->tab_license);
    else
        QDesktopServices::openUrl(link);
}

void InfoAboutDialog::on_label_warranty_linkActivated(const QString &link)
{
    if (link == "GPLv3")
        ui_->tabWidget->setCurrentWidget(ui_->tab_license);
}
