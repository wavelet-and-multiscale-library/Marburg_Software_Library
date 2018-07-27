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


#ifndef QOUTPUT_STREAM_H
#define QOUTPUT_STREAM_H


#include <QObject>
#include <QTextEdit>
#include <QMutex>


class QOutputStream : public QObject, std::basic_streambuf<char>
{
   Q_OBJECT

signals:
   void sendString(QString text);

public:
   QOutputStream(std::ostream& stream, QTextEdit* text_edit, QObject* parent = nullptr);

   ~QOutputStream();

protected:
   virtual int_type overflow(int_type v);

   virtual std::streamsize xsputn(const char* p, std::streamsize n);

private:
   std::ostream& m_stream;
   std::streambuf* m_old_buf;
   std::string m_string;
   QMutex mutex;
};

#endif // QOUTPUT_STREAM_H
