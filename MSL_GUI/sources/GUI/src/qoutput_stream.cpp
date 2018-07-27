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


#include "qoutput_stream.h"


QOutputStream::QOutputStream(std::ostream& stream, QTextEdit* text_edit, QObject* parent)
    : QObject(parent),
      m_stream(stream)
{
    m_old_buf = stream.rdbuf();
    stream.rdbuf(this);

    connect(this, SIGNAL(sendString(QString)), text_edit, SLOT(append(QString)));
}



QOutputStream::~QOutputStream()
{
    // output anything that is left
    if (!m_string.empty())
        emit sendString(m_string.c_str());
    m_stream.rdbuf(m_old_buf);
}



QOutputStream::int_type QOutputStream::overflow(int_type v)
{
    mutex.lock();
    if (v == '\n')
    {
        emit sendString(m_string.c_str());
        m_string.erase(m_string.begin(), m_string.end());
    }
    else
        m_string += v;

    mutex.unlock();
    return v;
}



std::streamsize QOutputStream::xsputn(const char *p, std::streamsize n)
{
    mutex.lock();

    m_string.append(p, p + n);
    int pos = 0;
    while (pos != std::string::npos)
    {
        pos = m_string.find('\n');
        if (pos != std::string::npos)
        {
            std::string tmp(m_string.begin(), m_string.begin() + pos);
            emit sendString(tmp.c_str());
            m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
        }
    }

    mutex.unlock();
    return n;
}
