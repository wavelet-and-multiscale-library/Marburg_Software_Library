// Copyright (c) 2015 Softwariness.com
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.

#ifndef _THROW_ASSERT_H
#define _THROW_ASSERT_H

#include <exception>
#include <string>
#include <sstream>
#include <iostream>

/// Exception type for assertion failures
class AssertionFailureException : public std::exception
{
private:
    const char* expression;
    const char* file;
    int line;
    std::string report;

public:
    // Construct an assertion failure exception
    AssertionFailureException(const char* expression, const char* file, int line)
        : expression(expression)
        , file(file)
        , line(line)
    {
        std::stringstream sStream;

        std::string expressionString = expression;

        if (expressionString == "false" || expressionString == "0" || expressionString == "FALSE")
        {
            sStream << "Unreachable code assertion";
        }
        else
        {
            sStream << "Assertion '" << expression << "'";
        }

        sStream << " failed in file '" << file << "' line " << line;
        report = sStream.str();

        std::cerr << report << std::endl;
    }

    // The assertion message
    virtual const char* what() const throw()
    {
        return report.c_str();
    }

    // The expression which was asserted to be true
    const char* Expression() const throw()
    {
        return expression;
    }

    // Source file
    const char* File() const throw()
    {
        return file;
    }

    // Source line
    int Line() const throw()
    {
        return line;
    }

    ~AssertionFailureException() throw()
    {
    }
};



#define _ASSERT_H_DECLS

#define __assert_fail(ASSERTION, FILE, LINE, FUNCTION) throw AssertionFailureException(ASSERTION, FILE, LINE)

#define __assert_perror_fail (ERRNUM, FILE, LINE, FUNCTION) throw AssertionFailureException(#ERRNUM, FILE, LINE)

#define __assert(ASSERTION, FILE, LINE) throw AssertionFailureException(ASSERTION, FILE, LINE)

// Assert that EXPRESSION evaluates to true, otherwise raise AssertionFailureException
//#define assert(EXPRESSION) if(!(EXPRESSION)) { throw AssertionFailureException(#EXPRESSION, __FILE__, __LINE__); }



inline void throwAssertionFailureException(const wchar_t* assertion, const wchar_t* file, unsigned line)
{
    char assertstring[100];
    char filestring[100];
    std::wcstombs(assertstring, assertion, 100);
    std::wcstombs(filestring, file, 100);
    throw AssertionFailureException(assertstring, filestring, line);
}

#define __ASSERT_H_

#define _assert(ASSERTION, FILE, LINE) throw AssertionFailureException(ASSERTION, FILE, LINE)

#define _wassert(ASSERTION, FILE, LINE) throwAssertionFailureException(ASSERTION, FILE, LINE)


#endif // _THROW_ASSERT_H

