/*
**********************************************************************************************
Author: Reza Bozorgpour, PhD Student
Biomedical Engineering Group, University of Wisconsin-Milwaukee
Original Code: Developed from scratch with no prior code provided
Programmer Contact Information:
 
 - Email: bozorgp2@uwm.edu
Date: 2022-10-14

For any modifications to the software, developers are strongly advised to contact the owner.
**********************************************************************************************
Copyright (c) 2022 by Mahsa Dabagh
All rights reserved

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of
   conditions and the following disclaimer in the documentation and/or other materials
   provided with the distribution.

THIS CODE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************
*/



#pragma once

#ifndef _PALSATILE_H_
#define _PALSATILE_H_

#include <cmath>

class Pulsatile {

public:

    double period = 0.9;

    inline double InletVelequation1(double x) {
        return (-8.579e+08 * pow(x, 8) + 4.059e+08 * pow(x, 7) - 7.817e+07 * pow(x, 6) +
            7.994e+06 * pow(x, 5) - 4.835e+05 * pow(x, 4) + 1.782e+04 * pow(x, 3) -
            321.6 * pow(x, 2) + 3.215 * x + 0.1423);
    }

    inline double InletVelequation2(double x) {
        return (-4.32e+06 * pow(x, 8) + 7.94e+06 * pow(x, 7) - 6.298e+06 * pow(x, 6) +
            2.82e+06 * pow(x, 5) - 7.809e+05 * pow(x, 4) + 1.37e+05 * pow(x, 3) -
            1.484e+04 * pow(x, 2) + 901.9 * x - 22.87);
    }

    inline double InletVelequation3(double x) {
        return (-5.123 * pow(x, 4) + 14.15 * pow(x, 3) - 14.21 * pow(x, 2) + 5.865 * x - 0.5696);
    }

};

#endif
