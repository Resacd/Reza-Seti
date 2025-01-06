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




#include "Utility.h"
#include <algorithm>
#include <iostream>
#include <fstream>


Utility::Utility() {


}

Utility::~Utility() {


}

bool Utility::isPointInBoundary(Point& targetPoint, const std::vector<Point>& boundaryPoints) {

    for (const Point& boundaryPoint : boundaryPoints) {

        const double epsilon = 1e-6; 

        if (std::abs(targetPoint.x - boundaryPoint.x) < epsilon &&
            std::abs(targetPoint.y - boundaryPoint.y) < epsilon &&
            std::abs(targetPoint.z - boundaryPoint.z) < epsilon) {
            return true;
        }
    }


    return false;
}


bool Utility::AreTrianglesEqual(const Triangle& a, const Triangle& b) {
    return (
        a.x1 == b.x1 && a.y1 == b.y1 && a.z1 == b.z1 &&
        a.x2 == b.x2 && a.y2 == b.y2 && a.z2 == b.z2 &&
        a.x3 == b.x3 && a.y3 == b.y3 && a.z3 == b.z3 &&
        a.nx == b.nx && a.ny == b.ny && a.nz == b.nz
        );
}

void Utility::RemoveDuplicateTriangles(std::vector<Triangle>& triangles) {
    triangles.erase(std::unique(triangles.begin(), triangles.end(), AreTrianglesEqual), triangles.end());
}


