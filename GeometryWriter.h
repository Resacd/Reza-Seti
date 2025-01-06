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
#pragma once
#ifndef _WRITER_H_
#define _WRITER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <array>
#include "Struct.h"
#include <unordered_map>
#include "Node.h"

class Writer {

public:

    Writer();
    ~Writer();
    void writePointoffToFile(const std::vector<Point>& points, const std::string& filename);
    void writePointsToFile(const std::vector<Point>& pointsBox, const std::string& filename);
    void writeBoxToFile(const std::vector<std::vector<std::vector<Point>>>& points, const std::string& filename);
    void writePointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename);
    void writeLabeledPointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename);
    void writeLineToVTK(const Line& line, const std::string& filename);
    void WriteLatticeToVTK(const std::string& filename, const std::vector<std::shared_ptr<Node>>& lattice, int start_x,
        int start_y, int start_z, double spacing);
    void WriteNotNullLatToVTKFile(const std::string& filename, std::vector<std::vector<std::vector<std::shared_ptr<Node>>>>& lattice);
};	

#endif 
