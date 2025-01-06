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



#include "VelocitySet.h"

VelocitySet::VelocitySet(std::vector<double> weights, std::vector<std::vector<int>> directions,
	std::vector<std::pair<std::vector<int>, std::vector<double>>> vectorOfPairs ,double cSquare)
	: m_weights(weights), m_directions(directions), m_UnkDisNor(vectorOfPairs), m_cSquare(cSquare), m_nrOfDimensions(directions[0].size())
{
}

VelocitySet::~VelocitySet()
{
}

double VelocitySet::GetWeight(int dir) const
{
	return m_weights[dir];
}

std::vector<int> VelocitySet::GetDirection(int dir) const
{
	return m_directions[dir];
}

double VelocitySet::Get_c_square(int dir) const
{
	return m_cSquare;
}

int VelocitySet::Get_nDimensions() const
{
	return m_nrOfDimensions;
}

int VelocitySet::Get_nDirections() const
{
	return m_directions.size();
}

int VelocitySet::OppositeDirection(int dir) const
{
	return 0;
}
