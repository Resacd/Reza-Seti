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
#include <iostream>
#include "VelocitySet.h"
#include "include.h"
#include "Struct.h"
#include "Macro.h"
class Node
{



public:
	std::vector<double> m_distributions;
	std::vector<double> m_newDistributions;
	std::vector<double> m_eqdistributions;
	std::shared_ptr<Node> m_neighbours[_nLatNodes];
	int x_position;
	int y_position;
	int z_position;
	NodeType NodeType_;
	BoundaryType BoundaryType_;
	double density;
	std::vector<double>velocity;
	double velMag = 0;
	std::vector<double> accumulated_shear_stress;
	double accumulated_shear_stress_magnitude;
	double osi;
	double tawss;
	double wall_shear_stress;
	std::vector<double> normal;
	Node() = default;
	Node(int coord_x, int coord_y, int coord_z, NodeType nt);
	~Node();

	double Equilibrium(std::shared_ptr<VelocitySet> velSet, int dir);
	void Stream(int dir);
	void Initialize(std::shared_ptr<VelocitySet> velSet, int dir);
	virtual double Density();
	virtual std::vector<double> Velocity(std::shared_ptr<VelocitySet> velSet);
	virtual void MacroscopicValue(std::shared_ptr<VelocitySet> velSet);
	double Pressure(std::shared_ptr<VelocitySet> velSet);
	void BounceBack(int dir, int opp_dir);
	virtual void ApplyBc(double& vel, std::vector<double>& normal);
	virtual void ApplyEqBc(double& vel, std::vector<double>& normal);

	struct Hash {
		std::size_t operator()(const std::shared_ptr<Node>& node) const {
			return std::hash<int>()(node->x_position) ^
				(std::hash<int>()(node->y_position) << 1) ^
				(std::hash<int>()(node->z_position) << 2);
		}
	};


	struct Equal {
		bool operator()(const std::shared_ptr<Node>& lhs, const std::shared_ptr<Node>& rhs) const {
			return lhs->x_position == rhs->x_position &&
				lhs->y_position == rhs->y_position &&
				lhs->z_position == rhs->z_position;
		}
	};
};

