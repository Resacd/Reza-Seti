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



#include "Node.h"

Node::Node(int coord_x, int coord_y, int coord_z, NodeType nt)
	: x_position(coord_x), y_position(coord_y), z_position(coord_z), NodeType_(nt)
{
	NodeType_ = nt;
	m_eqdistributions.resize(_nLatNodes);
        velocity.resize(_dim);
}


Node::~Node()
{
}

void Node::Stream(int dir)
{

   m_neighbours[dir]->m_distributions[dir] = m_newDistributions[dir];
}

double Node::Equilibrium(std::shared_ptr<VelocitySet> velSet, int dir)
{
	double du = Velocity(velSet)[0] * velSet->GetDirection(dir)[0] + Velocity(velSet)[1] * velSet->GetDirection(dir)[1]+ Velocity(velSet)[2] * velSet->GetDirection(dir)[2];
	double u_sqr = Velocity(velSet)[0] * Velocity(velSet)[0] + Velocity(velSet)[1] * Velocity(velSet)[1] + Velocity(velSet)[2] * Velocity(velSet)[2];

	return velSet->GetWeight(dir) * Density() * (1 + 3 * du + 9 / 2 * du * du - 3 / 2 * u_sqr);
}

void Node::Initialize(std::shared_ptr<VelocitySet> velSet, int dir)
{
	double rho = 1.0;   

	m_distributions.push_back(velSet->GetWeight(dir));
	m_newDistributions.push_back(0.0);
	velocity.resize(3);
}

double Node::Density()
{
	double density = 0;
	for (int dir = 0; dir < _nLatNodes; dir++)
	{
		density += this->m_distributions[dir];	
	}
	return density;
}



std::vector<double> Node::Velocity(std::shared_ptr<VelocitySet> velSet)
{

	for (int dir = 0; dir < _nLatNodes; dir++)
	{
		for (int dim = 0; dim < 3; dim++)
		{
			velocity[dim] += velSet->GetDirection(dir)[dim] * this->m_distributions[dir];
		}
	}

	for (int dim = 0; dim < 3; dim++)
	{
		velocity[dim] /= Density();
	}

	return velocity;
}



void Node::MacroscopicValue(std::shared_ptr<VelocitySet> velSet) {

	density = 0;
	velocity[0] = 0.0;
	velocity[1] = 0.0;
	velocity[2] = 0.0;

	for (int dir = 0; dir < _nLatNodes; dir++) {

		density += m_distributions[dir];
		velocity[0] += m_distributions[dir] * velSet->GetDirection(dir)[0];
		velocity[1] += m_distributions[dir] * velSet->GetDirection(dir)[1];
		velocity[2] += m_distributions[dir] * velSet->GetDirection(dir)[2];

	}

	velocity[0] /= density;
	velocity[1] /= density;
	velocity[2] /= density;
}


double Node::Pressure(std::shared_ptr<VelocitySet> velSet)
{
	double density = Density();
	double c_square = velSet->Get_c_square(0);
	double pressure = c_square * density;
	return pressure;
}


void Node::BounceBack(int dir, int opp_dir)
{

	m_distributions[opp_dir] = m_newDistributions[dir];
}


void Node::ApplyBc(double& vel, std::vector<double>& normal)
{

}


void Node::ApplyEqBc(double& vel, std::vector<double>& normal)
{

}
