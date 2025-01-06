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




#include "Boundary.h"
#include <cmath>
#include <iostream>

Boundary::Boundary(int coord_x, int coord_y, int coord_z, std::shared_ptr <VelocitySet> velset, NodeType nt)
    : Node(coord_x, coord_y, coord_z, nt) 
{

    m_inletVel.resize(_dim);
    velSet = velset;

}

Boundary::~Boundary() {


}

double Boundary::Density() {

	return 1;
}


void Boundary::SetInletVel(double ux, double uy, double uz)
{
    m_inletVel.push_back(ux);
    m_inletVel.push_back(uy);
    m_inletVel.push_back(uz);
}

std::vector<double> Boundary::Velocity(std::shared_ptr<VelocitySet> velSet)
{
    if (BoundaryType_ == BoundaryType::Inlet) {
        return m_inletVel;
    }
    if (BoundaryType_ == BoundaryType::Wall) {

        return { 0.0,0.0,0.0 };
    }
}



void Boundary::MacroscopicValue(std::shared_ptr<VelocitySet> velSet) {



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

void Boundary::ApplyBc(double& vel, std::vector<double> &normal)
{


    if (normal[0] > 0) {

        m_inletVel[0] = vel;
        m_inletVel[1] = 0.0;
        m_inletVel[2] = 0.0;

        double rhoN = (m_distributions[2] + m_distributions[3] + m_distributions[4] + m_distributions[5] +
            m_distributions[14] + m_distributions[15] + m_distributions[16] + m_distributions[17] + m_distributions[18] +
            2 * (m_distributions[0] + m_distributions[6] + m_distributions[7] + m_distributions[8] + m_distributions[9])) / (1 + m_inletVel[0]);

        double Nxy = (1.0 / 2) * (m_distributions[2] + m_distributions[14] + m_distributions[15]
            - m_distributions[3] - m_distributions[16] - m_distributions[17]) - (1.0 / 3) * (rhoN * m_inletVel[1]);

        double Nxz = (1.0 / 2) * (m_distributions[4] + m_distributions[16] + m_distributions[14]
            - m_distributions[5] - m_distributions[15] - m_distributions[17]) - (1.0 / 3) * (rhoN * m_inletVel[2]);


        m_distributions[1]  = m_distributions[0] - (1 / 3.0) * rhoN * m_inletVel[0];

        m_distributions[10] = m_distributions[7] + (1 / 6.0) * rhoN * (-m_inletVel[0] + m_inletVel[1]) - Nxy;

        m_distributions[11] = m_distributions[6] + (1 / 6.0) * rhoN * (-m_inletVel[0] - m_inletVel[1]) + Nxy;

        m_distributions[13] = m_distributions[8] + (1 / 6.0) * rhoN * (-m_inletVel[0] - m_inletVel[2]) + Nxz;

        m_distributions[12] = m_distributions[9] + (1 / 6.0) * rhoN * (-m_inletVel[0] + m_inletVel[2]) - Nxz;

    }


    if (normal[0] < 0) {

        m_inletVel[0] = vel;
        m_inletVel[1] = 0.0;
        m_inletVel[2] = 0.0;

        double rhoN = (m_distributions[2] + m_distributions[3] + m_distributions[4] + m_distributions[5] +
            m_distributions[14] + m_distributions[15] + m_distributions[16] + m_distributions[17] + m_distributions[18] +
            2 * (m_distributions[1] + m_distributions[10] + m_distributions[11] + m_distributions[12] + m_distributions[13])) / (1 - m_inletVel[0]);

        double Nxy = (1.0 / 2) * (m_distributions[2] + m_distributions[14] + m_distributions[15]
            - m_distributions[3] - m_distributions[16] - m_distributions[17]) - (1.0 / 3) * (rhoN * m_inletVel[1]);

        double Nxz = (1.0 / 2) * (m_distributions[4] + m_distributions[16] + m_distributions[14]
            - m_distributions[5] - m_distributions[15] - m_distributions[17]) - (1.0 / 3) * (rhoN * m_inletVel[2]);


        m_distributions[0] = m_distributions[1] + (1 / 3.0) * rhoN * m_inletVel[0];

        m_distributions[7] = m_distributions[10] + (1 / 6.0) * rhoN * (m_inletVel[0] - m_inletVel[1]) + Nxy;

        m_distributions[6] = m_distributions[11] + (1 / 6.0) * rhoN * (m_inletVel[0] + m_inletVel[1]) - Nxy;

        m_distributions[8] = m_distributions[13] + (1 / 6.0) * rhoN * (m_inletVel[0] + m_inletVel[2]) - Nxz;

        m_distributions[9] = m_distributions[12] + (1 / 6.0) * rhoN * (m_inletVel[0] - m_inletVel[2]) + Nxz;

    }


    if (normal[1] > 0) {

        m_inletVel[0] = 0.0;
        m_inletVel[1] = 1.0;
        m_inletVel[2] = 0.0;

        double rhoN = (m_distributions[0] + m_distributions[1] + m_distributions[4] + m_distributions[5] +
            m_distributions[8] + m_distributions[9] + m_distributions[12] + m_distributions[13] + m_distributions[18] +
            2 * (m_distributions[2] + m_distributions[6] + m_distributions[10] + m_distributions[14] + m_distributions[15])) / (1 + m_inletVel[1]);

        double Nyx = (1.0 / 2) * (m_distributions[0] + m_distributions[8] + m_distributions[9]
            - m_distributions[1] - m_distributions[12] - m_distributions[13]) - (1.0 / 3) * (rhoN * m_inletVel[0]);

        double Nyz = (1.0 / 2) * (m_distributions[4] + m_distributions[8] + m_distributions[12]
            - m_distributions[5] - m_distributions[9] - m_distributions[13]) - (1.0 / 3) * (rhoN * m_inletVel[2]);


        m_distributions[3]  = m_distributions[2]  - (1 / 3.0) * rhoN *  m_inletVel[1];

        m_distributions[11] = m_distributions[6]  + (1 / 6.0) * rhoN * (-m_inletVel[0] - m_inletVel[1]) + Nyx;

        m_distributions[7]  = m_distributions[10] + (1 / 6.0) * rhoN * (-m_inletVel[1] + m_inletVel[0]) - Nyx;

        m_distributions[17] = m_distributions[14] + (1 / 6.0) * rhoN * (-m_inletVel[1] - m_inletVel[2]) + Nyz;

        m_distributions[16] = m_distributions[15] + (1 / 6.0) * rhoN * (-m_inletVel[1] + m_inletVel[2]) - Nyz;
    }

    if (normal[1] < 0) {

        m_inletVel[0] = 0.0;
        m_inletVel[1] = 1.0;
        m_inletVel[2] = 0.0;

        double rhoN = (m_distributions[0] + m_distributions[1] + m_distributions[4] + m_distributions[5] +
            m_distributions[8] + m_distributions[9] + m_distributions[12] + m_distributions[13] + m_distributions[18] +
            2 * (m_distributions[3] + m_distributions[7] + m_distributions[11] + m_distributions[16] + m_distributions[17])) / (1 - m_inletVel[1]);

        double Nyx = (1.0 / 2) * (m_distributions[0] + m_distributions[8] + m_distributions[9]
            - m_distributions[1] - m_distributions[12] - m_distributions[13]) - (1.0 / 3) * (rhoN * m_inletVel[0]);

        double Nyz = (1.0 / 2) * (m_distributions[4] + m_distributions[8] + m_distributions[12]
            - m_distributions[5] - m_distributions[9] - m_distributions[13]) - (1.0 / 3) * (rhoN * m_inletVel[2]);


        m_distributions[2] = m_distributions[3]   + (1 / 3.0) * rhoN *  m_inletVel[1];

        m_distributions[6] = m_distributions[11]  + (1 / 6.0) * rhoN * ( m_inletVel[0] + m_inletVel[1]) - Nyx;

        m_distributions[10] = m_distributions[7]  + (1 / 6.0) * rhoN * ( m_inletVel[1] - m_inletVel[0]) + Nyx;

        m_distributions[14] = m_distributions[17] + (1 / 6.0) * rhoN * ( m_inletVel[1] + m_inletVel[2]) - Nyz;

        m_distributions[15] = m_distributions[16] + (1 / 6.0) * rhoN * ( m_inletVel[1] - m_inletVel[2]) + Nyz;
    }



    if (normal[2] > 0) {

        m_inletVel[0] = 0.0;
        m_inletVel[1] = 0.0;
        m_inletVel[2] = vel;

        double rhoN = (m_distributions[0] + m_distributions[1] + m_distributions[4] + m_distributions[5] +
            m_distributions[8] + m_distributions[9] + m_distributions[12] + m_distributions[13] + m_distributions[18] +
            2 * (m_distributions[3] + m_distributions[7] + m_distributions[11] + m_distributions[16] + m_distributions[17])) / (1 + m_inletVel[2]);

        double Nzx = (1.0 / 2) * ( m_distributions[0]  + m_distributions[6]   + m_distributions[7]
            - m_distributions[1] - m_distributions[10] - m_distributions[11]) - (1.0 / 3) * (rhoN * m_inletVel[0]);

        double Nzy = (1.0 / 2) * (m_distributions[2] + m_distributions[6] + m_distributions[10]
            - m_distributions[3] - m_distributions[7] - m_distributions[11]) - (1.0 / 3) * (rhoN * m_inletVel[1]);


        m_distributions[5] = m_distributions[4] - (1 / 3.0)   * rhoN *   m_inletVel[2];

        m_distributions[9] = m_distributions[12] + (1 / 6.0)  * rhoN * ( m_inletVel[0] - m_inletVel[2]) - Nzx;

        m_distributions[13] = m_distributions[8] + (1 / 6.0)  * rhoN * (-m_inletVel[0] - m_inletVel[2]) + Nzx;

        m_distributions[15] = m_distributions[16] + (1 / 6.0) * rhoN * ( m_inletVel[1] - m_inletVel[2]) - Nzy;

        m_distributions[17] = m_distributions[14] + (1 / 6.0) * rhoN * (-m_inletVel[1] - m_inletVel[2]) + Nzy;
    }


    if (normal[2]<0) {

        m_inletVel[0] = 0.0;
        m_inletVel[1] = 0.0;
        m_inletVel[2] = vel;

        double rhoN = (m_distributions[0] + m_distributions[1] + m_distributions[2] + m_distributions[3] +
            m_distributions[6] + m_distributions[7] + m_distributions[10] + m_distributions[11] + m_distributions[18] +
            2 * (m_distributions[5] + m_distributions[9] + m_distributions[13] + m_distributions[15] + m_distributions[17])) / (1 - m_inletVel[2]);

        double Nzx = (1.0 / 2) * (m_distributions[0] + m_distributions[6] + m_distributions[7]
            - m_distributions[1] - m_distributions[10] - m_distributions[11]) - (1.0 / 3) * (rhoN * m_inletVel[0]);

        double Nzy = (1.0 / 2) * (m_distributions[2] + m_distributions[6] + m_distributions[10]
            - m_distributions[3] - m_distributions[7] - m_distributions[11]) - (1.0 / 3) * (rhoN * m_inletVel[1]);


        m_distributions[4] = m_distributions[5] + (1 / 3.0) * rhoN * m_inletVel[2];

        m_distributions[8] = m_distributions[13] + (1 / 6.0) * rhoN * (m_inletVel[0] + m_inletVel[2]) - Nzx;

        m_distributions[12] = m_distributions[9] + (1 / 6.0) * rhoN * (-m_inletVel[0] + m_inletVel[2]) + Nzx;

        m_distributions[14] = m_distributions[17] + (1 / 6.0) * rhoN * (m_inletVel[1] + m_inletVel[2]) - Nzy;

        m_distributions[16] = m_distributions[15] + (1 / 6.0) * rhoN * (-m_inletVel[1] + m_inletVel[2]) + Nzy;
    }
}




void Boundary::ApplyEqBc(double& vel, std::vector<double> &norma)
{

        m_inletVel[0] = 0.0;
        m_inletVel[1] = 0.0;
        m_inletVel[2] = vel;

    for (int dir = 0; dir < _nLatNodes; dir++) {
        double rho = 1.0;
        double du = m_inletVel[0] * velSet->GetDirection(dir)[0] +
            m_inletVel[1] * velSet->GetDirection(dir)[1] +
            m_inletVel[2] * velSet->GetDirection(dir)[2];
        double u_sqr = m_inletVel[0] * m_inletVel[0] + m_inletVel[1] * m_inletVel[1] + m_inletVel[2] * m_inletVel[2];
        double equilibriumDist = velSet->GetWeight(dir) * rho *
            (1 + 3 * du + 9.0 / 2.0 * du * du - 3.0 / 2.0 * u_sqr);

        m_distributions[dir] = equilibriumDist;
    }
}
