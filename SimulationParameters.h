/*
**********************************************************************************************
Author: Reza Bozorgpour, PhD Student"
Biomedical Engineering Group, University of Wisconsin-Milwaukee
Original Code: Developed from scratch with no prior code provided
Programmer Contact Information:
 
 - Email: bozorgp2@uwm.edu
Date: 2022-10-14
For any modifications to the software, developers are strongly advised to contact the owner.
**********************************************************************************************#
Copyright (c) 2022 by Mahsa Dabagh
All rights reserved

 Redistribution and use in source and binary forms, with or without modification, are
 permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of
    conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of
    conditions and the following disclaimer in the documentation and/or other materials
    provided with the distribution.

 THIS CODE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************
*/



#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include <vector>
#include <string>
#include <tuple>
#include "Struct.h"

struct SimulationParameters {

    std::string stlFile;            // Input file name
    std::string FlowType;            // Pulsatile or Contant
    double dt=0.0;                  // Time resolution
    double rho = 0.0;               // Density
    double nu_blood = 0.0;          // Kinematic viscosity
    double tau = 0.0;               // Relaxation factor
    double dx = 0.0;                // Grid size
    double pressure = 0.0;          // Pressure
    double reynoldsNumber = 0.0;    // Reynolds number
    int num_Out=0;                  // Number of outlets
    double CteVel = 0.0;            // Inlet Constant Velocity
    double Systole = 0;             // Systole time
    double Diastole=0;              // Diastole time
    int NumCycle=0;                 // Number of Cardiac Cycle
    int MaxIter=0;                  // Maximum Iteration
    int IterToWrite=0;              // Interval to Write 
    std::vector<Equation> equations;  

    struct Inlet {
        std::string plane;                    // Plane (x, y, or z)
        double coordinate = 0.0;              // Coordinate of the plane
        std::vector<double> range_x;          // Range in y-direction
        std::vector<double> range_y;          // Range in y-direction
        std::vector<double> range_z;          // Range in z-direction
        std::vector<double> normal;           // Normal vector
    };

    struct Outlet {
        std::string plane;              // Plane (x, y, or z)
        double coordinate;              // Coordinate of the plane
        std::vector<double> range_x;    // Range in y-direction
        std::vector<double> range_y;    // Range in y-direction
        std::vector<double> range_z;    // Range in z-direction
        std::vector<double> normal;     // Normal vector
    };

    std::vector<Inlet> inlets;          // List of inlets
    std::vector<Outlet> outlets;        // List of outlets
};

#endif // SIMULATION_PARAMETERS_H

