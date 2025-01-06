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


#pragma once

#include "VelocitySet.h"
#include "Subdomain.h"
#include "OutputFile.h"
#include "Boundary.h"
#include <vector>
#include "pulsatile.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <filesystem>
#include <cstdlib> 



class Simulation
{

public:
	//OutputFile file1;
	int myrank_, num_procs_;
	SubDomain SubDomain_;
	Pulsatile Pulsatile_;
	double pressure_;
	double m_relaxation;
        double m_omega;
	double m_dt;
        double total_runtime;
        double Cs;
       	SimulationParameters simParams;  
	std::vector<double> wall_shear_stress;
	double InletVelocity_;
        Simulation ()=default;
	Simulation(SubDomain& SubDomain, SimulationParameters &simParams);
	~Simulation();
	void InitializeShearStressAccumulators();
	void Collision();
	void Streaming();
	void Communicate();
	void ApplyZouHeBc();
	void ApplyInBc();
	void ApplyOutletBc();
	void ApplyWallBc();
	void ComputeWallShearStress(int iteration);
	void ComputeOSIAndTAWSS(int iteration);
	double getInletVelocity(double x, double period);
	void Iterate();
	void SetRelaxation(double tau);
	void WriteOutput(int iter, const std::string& phase_name);


};


