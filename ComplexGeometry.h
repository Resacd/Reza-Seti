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
#ifndef _MESH_GENERATOR_H_
#define _MESH_GENERATOR_H_

#include "Struct.h"
#include <tuple>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include "Struct.h"
#include <unordered_set>
#include <mpi.h>
#include <omp.h>
#include <map>
#include "GeometryWriter.h"
#include "GeometryReader.h"
#include "mpi_call.h"
#include "Utility.h"
#include "SimulationParameters.h"

class Geometry {

public:

	Writer Writer;
	Reader Reader;
	Utility Utility;
	MPI_CALL<MPI_Datatype> mpi_caller;
	std::vector<std::vector<std::vector<Point>>> pointsBox;
	std::vector<Triangle>globalTriangles;
	std::vector<Triangle>localTriangles;
	BoundingBox boundingbox;
	std::vector<Point> boundaryPoints;
	std::pair<std::vector<Point>, std::vector<double>> InletboundaryPoints; 
	std::vector<std::pair<std::vector<Point>, std::vector<double>>> OutletboundaryPoints;
	std::vector<Point> WallboundaryPoints;

	SimulationParameters simParams;  
        std::string stlFile;
	int globalNumtris;
	std::vector<int> GeomSize;
	double spacing;
	double rayLength;
	int myrank, num_procs;
        int nx,ny,nz;
	std::shared_ptr<VelocitySet> velSet_;

	Geometry() = default;
	Geometry(std::shared_ptr<VelocitySet> velSet, SimulationParameters &simParams);
	~Geometry();
	void GetBoundingbox();
	void DistributeTriangles();
	void generatePointsInBoundingBox();
	bool doLinesIntersect(const Line& line1, const Line& line2);
        bool isPointInsideTriangle(double x, double y, double z, const Triangle& triangle);
        int countIntersectinsTriangle(Line& line, std::vector<Triangle>& posintersectedtri);  
	std::vector<Point> findBoundaryPoints();
	void adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch);
	void labelPoints();
	std::vector<Point> excludeInletOutlet();
	void getInOutBP();
	void writeBoundaryPoints();
};

#endif
