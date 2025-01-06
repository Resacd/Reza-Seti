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
#include "Node.h"
#include <fstream>
#include <string>
#include <mpi.h>
#include "Utility.h"
#include <filesystem>
#include <memory> 
#include <map>
#include <vector>
#include <tuple>
#include <cmath>
#include <filesystem>
#include <sys/stat.h>
#include <limits>




#ifndef _WIN32
#include <sys/stat.h>
#else
#include <windows.h>
#endif



#include "Struct.h"
#include "HaloLayer.h"
#include "Subdomain.h"

#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDelaunay3D.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkTriangleFilter.h>
#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangle.h>
#include <vtkPolygon.h>
#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>
#include <vtkContourFilter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPCANormalEstimation.h>
#include <vtkFillHolesFilter.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>




class OutputFile
{
private:

	std::mutex m_fileMutex;
	int m_outputIdx;
	std::mutex m_dataMutex;
	Utility Utility_;

public:

	OutputFile() = default;
	OutputFile(int outputIdx);
	~OutputFile();

    double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);
    double calculateMeanNeighboringValues(const std::vector<DataPoint>& data, const DataPoint& point, double searchRadius, char attribute);
    void smoothData(std::vector<DataPoint>& data, double neighborSearchRadius, const std::string& flowType);
	void printElapsedTime(std::chrono::steady_clock::time_point startTime);
	void processFlowData(const std::string& inputFileName, const std::string& outputFileName, double neighborSearchRadius, const std::string& flowType);
	void WriteVTK_LocalData(SubDomain &subdomain);
	void WriteVTK_Collective(SubDomain& subdomain);
	void WriteSurfaceAttribute(SubDomain &subdomain,int iter);

};


