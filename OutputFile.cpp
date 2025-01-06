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



#include "OutputFile.h"




OutputFile::OutputFile(int outputIdx)
    : m_outputIdx(outputIdx)
{
}

OutputFile::~OutputFile()
{
}





void OutputFile::WriteVTK_Collective(SubDomain& subdomain)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string folder_name = "Fluid_Variables";

#ifdef _WIN32
    std::filesystem::create_directory(folder_name);
#else
    mkdir(folder_name.c_str(), 0777);
#endif

    std::stringstream file_name_stream;
    file_name_stream << folder_name << "/" << std::to_string(m_outputIdx) << ".vtu";
    std::string file_name = file_name_stream.str();

    size_t local_size = 0;
    for (int z = 0; z < subdomain.my_Nz_; ++z) {
        for (int y = 0; y < subdomain.my_Ny_; ++y) {
            for (int x = 0; x < subdomain.my_Nx_; ++x) {
                auto& node = subdomain.m_lattice[x][y][z];
                if (node != nullptr && node->NodeType_ != NodeType::Halo) {
                    local_size++;
                }
            }
        }
    }

    std::vector<double> local_data_send(local_size * 8);
    size_t index = 0;

    for (int z = 0; z < subdomain.my_Nz_; ++z) {
        for (int y = 0; y < subdomain.my_Ny_; ++y) {
            for (int x = 0; x < subdomain.my_Nx_; ++x) {
                auto& node = subdomain.m_lattice[x][y][z];
                if (node != nullptr && node->NodeType_ != NodeType::Halo) {
                    double velMag = sqrt(node->velocity[0] * node->velocity[0] +
                        node->velocity[1] * node->velocity[1] +
                        node->velocity[2] * node->velocity[2]);
                    std::vector<int> globInD = subdomain.Local2Global(x, y, z);
                    local_data_send[index++] = static_cast<double>(globInD[0]);
                    local_data_send[index++] = static_cast<double>(globInD[1]);
                    local_data_send[index++] = static_cast<double>(globInD[2]);
                    local_data_send[index++] = node->velocity[0];
                    local_data_send[index++] = node->velocity[1];
                    local_data_send[index++] = node->velocity[2];
                    local_data_send[index++] = node->density;
                    local_data_send[index++] = velMag;
                }
            }
        }
    }

    std::vector<int> local_sizes(size);
    MPI_Gather(&local_size, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displacements(size, 0);
    std::vector<int> recv_counts(size);
    if (rank == 0) {
        for (int i = 1; i < size; ++i) {
            displacements[i] = displacements[i - 1] + local_sizes[i - 1] * 8;
        }

        for (int i = 0; i < size; ++i) {
            recv_counts[i] = local_sizes[i] * 8;
        }
    }

    int total_size = std::accumulate(local_sizes.begin(), local_sizes.end(), 0) * 8;

    std::vector<double> all_data(total_size);

    MPI_Gatherv(local_data_send.data(), local_size * 8, MPI_DOUBLE,
        all_data.data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, comm);

    if (rank == 0) {
        vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkDoubleArray> velocityArray = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> densityArray = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> velMagArray = vtkSmartPointer<vtkDoubleArray>::New();

        velocityArray->SetNumberOfComponents(3);
        velocityArray->SetName("Velocity");

        densityArray->SetNumberOfComponents(1);
        densityArray->SetName("Density");

        velMagArray->SetNumberOfComponents(1);
        velMagArray->SetName("VelocityMagnitude");

        std::map<std::vector<int>, vtkIdType> pointIDMap;
        for (size_t i = 0; i < total_size / 8; ++i) {
            size_t dataIndex = i * 8;

            double x = all_data[dataIndex];
            double y = all_data[dataIndex + 1];
            double z = all_data[dataIndex + 2];

            vtkIdType pid = points->InsertNextPoint(x, y, z);
            std::vector<int> globalIndex = { static_cast<int>(x), static_cast<int>(y), static_cast<int>(z) };
            pointIDMap[globalIndex] = pid;

            double velocity[3] = { all_data[dataIndex + 3], all_data[dataIndex + 4], all_data[dataIndex + 5] };
            velocityArray->InsertNextTuple(velocity);

            double density = all_data[dataIndex + 6];
            densityArray->InsertNextValue(density);

            double velMag = all_data[dataIndex + 7];
            velMagArray->InsertNextValue(velMag);
        }

        int global_Nx = subdomain.Geometry_.nx;
        int global_Ny = subdomain.Geometry_.ny;
        int global_Nz = subdomain.Geometry_.nz;

        for (int z = 0; z < global_Nz - 1; ++z) {
            for (int y = 0; y < global_Ny - 1; ++y) {
                for (int x = 0; x < global_Nx - 1; ++x) {
                    std::vector<vtkIdType> cellPoints(8);

                    std::vector<int> indices[8] = {
                        {x, y, z},
                        {x + 1, y, z},
                        {x + 1, y + 1, z},
                        {x, y + 1, z},
                        {x, y, z + 1},
                        {x + 1, y, z + 1},
                        {x + 1, y + 1, z + 1},
                        {x, y + 1, z + 1}
                    };

                    bool validCell = true;
                    for (int i = 0; i < 8; ++i) {

                        if (pointIDMap.find(indices[i]) == pointIDMap.end()) {
                            validCell = false;
                            break;
                        }
                        cellPoints[i] = pointIDMap[indices[i]];
                    }

                    if (validCell) {
                        vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
                        for (int i = 0; i < 8; ++i) {
                            hex->GetPointIds()->SetId(i, cellPoints[i]);
                        }
                        uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
                    }
                }
            }
        }



        uGrid->SetPoints(points);
        uGrid->GetPointData()->AddArray(velocityArray);
        uGrid->GetPointData()->AddArray(densityArray);
        uGrid->GetPointData()->AddArray(velMagArray);

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(file_name.c_str());
        writer->SetInputData(uGrid);
        writer->Write();
    }
}



void OutputFile::WriteSurfaceAttribute(SubDomain& subdomain, int iter) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string file_name_prefix = "_FlowHemodynamics.csv";
    std::string file_name = std::to_string(iter) + file_name_prefix;

    size_t local_size = subdomain.nextToTheWall.size();
    std::vector<double> local_data_send;

    bool isPulsatile = (subdomain.Geometry_.simParams.FlowType == "Pulsatile");
    size_t data_elements = isPulsatile ? 7 : 5;
    local_data_send.resize(local_size * data_elements);

    size_t index = 0;
    for (const auto& node : subdomain.nextToTheWall) {
        std::vector<int> ret = subdomain.Local2Global(node->x_position, node->y_position, node->z_position);
        double x = subdomain.Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].x;
        double y = subdomain.Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].y;
        double z = subdomain.Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].z;

        local_data_send[index++] = x;
        local_data_send[index++] = y;
        local_data_send[index++] = z;
        local_data_send[index++] = node->wall_shear_stress;

        if (isPulsatile) {
            local_data_send[index++] = node->osi;
            local_data_send[index++] = node->tawss;
        }
        local_data_send[index++] = node->velMag; 
    }

    std::vector<int> local_sizes(size);
    MPI_Gather(&local_size, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displacements(size, 0);
    std::vector<int> recv_counts(size);
    if (rank == 0) {
        for (int i = 1; i < size; ++i) {
            displacements[i] = displacements[i - 1] + local_sizes[i - 1] * data_elements;
        }

        for (int i = 0; i < size; ++i) {
            recv_counts[i] = local_sizes[i] * data_elements;
        }
    }

    int total_size = std::accumulate(local_sizes.begin(), local_sizes.end(), 0) * data_elements;

    std::vector<double> all_data(total_size);

    MPI_Gatherv(local_data_send.data(), local_size * data_elements, MPI_DOUBLE,
        all_data.data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, comm);

    if (rank == 0) {
        std::ofstream outfile(file_name);

        if (outfile.is_open()) {

            outfile << "X,Y,Z,WSS";
            if (isPulsatile) {
                outfile << ",OSI,TAWSS";
            }
            outfile << ",VelMag\n";

            for (size_t i = 0; i < total_size; i += data_elements) {
                outfile << all_data[i] << "," << all_data[i + 1] << "," << all_data[i + 2] << ","
                    << all_data[i + 3];
                if (isPulsatile) {
                    outfile << "," << all_data[i + 4] << "," << all_data[i + 5];
                }
                outfile << "," << all_data[i + (data_elements - 1)] << "\n";
            }
            outfile.close();
        }
        else {
            std::cerr << "Unable to open file " << file_name << std::endl;
        }
    }
}






double OutputFile::calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
 
   return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2) + std::pow(z2 - z1, 2));
}


double OutputFile::calculateMeanNeighboringValues(const std::vector<DataPoint>& data, const DataPoint& point, double searchRadius, char attribute) {
    double sum = 0.0;
    int count = 0;

    for (const auto& neighbor : data) {
        double distance = calculateDistance(point.x, point.y, point.z, neighbor.x, neighbor.y, neighbor.z);
        if (distance > 0 && distance <= searchRadius) {
            switch (attribute) {
            case 'W':
                sum += neighbor.WSS;
                break;
            case 'O':
                sum += neighbor.OSI;
                break;
            case 'T':
                sum += neighbor.TAWSS;
                break;
            }
            count++;
        }
    }

    if (count > 0) {
        return sum / count;
    }
    else {
        switch (attribute) {
        case 'W':
            return point.WSS;
        case 'O':
            return point.OSI;
        case 'T':
            return point.TAWSS;
        }
    }

    return 0.0;
}


void OutputFile::smoothData(std::vector<DataPoint>& data, double neighborSearchRadius, const std::string& flowType) {
    std::vector<DataPoint> modifiedData = data;

#pragma omp parallel for
    for (size_t i = 0; i < data.size(); ++i) {
        const auto& point = data[i];

        double meanWSS = calculateMeanNeighboringValues(data, point, neighborSearchRadius, 'W');
        modifiedData[i].WSS = meanWSS;

        if (flowType == "Pulsatile") {
            double meanOSI = calculateMeanNeighboringValues(data, point, neighborSearchRadius, 'O');
            double meanTAWSS = calculateMeanNeighboringValues(data, point, neighborSearchRadius, 'T');
            modifiedData[i].OSI = meanOSI;
            modifiedData[i].TAWSS = meanTAWSS;
        }

        modifiedData[i].velMag = point.velMag;
    }

    data = modifiedData;
}




void OutputFile::processFlowData(const std::string& inputFileName, const std::string& outputFileName,
    double neighborSearchRadius, const std::string& flowType) {
    std::string folderPath = outputFileName.substr(0, outputFileName.find_last_of('/'));

#ifdef _WIN32
    std::filesystem::create_directory(folderPath);
#else
    mkdir(folderPath.c_str(), 0777);
#endif

    std::vector<DataPoint> data;
    std::ifstream inputFile(inputFileName);
    std::ofstream outputFile(outputFileName);
    std::string line;

    if (std::getline(inputFile, line)) {
        outputFile << line << std::endl;
    }
    else {
      //  std::cerr << "Error: Input file is empty or missing header: " << inputFileName << std::endl;
        return;
    }

    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string temp;
        DataPoint point;

        try {
            std::getline(ss, temp, ',');
            point.x = std::stod(temp);

            std::getline(ss, temp, ',');
            point.y = std::stod(temp);

            std::getline(ss, temp, ',');
            point.z = std::stod(temp);

            std::getline(ss, temp, ',');
            point.WSS = std::stod(temp);

            if (flowType == "Pulsatile") {
                std::getline(ss, temp, ',');
                point.OSI = std::stod(temp);

                std::getline(ss, temp, ',');
                point.TAWSS = std::stod(temp);
            }

            std::getline(ss, temp, ',');
            point.velMag = std::stod(temp);

            data.push_back(point);
        }
        catch (const std::exception& e) {
         //   std::cerr << "Error parsing line: " << line << " (" << e.what() << ")" << std::endl;
        }
}

    smoothData(data, neighborSearchRadius, flowType);

    for (const auto& point : data) {
        outputFile << point.x << "," << point.y << "," << point.z << ","
            << point.WSS;

        if (flowType == "Pulsatile") {
            outputFile << "," << point.OSI << "," << point.TAWSS;
        }

        outputFile << "," << point.velMag << std::endl;
    }

    inputFile.close();
    outputFile.close();





#ifdef _WIN32
    std::string command = "del " + inputFileName;
#else
    std::string command = "rm -f " + inputFileName;
#endif
//    int ret = system(command.c_str());
//    if (ret != 0) {
     //   std::cerr << "Error: Failed to remove temporary file: " << inputFileName << std::endl;
//    }
}

