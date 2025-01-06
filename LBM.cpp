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



#include <iostream>
#include "LatticeModel.h"
#include "Simulation.h"
#include "OutputFile.h"
#include <mpi.h>
#include "Subdomain.h"
#include "Input.h"
#include "ComplexGeometry.h"
#include <chrono>
#include <thread>
#include <iomanip>
#include <omp.h>
#include "Utility.h"
#include "GeometryReader.h"
#include "LBM.h"



int main(int argc, char* argv[]) {


    auto start = std::chrono::high_resolution_clock::now();
    int myrank, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);




    if (myrank == 0) {
        
        omp_set_num_threads(omp_get_max_threads());

        bool stopChrono = false;
        std::thread chronoThread([&]() {
            while (!stopChrono) {
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = now - start;
                std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
                    << std::setw(8) << elapsed.count() << " seconds" << std::flush;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            });
        chronoThread.detach();
    }

    SimulationParameters simParams;



    try {
        parseInputFile("config.txt", simParams);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input file: " << e.what() << std::endl;
     
    }
    if (myrank == 0) {

        printSimulationParameters(simParams);

    }

    std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();

    Geometry Geometry(velSet, simParams);
    SubDomain SubDomain(Geometry);
    if (myrank==0){std::cout<<"subdomain is done"<<std::endl;}
    Simulation Simulation(SubDomain, simParams);

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (myrank == 0) {
        std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
            << std::setw(8) << elapsed.count() << " seconds." << std::endl;

    }

    return 0;
}
