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




#ifndef LBM_H
#define LBM_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <iomanip>
#include <limits>

#include "SimulationParameters.h"

std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t");
    size_t last = str.find_last_not_of(" \t");
    return (first == std::string::npos || last == std::string::npos) ? "" : str.substr(first, last - first + 1);
}

std::vector<double> parseRange(const std::string& value) {
    std::vector<double> result;
    size_t start = value.find('[');
    size_t end = value.find(']');
    if (start != std::string::npos && end != std::string::npos) {
        std::istringstream rangeStream(value.substr(start + 1, end - start - 1));
        std::string number;
        while (std::getline(rangeStream, number, ',')) {
            result.push_back(std::stod(trim(number)));
        }
    }
    return result;
}

void parseInputFile(const std::string& fileName, SimulationParameters& simParams) {
    std::ifstream file(fileName);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open input file: " + fileName);
    }

    std::string line;
    std::string currentSection;

    while (std::getline(file, line)) {
        line = trim(line);

        if (line.empty() || line[0] == '#') continue;

        if (line.back() == ':') {
            currentSection = line.substr(0, line.size() - 1);
            continue;
        }

        if (currentSection == "parameters") {
            std::istringstream iss(line);
            std::string key, value;
            if (std::getline(iss, key, ':') && std::getline(iss, value)) {
                key = trim(key);
                value = trim(value);
                if (key == "time_resolution") simParams.dt = std::stod(value);
                else if (key == "density") simParams.rho = std::stod(value);
                else if (key == "kinematic_viscosity") simParams.nu_blood = std::stod(value);
                else if (key == "relaxation_factor") simParams.tau = std::stod(value);
                else if (key == "grid_size") simParams.dx = std::stod(value);
                else if (key == "pressure") simParams.pressure = std::stod(value);
                else if (key == "CteVel") simParams.CteVel = std::stod(value);
                else if (key == "reynolds_number") simParams.reynoldsNumber = std::stod(value);
                else if (key == "number_of_outlets") simParams.num_Out = std::stoi(value);
                else if (key == "IterToWrite") simParams.IterToWrite = std::stoi(value);
                else if (key == "MaxIter") simParams.MaxIter = std::stoi(value);
                else if (key == "Systole") simParams.Systole = std::stod(value);
                else if (key == "Diastole") simParams.Diastole = std::stod(value);
                else if (key == "NumCycle") simParams.NumCycle = std::stod(value);
                else if (key == "FlowType") {
                    value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
                    simParams.FlowType = value;
                }
                else if (key == "stl_file") {
                    value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
                    simParams.stlFile = value;
                }
            }
        } else if (currentSection == "inlets" || currentSection == "outlets") {
            std::string plane;
            double coordinate;
            std::vector<double> rangeX, rangeY, rangeZ, normal;

            do {
                line = trim(line);
                if (line.empty()) continue;

                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    std::string key = trim(line.substr(0, pos));
                    std::string value = trim(line.substr(pos + 1));

                    if (key == "plane") {
                        plane = value;
                    } else if (key == "coordinate") {
                        coordinate = std::stod(value);
                    } else if (key == "range_x") {
                        rangeX = parseRange(value);
                    } else if (key == "range_y") {
                        rangeY = parseRange(value);
                    } else if (key == "range_z") {
                        rangeZ = parseRange(value);
                    } else if (key == "normal") {
                        normal = parseRange(value);
                    }
                }
            } while (std::getline(file, line) && !line.empty());

            if (rangeX.empty()) rangeX = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
            if (rangeY.empty()) rangeY = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
            if (rangeZ.empty()) rangeZ = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};

            if (currentSection == "inlets") {
                SimulationParameters::Inlet inlet = {plane, coordinate, rangeX, rangeY, rangeZ, normal};
                simParams.inlets.push_back(inlet);
            } else if (currentSection == "outlets") {
                SimulationParameters::Outlet outlet = {plane, coordinate, rangeX, rangeY, rangeZ, normal};
                simParams.outlets.push_back(outlet);
            }
        } else if (currentSection == "equations") {
            std::string equationName;
            std::vector<double> coeffs;

            size_t pos = line.find(':');
            if (pos != std::string::npos) {
                equationName = trim(line.substr(0, pos));
                std::string value = trim(line.substr(pos + 1));


                coeffs = parseRange(value);
                Equation equation = { equationName, coeffs };
                simParams.equations.push_back(equation);
            }
        }
    }

    file.close();
}

void printSimulationParameters(const SimulationParameters& simParams) {
    std::cout << std::fixed << std::setprecision(10);

    std::cout << "Simulation Parameters:\n";
    std::cout << "STL File: " << simParams.stlFile << "\n";
    std::cout << "Time Resolution (dt): " << simParams.dt << "\n";
    std::cout << "Density (rho): " << simParams.rho << "\n";
    std::cout << "Kinematic Viscosity (nu_blood): " << simParams.nu_blood << "\n";
    std::cout << "Relaxation Factor (tau): " << simParams.tau << "\n";
    std::cout << "Grid Size (dx): " << simParams.dx << "\n";
    std::cout << "Pressure: " << simParams.pressure << "\n";
    std::cout << "CteVel (u): " << simParams.CteVel << "\n";
    std::cout << "Reynolds Number: " << simParams.reynoldsNumber << "\n";
    std::cout << "Number of Outlets: " << simParams.num_Out << "\n";
    std::cout << "Systole: " << simParams.Systole << "\n";
    std::cout << "Diastole: " << simParams.Diastole << "\n";
    std::cout << "NumCycle: " << simParams.NumCycle << "\n";
    std::cout << "IterToWrite: " << simParams.IterToWrite << "\n";
    std::cout << "FlowType: " << simParams.FlowType << "\n";
    std::cout << "MaxIter: " << simParams.MaxIter << "\n";

    std::cout << "\nInlets:\n";
    for (const auto& inlet : simParams.inlets) {
        std::cout << "  Plane: " << inlet.plane << "\n";
        std::cout << "  Coordinate: " << inlet.coordinate << "\n";
        std::cout << "  Range X: ";
        for (const auto& x : inlet.range_x) std::cout << x << " ";
        std::cout << "\n  Range Y: ";
        for (const auto& y : inlet.range_y) std::cout << y << " ";
        std::cout << "\n  Range Z: ";
        for (const auto& z : inlet.range_z) std::cout << z << " ";
        std::cout << "\n  Normal: ";
        for (const auto& n : inlet.normal) std::cout << n << " ";
        std::cout << "\n";
    }

    std::cout << "\nOutlets:\n";
    for (const auto& outlet : simParams.outlets) {
        std::cout << "  Plane: " << outlet.plane << "\n";
        std::cout << "  Coordinate: " << outlet.coordinate << "\n";
        std::cout << "  Range X: ";
        for (const auto& x : outlet.range_x) std::cout << x << " ";
        std::cout << "\n  Range Y: ";
        for (const auto& y : outlet.range_y) std::cout << y << " ";
        std::cout << "\n  Range Z: ";
        for (const auto& z : outlet.range_z) std::cout << z << " ";
        std::cout << "\n  Normal: ";
        for (const auto& n : outlet.normal) std::cout << n << " ";
        std::cout << "\n";
    }

    std::cout << "\nEquations:\n";
    for (const auto& equation : simParams.equations) {
        std::cout << "  " << equation.name << ": ";
        for (const auto& coeff : equation.coeffs) {
            std::cout << coeff << " ";
        }
        std::cout << "\n";
    }
}

#endif




