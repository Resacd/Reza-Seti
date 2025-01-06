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



#include "ComplexGeometry.h"
#include <queue>
#include <unordered_set>
#include <map>


Geometry::Geometry(std::shared_ptr<VelocitySet> velSet, SimulationParameters &simParams)
{

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    velSet_ = velSet;
    this->simParams = simParams;
    spacing=simParams.dx;
    stlFile=simParams.stlFile;
    GetBoundingbox();
    DistributeTriangles();
    generatePointsInBoundingBox();
    labelPoints();
//    Writer.writeLabeledPointsToFile(pointsBox, "LabeledPointsFromGeometry.vtk");
    boundaryPoints = findBoundaryPoints();
    getInOutBP(); 
    WallboundaryPoints = excludeInletOutlet();
    if (myrank == 0) {

        writeBoundaryPoints();
        std::cout << std::endl;
        std::cout << "Geometry Reconstruction is done!" << std::endl;

    }

}




Geometry::~Geometry() {


}

void Geometry::GetBoundingbox() {

    globalTriangles = Reader.readSTLASCII(stlFile);

    if (myrank == 0) {

        globalNumtris = globalTriangles.size();
        Point minPoint = { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
        Point maxPoint = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };


        for (int t = 0; t < globalTriangles.size(); ++t) {
            const auto& triangle = globalTriangles[t];
            {
                minPoint.x = std::min(minPoint.x, std::min(triangle.x1, std::min(triangle.x2, triangle.x3)));
                minPoint.y = std::min(minPoint.y, std::min(triangle.y1, std::min(triangle.y2, triangle.y3)));
                minPoint.z = std::min(minPoint.z, std::min(triangle.z1, std::min(triangle.z2, triangle.z3)));
                maxPoint.x = std::max(maxPoint.x, std::max(triangle.x1, std::max(triangle.x2, triangle.x3)));
                maxPoint.y = std::max(maxPoint.y, std::max(triangle.y1, std::max(triangle.y2, triangle.y3)));
                maxPoint.z = std::max(maxPoint.z, std::max(triangle.z1, std::max(triangle.z2, triangle.z3)));
            }
        }

        boundingbox = { minPoint.x,maxPoint.x, minPoint.y,maxPoint.y, minPoint.z, maxPoint.z };
    }
    MPI_Bcast(&globalNumtris, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Datatype mpiBoundingBoxType = mpi_caller.createBoundingBoxMPIType();
    MPI_Bcast(&boundingbox, 1, mpiBoundingBoxType, 0, MPI_COMM_WORLD);


    nx = std::ceil((boundingbox.max_x - boundingbox.min_x) / spacing);
    ny = std::ceil((boundingbox.max_y - boundingbox.min_y) / spacing);
    nz = std::ceil((boundingbox.max_z - boundingbox.min_z) / spacing);

    GeomSize = { nx + 1, ny + 1, nz + 1 };
}



void Geometry::generatePointsInBoundingBox() {

    nx = std::ceil((boundingbox.max_x - boundingbox.min_x) / spacing);
    ny = std::ceil((boundingbox.max_y - boundingbox.min_y) / spacing);
    nz = std::ceil((boundingbox.max_z - boundingbox.min_z) / spacing);

    GeomSize = { nx + 1, ny + 1, nz + 1 };

    pointsBox.resize(nx + 1);
    for (int i = 0; i <= nx; ++i) {
        pointsBox[i].resize(ny + 1);
        for (int j = 0; j <= ny; ++j) {
            pointsBox[i][j].resize(nz + 1);
        }
    }

    double x, y, z;

    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                x = boundingbox.min_x + i * spacing;
                y = boundingbox.min_y + j * spacing;
                z = boundingbox.min_z + k * spacing;
                Point p{ x, y, z };
                pointsBox[i][j][k] = p;
            }
        }
    }
}


void Geometry::DistributeTriangles() {
    int localNumTriangles = globalNumtris / num_procs;
    int remainingTriangles = globalNumtris % num_procs;
    spacing*=2;
    std::vector<int> sendcounts(num_procs, localNumTriangles);

    for (int i = 0; i < remainingTriangles; ++i) {
        sendcounts[i]++;
    }

    std::vector<int> displacements(num_procs, 0);

    for (int i = 1; i < num_procs; ++i) {
        displacements[i] = displacements[i - 1] + sendcounts[i - 1];
    }


    localTriangles.resize(sendcounts[myrank]);
    MPI_Datatype mpiTriangleType = mpi_caller.createTriangleMPIType();
    MPI_Scatterv(globalTriangles.data(), sendcounts.data(), displacements.data(), mpiTriangleType,
        localTriangles.data(), sendcounts[myrank], mpiTriangleType, 0, MPI_COMM_WORLD);

}



bool Geometry::doLinesIntersect(const Line& line1, const Line& line2) {

    double x1 = line1.start.x;
    double y1 = line1.start.y;
    double z1 = line1.start.z;
    double x2 = line1.end.x;
    double y2 = line1.end.y;
    double z2 = line1.end.z;

    double x3 = line2.start.x;
    double y3 = line2.start.y;
    double z3 = line2.start.z;
    double x4 = line2.end.x;
    double y4 = line2.end.y;
    double z4 = line2.end.z;

    double ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));
    double ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));

    const double epsilon = 1e-7;
    if (ua >= -epsilon && ua <= 1 + epsilon && ub >= -epsilon && ub <= 1 + epsilon) {
        return true;
    }

    return false;
}


bool Geometry::isPointInsideTriangle(double x, double y, double z, const Triangle& triangle) {

    double x1 = triangle.x1, y1 = triangle.y1, z1 = triangle.z1;
    double x2 = triangle.x2, y2 = triangle.y2, z2 = triangle.z2;
    double x3 = triangle.x3, y3 = triangle.y3, z3 = triangle.z3;


    double dirX = 0.0, dirY = 0.0, dirZ = 1.0;


    double uX = x2 - x1, uY = y2 - y1, uZ = z2 - z1;
    double vX = x3 - x1, vY = y3 - y1, vZ = z3 - z1;


    double nX = uY * vZ - uZ * vY;
    double nY = uZ * vX - uX * vZ;
    double nZ = uX * vY - uY * vX;


    double nLengthSq = nX * nX + nY * nY + nZ * nZ;
    const double tol = 1e-16;
    if (nLengthSq < tol * tol) return false;


    double w0X = x - x1, w0Y = y - y1, w0Z = z - z1;


    double a = -(nX * w0X + nY * w0Y + nZ * w0Z);
    double b = nX * dirX + nY * dirY + nZ * dirZ;
    if (fabs(b) < tol) return false;


    double pos = a / b;


    double ix = x + pos * dirX;
    double iy = y + pos * dirY;
    double iz = z + pos * dirZ;


    double uu = uX * uX + uY * uY + uZ * uZ;
    double uv = uX * vX + uY * vY + uZ * vZ;
    double vv = vX * vX + vY * vY + vZ * vZ;
    double wX = ix - x1, wY = iy - y1, wZ = iz - z1;
    double wu = wX * uX + wY * uY + wZ * uZ;
    double wv = wX * vX + wY * vY + wZ * vZ;
    double D = uv * uv - uu * vv;

    double s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0) return false;
    double t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0) return false;

    return true;
}


int Geometry::countIntersectinsTriangle(Line& line, std::vector<Triangle>& posintersectedtri) {
    const double threshold = 1e-16;
    int count = 0;

    for (const auto& triangle : posintersectedtri) {
        double px = line.start.x;
        double py = line.start.y;
        double pz = line.start.z;
        double dx = line.end.x - line.start.x;
        double dy = line.end.y - line.start.y;
        double dz = line.end.z - line.start.z;


        double nx = triangle.nx;
        double ny = triangle.ny;
        double nz = triangle.nz;
        double qx = triangle.x1;
        double qy = triangle.y1;
        double qz = triangle.z1;


        double numerator = (qx - px) * nx + (qy - py) * ny + (qz - pz) * nz;
        double denominator = dx * nx + dy * ny + dz * nz;


        if (fabs(denominator) < threshold) {
            continue;
        }

        double t = numerator / denominator;


        if (t >= 0.0 && t <= 1.0) {

            double ix = px + t * dx;
            double iy = py + t * dy;
            double iz = pz + t * dz;

            if (isPointInsideTriangle(ix, iy, iz, triangle)) {
                count++;
            }
        }
    }

    return count;
}



void Geometry::adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch) {
    switch (directionSwitch) {
    case 0:  // positive z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
        }

        if (line.end.z == boundingbox.max_z) {
            line.end.z -= rayOffset;
        }

        break;

    case 1:  // negative z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }


        if (line.start.z == boundingbox.min_z) {
            line.start.z -= rayOffset;
        }

        break;

    case 2:  // positive x-direction

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        break;


    case 3:  // negative x-direction


    case 4:  // positive y-direction
    case 5:  // negative y-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }
        break;
    }
}



void Geometry::labelPoints() {
    const double rayOffset = 0.01;
    std::vector<Triangle> possiblelocalIntersectedTriangles;
    int directionSwitch;
    int local_intersected_triangles;
    int total_intersected_triangles;
    int intersectingCount;
    rayLength = boundingbox.max_z - boundingbox.min_z;
    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                Line line;
                line.start.x = pointsBox[i][j][k].x;
                line.start.y = pointsBox[i][j][k].y;
                line.start.z = pointsBox[i][j][k].z;

                directionSwitch = 0;
                local_intersected_triangles = 0;
                total_intersected_triangles = 0;
                bool foundNonZeroIntersections = false;

  //              do {
                switch (directionSwitch) {
                case 0:  // positive z-direction
                    rayLength = boundingbox.max_z - boundingbox.min_z;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z + rayLength + 0.1;
                    break;
                case 1:  // negative z-direction
                    rayLength = boundingbox.max_z - boundingbox.min_z;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z - (rayLength + 0.1);
                    break;
                case 2:  // positive x-direction
                    rayLength = boundingbox.max_x - boundingbox.min_x;
                    line.end.x = pointsBox[i][j][k].x + rayLength + 0.1;
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                case 3:  // negative x-direction
                    rayLength = boundingbox.max_x - boundingbox.min_x;
                    line.end.x = pointsBox[i][j][k].x - (rayLength + 0.1);
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                case 4:  // positive y-direction
                    rayLength = boundingbox.max_y - boundingbox.min_y;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y + rayLength + 0.1;
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                case 5:  // negative y-direction
                    rayLength = boundingbox.max_y - boundingbox.min_y;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y - (rayLength + 0.1);
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                }

                adjustLineCoordinates(line, boundingbox, rayOffset, directionSwitch);
                local_intersected_triangles = countIntersectinsTriangle(line, localTriangles);
                MPI_Allreduce(&local_intersected_triangles, &total_intersected_triangles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

                if (total_intersected_triangles % 2 == 1) {
                    pointsBox[i][j][k].Type = NodeType::Fluid;
                }
                else {
                    pointsBox[i][j][k].Type = NodeType::NotFluid;
                }
            }
        }
    }
}


std::vector<Point> Geometry::findBoundaryPoints() {
    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                Point& point = pointsBox[i][j][k];

                if (point.Type == NodeType::Fluid) {
                    bool isBoundaryPoint = false;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                if (i + dx >= 0 && i + dx < pointsBox.size() &&
                                    j + dy >= 0 && j + dy < pointsBox[i].size() &&
                                    k + dz >= 0 && k + dz < pointsBox[i][j].size()) {
                                    Point& neighbor = pointsBox[i + dx][j + dy][k + dz];
                                    if (neighbor.Type == NodeType::NotFluid) {
                                        isBoundaryPoint = true;
                                        break;
                                    }
                                }
                            }
                            if (isBoundaryPoint) {
                                break;
                            }
                        }
                        if (isBoundaryPoint) {
                            break;
                        }
                    }
                    bool isAtBoundary = (i == 0 || i == pointsBox.size() - 1) ||
                        (j == 0 || j == pointsBox[i].size() - 1) ||
                        (k == 0 || k == pointsBox[i][j].size() - 1);

                    if (isAtBoundary && !isBoundaryPoint && point.Type == NodeType::Fluid) {
                        isBoundaryPoint = true;
                    }


                    if (isBoundaryPoint) {
                        boundaryPoints.push_back(point);
                    }
                }
            }
        }
    }

    return boundaryPoints;
}


void Geometry::getInOutBP() {
    const long double tolerance = 2.5e-2;

    InletboundaryPoints.first.clear();
    InletboundaryPoints.second.clear();
    OutletboundaryPoints.clear();
    OutletboundaryPoints.resize(simParams.outlets.size());


    for (const auto& inlet : simParams.inlets) {
        const std::string& plane = inlet.plane;
        double coord = inlet.coordinate;
        const std::vector<double>& normal = inlet.normal;
        const std::vector<double>& rangeX = inlet.range_x;
        const std::vector<double>& rangeY = inlet.range_y;
        const std::vector<double>& rangeZ = inlet.range_z;

        for (const auto& point : boundaryPoints) {
            bool isInlet = false;

            if (plane == "x" && std::abs(point.x - coord) <= tolerance) {
                isInlet = (point.y >= rangeY[0] && point.y <= rangeY[1] &&
                    point.z >= rangeZ[0] && point.z <= rangeZ[1]);
            }
            else if (plane == "y" && std::abs(point.y - coord) <= tolerance) {
                isInlet = (point.x >= rangeX[0] && point.x <= rangeX[1] &&
                    point.z >= rangeZ[0] && point.z <= rangeZ[1]);
            }
            else if (plane == "z" && std::abs(point.z - coord) <= tolerance) {
                isInlet = (point.x >= rangeX[0] && point.x <= rangeX[1] &&
                    point.y >= rangeY[0] && point.y <= rangeY[1]);
            }

            if (isInlet) {
                InletboundaryPoints.first.push_back(point);
            }
        }

        InletboundaryPoints.second = normal;
    }


    for (size_t i = 0; i < simParams.outlets.size(); ++i) {
        const auto& outlet = simParams.outlets[i];
        const std::string& plane = outlet.plane;
        double coord = outlet.coordinate;
        const std::vector<double>& normal = outlet.normal;
        const std::vector<double>& rangeX = outlet.range_x;
        const std::vector<double>& rangeY = outlet.range_y;
        const std::vector<double>& rangeZ = outlet.range_z;

        for (const auto& point : boundaryPoints) {
            bool isOutlet = false;

            if (plane == "x" && std::abs(point.x - coord) <= tolerance) {
                isOutlet = (point.y >= rangeY[0] && point.y <= rangeY[1] &&
                    point.z >= rangeZ[0] && point.z <= rangeZ[1]);
            }
            else if (plane == "y" && std::abs(point.y - coord) <= tolerance) {
                isOutlet = (point.x >= rangeX[0] && point.x <= rangeX[1] &&
                    point.z >= rangeZ[0] && point.z <= rangeZ[1]);
            }
            else if (plane == "z" && std::abs(point.z - coord) <= tolerance) {
                isOutlet = (point.x >= rangeX[0] && point.x <= rangeX[1] &&
                    point.y >= rangeY[0] && point.y <= rangeY[1]);
            }

            if (isOutlet) {
                OutletboundaryPoints[i].first.push_back(point);
            }
        }

        OutletboundaryPoints[i].second = normal;
    }
}



std::vector<Point> Geometry::excludeInletOutlet() {
    std::vector<Point> remainingPoints;

    std::unordered_set<Point> inletSet(InletboundaryPoints.first.begin(), InletboundaryPoints.first.end());
    std::unordered_set<Point> outletSet;

    for (const auto& outlet : OutletboundaryPoints) {
        outletSet.insert(outlet.first.begin(), outlet.first.end());
    }

    for (const Point& p : boundaryPoints) {
        if (inletSet.find(p) == inletSet.end() && outletSet.find(p) == outletSet.end()) {
            remainingPoints.push_back(p);
        }
    }

    return remainingPoints;
}



void Geometry::writeBoundaryPoints() {

    if (!InletboundaryPoints.first.empty()) {
        Writer.writePointsToFile(InletboundaryPoints.first, "InboundaryPointsFromGeometry.vtk");
    }
    else {
        std::cerr << "Warning: No inlet boundary points to write." << std::endl;
    }

    for (size_t i = 0; i < OutletboundaryPoints.size(); ++i) {
        const auto& outletPoints = OutletboundaryPoints[i].first;
        if (!outletPoints.empty()) {
            std::string fileName = "OutboundaryPointsFromGeometry" + std::to_string(i + 1) + ".vtk";
            Writer.writePointsToFile(outletPoints, fileName);
        }
        else {
            std::cerr << "Warning: Outlet " << i + 1 << " has no boundary points to write." << std::endl;
        }
    }

    if (!WallboundaryPoints.empty()) {
        Writer.writePointsToFile(WallboundaryPoints, "WallboundaryPointsFromGeometry.vtk");
    }
    else {
        std::cerr << "Warning: No wall boundary points to write." << std::endl;
    }
}
