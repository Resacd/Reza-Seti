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

#ifndef _USER_DEFINED_STRUCTS_H_
#define _USER_DEFINED_STRUCTS_H_
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <string>


enum class BoundaryType {
    Inlet,
    Outlet,
    Wall,
    None
};


enum class NodeType {

    Fluid,
    Boundary,
    NotFluid,
    Halo
};



inline std::string nodeTypeToString(NodeType type) {
    switch (type) {
    case NodeType::Fluid: return "Fluid";
    case NodeType::Boundary: return "InBoundary";
    case NodeType::NotFluid: return "NotFluid";
    case NodeType::Halo: return "Halo";
    default: return "Unknown";
    }
}


struct Point {

    double x;
    double y;
    double z;


    NodeType Type;
    Point operator-(const Point& other) const {
        return { x - other.x, y - other.y, z - other.z };
    }

    Point operator+(const Point& other) const {
        return { x + other.x, y + other.y, z + other.z };
    }

    Point operator*(const double scalar) const {
        return { x * scalar, y * scalar, z * scalar };
    }

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};




struct Equation {
 
       std::string name; 
       std::vector<double> coeffs; 
};



struct Line {
    Point start;
    Point end;

    Point direction() const {
        return end - start;
    }
};

struct Triangle {
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    double nx, ny, nz;
};


struct BoundingBox {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
};

struct Edge {
    double x1, y1, z1;
    double x2, y2, z2;

    Edge(double x1, double y1, double z1, double x2, double y2, double z2)
        : x1(x1), y1(y1), z1(z1), x2(x2), y2(y2), z2(z2) {}

    bool operator==(const Edge& other) const {
        return (x1 == other.x1 && y1 == other.y1 && z1 == other.z1 && x2 == other.x2 && y2 == other.y2 && z2 == other.z2) ||
            (x1 == other.x2 && y1 == other.y2 && z1 == other.z2 && x2 == other.x1 && y2 == other.y1 && z2 == other.z1);
    }
};


namespace std {
    template <>
    struct hash<Edge> {
        std::size_t operator()(const Edge& e) const {
            std::size_t h1 = std::hash<double>()(e.x1);
            std::size_t h2 = std::hash<double>()(e.y1);
            std::size_t h3 = std::hash<double>()(e.z1);
            std::size_t h4 = std::hash<double>()(e.x2);
            std::size_t h5 = std::hash<double>()(e.y2);
            std::size_t h6 = std::hash<double>()(e.z2);
            return h1 ^ h2 ^ h3 ^ h4 ^ h5 ^ h6;
        }
    };
}


namespace std {
    template <>
    struct hash<Point> {
        std::size_t operator()(const Point& p) const {
            
            std::size_t hx = std::hash<double>{}(p.x);
            std::size_t hy = std::hash<double>{}(p.y);
            std::size_t hz = std::hash<double>{}(p.z);
            return hx ^ (hy << 1) ^ (hz << 2); 
        }
    };
}


struct NodeData {
    int global_x;
    int global_y;
    int global_z;
    double velX;
    double velY;
    double velZ;
    double velMag;
};

struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};


struct Neighbors {
    int rank;
    int coords[3];
};


struct DataPoint {
 
         double x, y, z;
         double OSI, TAWSS;
         double WSS,velMag;
};


#endif 
