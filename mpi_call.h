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
#ifndef _MPI_CALL_H
#define _MPI_CALL_H

#include "include.h"
#include "typedef.h"


template<class T>

class MPI_CALL {

public:
	int sendto_;
	int receive_;
	char name;
	MPI_Datatype mpi_type_;
	MPI_Status status_;
	size_t Size_;
	size_t vec_sizes;

	template <typename T>
	std::vector<T> flatten(const std::vector<std::vector<T>>& vec) {
		std::vector<T> result;
		for (const auto& v : vec)
			result.insert(result.end(), v.begin(), v.end());
		return result;
	}

	void Send_to(size_t sendto, std::vector<T>& Vector);
	void Recv_from(size_t receiv, std::vector<T>& Vector);

	void Send_to(size_t sendto, std::vector<std::vector<T>>& VecVec);
	void Recv_from(size_t receive, std::vector<std::vector<T>>& VecVec);

	void Send_to(size_t sendto, std::vector<std::vector<T>>& VecVec, size_t Vec_size);
	void Recv_from(size_t receive, std::vector<std::vector<T>>& VecVec, size_t vec_size);

	void Broadcast(int myrank, std::vector<T>& Vec);

	MPI_Datatype createBoundingBoxMPIType() {
		int blocklengths[6] = { 1, 1, 1, 1, 1, 1 };
		MPI_Datatype types[6] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
		MPI_Aint displacements[6];

		displacements[0] = offsetof(BoundingBox, min_x);
		displacements[1] = offsetof(BoundingBox, max_x);
		displacements[2] = offsetof(BoundingBox, min_y);
		displacements[3] = offsetof(BoundingBox, max_y);
		displacements[4] = offsetof(BoundingBox, min_z);
		displacements[5] = offsetof(BoundingBox, max_z);

		MPI_Datatype mpiBoundingBoxType;
		MPI_Type_create_struct(6, blocklengths, displacements, types, &mpiBoundingBoxType);
		MPI_Type_commit(&mpiBoundingBoxType);

		return mpiBoundingBoxType;
	}

	MPI_Datatype createTriangleMPIType() {
		MPI_Datatype mpiTriangleType;

		MPI_Datatype types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
		int blockcounts[4] = { 3, 3, 3, 3 };  
		MPI_Aint offsets[4];

		offsets[0] = offsetof(Triangle, x1);
		offsets[1] = offsetof(Triangle, x2);
		offsets[2] = offsetof(Triangle, x3);
		offsets[3] = offsetof(Triangle, nx);

		MPI_Type_create_struct(4, blockcounts, offsets, types, &mpiTriangleType);  
		MPI_Type_commit(&mpiTriangleType);

		return mpiTriangleType;
	}


	//void SumVectors_in_all_procs(std::vector<T>& Vec);
	//void SumVectors_in_myproc(int myproc, std::vector<T>& Vec);

};

template<class T>
void MPI_CALL<T>::Send_to(size_t sendto, std::vector<T>& Vector)
{

	Size_ = Vector.size();
	MPI_Send(&Size_, 1, MPI_UNSIGNED, sendto, 0, MPI_COMM_WORLD);
	MPI_Send(Vector.data(), Size_, mpi_type_, sendto, 0, MPI_COMM_WORLD);

}
template<class T>
void MPI_CALL<T>::Recv_from(size_t receive, std::vector<T>& Vector)
{

	MPI_Recv(&Size_, 1, MPI_UNSIGNED, receive, 0, MPI_COMM_WORLD, &status_);
	Vector.resize(Size_, 0);
	MPI_Recv(Vector.data(), Size_, mpi_type_, receive, 0, MPI_COMM_WORLD, &status_);
}

template<class T>
void MPI_CALL<T>::Send_to(size_t sendto, std::vector<std::vector<T>>& VecVec)
{
	std::vector<T> Aux;
	std::vector<int> vec_sizes;

	for (size_t i = 0; i < VecVec.size(); i++) {
		vec_sizes.push_back(VecVec[i].size());
		for (size_t j = 0; j < VecVec[i].size(); j++) {
			Aux.push_back(VecVec[i][j]);
		}
	}

	int Size_ = vec_sizes.size();
	MPI_Send(&Size_, 1, MPI_INT, sendto, 0, MPI_COMM_WORLD);
	MPI_Send(vec_sizes.data(), Size_, MPI_INT, sendto, 0, MPI_COMM_WORLD);
	int SizeAux_ = Aux.size();
	MPI_Send(&SizeAux_, 1, mpi_type_, sendto, 0, MPI_COMM_WORLD);
	MPI_Send(Aux.data(), SizeAux_, mpi_type_, sendto, 0, MPI_COMM_WORLD);
}

template<class T>
void MPI_CALL<T>::Recv_from(size_t receive, std::vector<std::vector<T>>& VecVec)
{
	int Size_ = 0;
	int SizeAux_ = 0;
	MPI_Recv(&Size_, 1, MPI_INT, receive, 0, MPI_COMM_WORLD, &status_);
	std::vector<int> vec_sizes(Size_);
	MPI_Recv(vec_sizes.data(), Size_, MPI_INT, receive, 0, MPI_COMM_WORLD, &status_);
	MPI_Recv(&SizeAux_, 1, mpi_type_, receive, 0, MPI_COMM_WORLD, &status_);
	std::vector<T> Aux(SizeAux_);
	MPI_Recv(Aux.data(), SizeAux_, mpi_type_, receive, 0, MPI_COMM_WORLD, &status_);

	std::vector<T>Vec;

	for (size_t i = 0; i < vec_sizes.size(); i++)
	{
		VecVec.push_back(std::vector<T>());
		for (size_t j = 0; j < vec_sizes[i]; j++)
		{
			VecVec[i].push_back(Aux[i * vec_sizes[i] + j]);
		}
	}
}

template<class T>
void MPI_CALL<T>::Send_to(size_t sendto, std::vector<std::vector<T>>& VecVec, size_t vec_size)
{
	Size_ = vec_size * VecVec.size();
	size_t row_size = VecVec.size();
	std::vector<T> Aux;

	for (size_t i = 0; i < VecVec.size(); i++) {
		for (size_t j = 0; j < vec_size; j++) {
			Aux.push_back(VecVec[i][j]);

		}
	}

	MPI_Send(&Size_, 1, MPI_UNSIGNED, sendto, 0, MPI_COMM_WORLD);
	MPI_Send(&row_size, 1, MPI_UNSIGNED, sendto, 0, MPI_COMM_WORLD);
	MPI_Send(&Aux[0], Size_, mpi_type_, sendto, 0, MPI_COMM_WORLD);
}


template<class T>
void MPI_CALL<T>::Recv_from(size_t receive, std::vector<std::vector<T>>& VecVec, size_t vec_size)
{
	size_t Size_;
	size_t row_size;
	MPI_Recv(&Size_, 1, MPI_UNSIGNED, receive, 0, MPI_COMM_WORLD, &status_);
	MPI_Recv(&row_size, 1, MPI_UNSIGNED, receive, 0, MPI_COMM_WORLD, &status_);
	std::vector<T> Aux(Size_);
	MPI_Recv(&Aux[0], Size_, mpi_type_, receive, 0, MPI_COMM_WORLD, &status_);

	for (int i = 0; i < row_size; i++)
	{
		VecVec.push_back(std::vector<T>());
		for (int j = 0; j < vec_size; j++)
		{
			VecVec[i].push_back(Aux[i * vec_size + j]);
		}

	}
}

template<class T>
void MPI_CALL<T>::Broadcast(int myrank, std::vector<T>& Vec) {

	MPI_Bcast(Vec.data(), Vec.size(), mpi_type_, myrank, MPI_COMM_WORLD);
}
#endif

