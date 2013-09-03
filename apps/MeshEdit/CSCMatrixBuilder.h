/*
 *  CSCMatrixBuilder.h
 *  GEL
 *	Added to simplify sequetnial building of sparse matrices in CSC format
 *  At this version this doesn't check if we are outside matrix (no matrix size)
 *  1. Insert the entries column by column (for symmetric matrices only up triangle) by using:
 *				insert_entry(row_number,entry);
 *				next_column()
 *				next_column_nonsort() - can be used instead of next_column where the columns were also filled
 *											sequentially (in their order)
 *		You can also use LoadDense or LoadSymmetricDense if you have dense representation of the matrix in
 *      LinAlg::CMatrix type.
 *  2. Call get_Matrix() to get the representation of the matrix as a ARluSymMatrix<Type> 
 *		(you can pass this type directly to the symmetric problem constructor).
 *  Created by Katarzyna Gebal on 18/11/08
 *	
 */
#ifndef __MESHEDIT_CSCMATRIXBUILDER_H__
#define __MESHEDIT_CSCMATRIXBUILDER_H__

#include <vector>
#include <LinAlg/Matrix.h>
#include "arlsmat.h"

template<class T>
class CSCMatrixBuilder
{
	typedef T Type;
	std::vector<Type> vA;
	std::vector<int> vpcol;
	std::vector<int> virow;

	ARluSymMatrix<Type> mat;
	int *irow;
	int *pcol;
	Type *A;
	bool mat_created;

	int get_nnz()
	{
		return vA.size();
	}

	int get_ncol()
	{
		return vpcol.size()-1;
	}

	int* get_pcol()
	{
		int *pcol = new int[vpcol.size()];
		std::vector<int>::iterator cit = vpcol.begin();
		for(int i = 0; i < vpcol.size(); i++)
		{
			pcol[i] = *cit;
			++cit;
		}
		return pcol;
	}

	int* get_irow()
	{
		int nnz = vA.size();
		std::vector<int>::iterator rit = virow.begin();
		int *irow = new int[nnz];
		for(int i = 0; i < nnz; i++)
		{
			irow[i] = *rit;
			++rit;
		}
		return irow;
	}

	
	Type* get_A()
	{
		int nnz = vA.size();
		Type *A = new Type[nnz];
		for(int i = 0; i < nnz; i++)
			A[i] = vA[i];
		return A;
	}

	void create_mat()
	{
		if(mat_created)
			return;
		mat_created = true;
		irow = get_irow();               
		pcol = get_pcol();     
		A = get_A();
		mat = ARluSymMatrix<Type>(get_ncol(), get_nnz(), A, irow, pcol);	

	}
public:
	ARluSymMatrix<Type>& get_Matrix()
	{
		create_mat();
		return mat;
	}
	CSCMatrixBuilder()
	{
		mat_created = false;
		vpcol.push_back(0);
	}
	~CSCMatrixBuilder()
	{
		if(mat_created)
		{
			delete [] irow;
			delete [] pcol;
			delete [] A;
		}
	}

	void sort_entries(int id)
	{
		int start = vpcol[id];
		int stop = vpcol[id+1];
		Type dpom;
		int ipom;
		for(int i = start; i < stop; i++)
			for(int j = i+1; j < stop; j++)
			{
				if(virow[i] > virow[j])
				{
					ipom = virow[i]; virow[i] = virow[j]; virow[j] = ipom;
					dpom = vA[i]; vA[i] = vA[j]; vA[j] = dpom;
				}
			}
	}

	void insert_entry(int row, Type entry)
	{
		vA.push_back(entry);
		virow.push_back(row);
	}

	void next_column_nonsort()
	{
		vpcol.push_back(vA.size());
	}

	void next_column()
	{
		vpcol.push_back(vA.size());
		sort_entries(vpcol.size()-2);
	}

	void LoadSymmetricDense(LinAlg::CMatrix Q, Type eps)
	{
		 for(int i = 0; i < Q.Cols(); i++)
		 {
			 for(int j = i; j < Q.Rows(); j++)
			 {
				 if(Q[i][j] > eps || Q[i][j] < -eps)
					 insert_entry(j,Q[i][j]);
			 }
			 next_column();
		 }
	}

	void LoadDense(LinAlg::CMatrix Q, Type eps)
	{
		 for(int i = 0; i < Q.Cols(); i++)
		 {
			 for(int j = 0; j < Q.Rows(); j++)
			 {
				 if(Q[i][j] > eps || Q[i][j] < -eps)
					 insert_entry(j,Q[i][j]);
			 }
			 next_column();
		 }
	}

};

#endif