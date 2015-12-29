/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef __SPARSEMATRIX_HPP
#define __SPARSEMATRIX_HPP

#include "Vector.h"
#include <vector>

template <class T>
struct MatrixEntry
{
	MatrixEntry( void )		{ N =-1; Value = 0; }
	MatrixEntry( int i )	{ N = i; Value = 0; }
	MatrixEntry( int i , T v ) { N = i ; Value = v; }
	int N;
	T Value;
};
template <class T>
struct MatrixEntry2
{
	MatrixEntry2( void )	{ inN = outN = -1; value = 0; }
	int inN,outN;
	T value;
};

template<class T> class SparseMatrix
{
	MatrixEntry<T>** m_ppElements;
public:
	bool rowMajor;
	int groups;
	int* groupSizes;

	MatrixEntry< T >* operator[] ( int idx );
	const MatrixEntry< T >* operator[] ( int idx ) const;

	SparseMatrix();
	SparseMatrix(const SparseMatrix& M);
	~SparseMatrix();
	SparseMatrix<T>& operator = (const SparseMatrix<T>& M);

	void Transpose(void);
	SparseMatrix( int groups , bool rowMajor=true );
	void Resize	( int groups , bool rowMajor=true);
	void SetGroupSize( int group , int count );
	int Entries(void);

	template<class T2>
	Vector<T2> operator * (const Vector<T2>& V) const;

	template<class T2>
	void Multiply			( const Vector<T2>& In, Vector<T2>& Out ) const;
	template<class T2>
	void MultiplyTranspose	( const Vector<T2>& In, Vector<T2>& Out ) const;
	template<class T2>
	bool Multiply			( const Vector<T2>& In, Vector<T2>& Out, int startRow, int stopRow ) const;

	template<class T2>
	static bool Multiply	( const SparseMatrix<T>& M, const Vector<MatrixEntry<T> >& D, const Vector<T2>& In, Vector<T2>& Out, int startRow, int stopRow );
	template<class T2>
	static bool Multiply	( const SparseMatrix<T>& M, const Vector<MatrixEntry2<T> >& D, const Vector<T2>& In, Vector<T2>& Out, int startRow, int stopRow );

	template<class T2>
	static int SolveJacobi				(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& Diagonal, const Vector<T2>& b, int Iterations, Vector<T2>& Solution, bool ZeroSolutionVector = true);

	// Note that the solver has been modified so that the default assumption is that the diagonal has _not_ been zeroed out!
	template<class T2> static int SolveGaussSeidel( const SparseMatrix<T>& M , const Vector<MatrixEntry<T> >& Diagonal , const Vector<T2>& b , int Iterations , Vector<T2>& Solution , bool reverseIndex = false , bool zeroDiagonal = false );
	template<class T2>
	static int SolveGaussSeidel( const SparseMatrix<T>& M , const Vector<MatrixEntry<T> >& Diagonal , const Vector<T2>& b , int Iterations , Vector<T2>& Solution , int start , int stop , bool reverseIndex = false , bool zeroDiagonal = false );
	template<class T2>
	static int SolveGaussSeidel(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& Diagonal, const Vector<T2>& b, int Iterations, const Vector<T2>& sorWeights , Vector<T2>& Solution, bool reverseIndex = false );
	template<class T2>
	static int SolveGaussSeidel(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& Diagonal, const Vector<T2>& b, int Iterations, const Vector<T2>& sorWeights , Vector<T2>& Solution, int start, int stop, bool reverseIndex = false );
	template<class T2>
	static int SolveGaussSeidel(const SparseMatrix<T>& M,const Vector<MatrixEntry2<T> >& Diagonal, const Vector<T2>& b, int Iterations, Vector<T2>& Solution, bool reverseIndex = false );
	template<class T2>
	static int SolveGaussSeidel(const SparseMatrix<T>& M,const Vector<MatrixEntry2<T> >& Diagonal, const Vector<T2>& b, int Iterations, Vector<T2>& Solution, int start, int stop, bool reverseIndex = false );


	template<class T2>
	static int SolveGaussSeidel2(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& Diagonal, const Vector<T2>& b, int Iterations, Vector<T2>& Solution, bool reverseIndex = false );
	template<class T2>
	static int SolveGaussSeidel2(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& Diagonal, const Vector<T2>& b, int Iterations, Vector<T2>& Solution, int start, int stop, bool reverseIndex = false );

	bool copyColumnMajor(SparseMatrix<T>& out,int inDim,int outDim);
	bool copyRowMajor	(SparseMatrix<T>& out,int inDim,int outDim);
};
// In order to be supported:
// 1] The Matrix class has to support the Multiply(Vector<Data>,Vector<Data>) method.
// 2] The data class has to support Data*Data returning a value castable to a double

template <class Matrix>
class PivotSymmetricMatrix
{
public:
	std::vector<std::pair<Matrix,std::pair<int,int> > > rightMatrices;
	std::pair<Matrix,std::pair<int,int> > pivot;

	template <class Data>
	Vector<Data> operator * (const Vector<Data>& in) const;

	template <class Data>
	void Multiply(const Vector<Data>& in,Vector<Data>& out) const;

	void setPivot(const Matrix& M,int inDim,int outDim);
	void push(const Matrix& M,int inDim,int outDim);
	void pop(void);
};
template <class Matrix,class Data>
static int SolveConjugateGradient(const Matrix& SPD,const Vector<Data>& b,const int& iters,Vector<Data>& solution,const double eps=1e-8 , bool printError=false );
#include "SparseMatrix.inl"

#endif // __SPARSEMATRIX_HPP