#ifndef LAPLACIAN_MATRIX_2D_INCLUDED
#define LAPLACIAN_MATRIX_2D_INCLUDED
#include "LaplacianMatrix/LaplacianMatrix1D.h"

//////////////////////
// FiniteElements2D //
//////////////////////
template<class Real,int Type1,int Degree1,int Type2=Type1,int Degree2=Degree1>
class FiniteElements2D
{
public:
	template<class Data>
	static bool UpSample	( Vector<Data>& in , Vector<Data>& out , int lowD1 , int lowD2 , int& highD1 , int& highD2 );
	template<class Data>
	static bool DownSample	( Vector<Data>& in , Vector<Data>& out , int highD1 , int highD2 , int& lowD1 , int& lowD2 );
	template<class Data>
	static bool UpSample	( const Data* in , Data* out , int lowD1 , int lowD2 , int& highD1 , int& highD2 );
	template<class Data>
	static bool DownSample	( const Data* in , Data* out , int highD1 , int highD2 , int& lowD1 , int& lowD2 );

	static bool UpSample	( const Real* in , Real* out , int lowD1 , int lowD2 , int& highD1 , int& highD2 , int Channels );
	static bool DownSample	( const Real* in , Real* out , int highD1 , int highD2 , int& lowD1 , int& lowD2 , int Channels );

	static bool IsUpSamplable(int lowD1,int lowD2,int& highD1,int& highD2);
	static bool IsDownSamplable(int highD1,int highD2,int& lowD1,int& lowD2);
	static bool UpSampleMatrix(SparseMatrix<Real>& M,int lowD1,int lowD2,int& highD1,int& highD2);
	static bool DownSampleMatrix(SparseMatrix<Real>& M,int highD1,int highD2,int& lowD1,int& lowD2);

	static bool TensorMatrix(int inDim1,int outDim1,const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,SparseMatrix<Real>& m);
	template<class Data>
	static bool TensorMatrixMultiply(int inDim1,int outDim1,int inDim2,int outDim2,const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,const Vector<Data>& in,Vector<Data>& out);
	template<class Data>
	static bool TensorMatrixMultiply(int inDim1,int outDim1,int inDim2,int outDim2,const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,const Data* in,Data* out);
	static bool TensorMatrixMultiply(int inDim1,int outDim1,int inDim2,int outDim2,const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,const Real* in,Real* out,int Channels);
	static bool DotProductMatrix(int dim1,int dim2,SparseMatrix<Real>& m);
	static bool DotProductMatrix(int dim1,int dim2,SparseMatrix<Real>& m,Real min1,Real max1,Real min2,Real max2);
	static bool LaplacianMatrix(int dim1,int dim2,SparseMatrix<Real>& m,bool weakForm,bool negate=false);
	static bool LaplacianMatrix(int dim1,int dim2,SparseMatrix<Real>& m,Real min1,Real max1,Real min2,Real max2,bool weakForm,bool negate=false);
	static bool LaplacianMatrix(int dim1,int dim2,double iWeight , double gWeight , SparseMatrix<Real>& m,bool weakForm,bool negate=false);
	static bool LaplacianMatrix(int dim1,int dim2,double iWeight , double gWeight , SparseMatrix<Real>& m,Real min1,Real max1,Real min2,Real max2,bool weakForm,bool negate=false);
	static void StripDiagonal(const SparseMatrix<Real>& M,Vector<MatrixEntry<Real> >& D,bool clearDiagonal);

	template<class Data>
	static bool Gradient(const Vector<Data>& in,Vector<Data>& dX,Vector<Data>& dY,int w,int h,int& ww,int& hh);
	template<class Data>
	static bool Laplacian(const Vector<Data>& dX,const Vector<Data>& dY,Vector<Data>& out,int w,int h,int ww,int hh);

	class FullDivergenceStencil
	{
	public:
		class DivergenceStencil
		{
		public:
			Real values1[2*Degree1][2*Degree2+1];
			Real values2[2*Degree1+1][2*Degree2];
		};
		DivergenceStencil caseTable[2*Degree1+1][2*Degree2+2];
	};
	static bool DivergenceStencil(int dim1,int dim2,FullDivergenceStencil& s);

	class FullMatrixStencil
	{
	public:
		class MatrixStencil
		{
		public:
			Real values[2*Degree1+1][2*Degree2+1];
		};
		MatrixStencil caseTable[2*Degree1+1][2*Degree2+1];
	};
	class FullPaddedMatrixStencil
	{
	public:
		class PaddedMatrixStencil
		{
		public:
			Real values[2*Degree1+1][2*Degree2+1];
		};
		PaddedMatrixStencil caseTable[2*Degree1+3][2*Degree2+3];
	};
	static bool BiLaplacianStencil( int dim1 , int dim2 , FullMatrixStencil& s , bool weakForm , bool negate=false );
	static bool LaplacianStencil  ( int dim1 , int dim2 , FullMatrixStencil& s , bool weakForm , bool negate=false );
	static bool DotProductStencil ( int dim1 , int dim2 , FullMatrixStencil& s );

	class FullProlongationStencil
	{
	public:
		class ProlongationStencil
		{
		public:
			Real values[Degree1+2][Degree2+2];
		};
		ProlongationStencil caseTable[2*Degree1+1][2*Degree2+1];
	};
	static bool ProlongationStencil(int lowD1,int lowD2,FullProlongationStencil &s,int& highD1,int& highD2);
#if NEW_LAPLACIAN_CODE
	class FullRestrictionStencil
	{
	public:
		class RestrictionStencil
		{
		public:
			Real values[(Degree1+3)>>1][(Degree2+3)>>1];
		};
		RestrictionStencil caseTable[2*Degree1+2][2*Degree2+2];
	};
	static bool RestrictionStencil(int highD1,int highD2,FullRestrictionStencil &s,int& lowD1,int& lowD2);
#endif // NEW_LAPLACIAN_CODE
};
#include "LaplacianMatrix/LaplacianMatrix2D.inl"
#endif // LAPLACIAN_MATRIX_INCLUDED