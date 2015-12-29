#ifndef LAPLACIAN_MATRIX_1D_INCLUDED
#define LAPLACIAN_MATRIX_1D_INCLUDED
#include "LinearAlgebra/SparseMatrix.h"
#include "FunctionBasis/PPolynomial.h"

#define MISHA_CODE 0

#define NEW_LAPLACIAN_CODE 1
enum
{
	ZERO_VALUE,
	ZERO_DERIVATIVE,
	PERIODIC_BOUNDARY,
	TOTAL_BOUNDARY_TYPES
};
#define DERIVATIVE(Type) ( ((~(Type&1))&1) | (Type&2) )

static const char* BoundaryTypes[]=
{
	"Zero Boundary Values",
	"Zero Boundary Derivatives",
	"Periodic Boundary"
};

//////////////////////////
// BinomialCoefficients //
//////////////////////////
template<int Degree>
class BinomialCoefficients
{
public:
	int coeffs[Degree+1];
	BinomialCoefficients(void);
};

/////////////////////////
// FiniteDifferences1D //
/////////////////////////
template<class Real,int Type>
class FiniteDifferences1D
{
	static void SetValues(int dim,int i,MatrixEntry<Real>* values,bool makePositive);
public:
	static void GetMatrix		(const int& dim,SparseMatrix<Real>& m,bool makePositive);
	static void StripDiagonal	(const SparseMatrix<Real>& M,Vector<MatrixEntry<Real> >& D,bool clearDiagonal);
};


//////////////////////
// FiniteElements1D //
//////////////////////
template< class Real , int Type , int Degree >
class FiniteElements1D
{
	template<class Data>
	static bool UpSample(Vector<Data>& in,Vector<Data>& out,const Real* weights,int weightCount);
	template<class Data>
	static bool DownSample(Vector<Data>& in,Vector<Data>& out,const Real* weights,int weightCount);
	static bool UpSampleMatrix(SparseMatrix<Real>& M,int lowD,const Real* weights,int weightCount,int& highD);
	static bool DownSampleMatrix(SparseMatrix<Real>& M,int highD,const Real* weights,int weightCount,int& lowD);
public:
	static int DimensionOffset(void);
	static int DomainSize(int dim);
	static int Dimension(int domainSize);

	static int DerivativeMatrix(int dim,SparseMatrix<Real>& m);

	template<class Data>
	static bool UpSample(Vector<Data>& in,Vector<Data>& out);
	template<class Data>
	static bool DownSample(Vector<Data>& in,Vector<Data>& out);

	static bool IsDownSamplable(int highD,int& lowD);
	static bool IsUpSamplable(int lowD,int& highD);
	static bool UpSampleMatrix(SparseMatrix<Real>& M,int lowD,int& highD);
	static bool DownSampleMatrix(SparseMatrix<Real>& M,int highD,int& lowD);

	static bool IdentityMatrix(int dim,SparseMatrix<Real>& m,bool rowMajor=true);
	static bool DotProductMatrix(int dim,SparseMatrix<Real>& m);
	static bool DotProductMatrix(int dim,SparseMatrix<Real>& m,Real min,Real max);
	static bool LaplacianMatrix(int dim,SparseMatrix<Real>& m,bool weakForm,bool negate=false);
	static bool LaplacianMatrix(int dim,SparseMatrix<Real>& m,Real min,Real max,bool weakForm,bool negate=false);
	static void StripDiagonal(const SparseMatrix<Real>& M,Vector<MatrixEntry<Real> >& D,bool clearDiagonal);


	template<int DotType=Type,int DotDegree=Degree>
	class DotProduct
	{
	public:
		class Helper
		{
			double shift,dotShift;
			static const int Start,Stop;
			PPolynomial< Degree    , double > F1;
			PPolynomial< DotDegree , double > F2;
			PPolynomial< Degree+DotDegree , double > F1F2[Degree+DotDegree+1];
			double min[Degree+DotDegree+1],max[Degree+DotDegree+1];
			double fullValues[Degree+DotDegree+1];
		public:
			Helper(void);
			void setDerivatives(int d1,int d2);
			double GetValue(int i,int j,double min,double max) const;
			void SetValues(int dim,int i,MatrixEntry<Real>* values,Real min,Real max);

			static int StartOffset(void);
			static int StopOffset(void);
		};
		static bool DerivativeMatrix(int dim,SparseMatrix<Real>& m,int d1,int d2,bool negate=false);
		static bool DerivativeMatrix(int dim,SparseMatrix<Real>& m,Real min,Real max,int d1,int d2,bool negate=false);
		class FullDotProductStencil
		{
		public:
			class DotProductStencil
			{
			public:
				Real values[Degree+DotDegree+1];
			};
			DotProductStencil caseTable[2*DotDegree+1];
		};
		class FullPaddedDotProductStencil
		{
		public:
			class PaddedDotProductStencil
			{
			public:
				Real values[Degree+DotDegree+1];
			};
			PaddedDotProductStencil caseTable[2*DotDegree+3];
		};
		static bool DotProductStencil( int dim , FullDotProductStencil& s , int d1 , int d2 , bool negate=false );
		static bool PaddedDotProductStencil( int dim , FullPaddedDotProductStencil& s , int d1 , int d2 , bool negate=false );
	};
	class FullProlongationStencil
	{
	public:
		class ProlongationStencil
		{
		public:
			static const int Size=Degree+2;
			static int Start(int i) { return (i<<1)-((Degree+1)>>1); }
			Real values[Size];
		};
		ProlongationStencil caseTable[2*Degree+1];
		void Normalize(int lowD);
	};
	static bool ProlongationStencil(int lowD,FullProlongationStencil &s,int& highD);

#if NEW_LAPLACIAN_CODE
	// BADNESS!!! This is only going to work for Type==ZERO_DERIVATIVE indexing
	class FullRestrictionStencil
	{
	public:
		class RestrictionStencil
		{
		public:
			static const int Size=((Degree+1)>>1)<<1;
//			static const int Size=(Degree+3)>>1;
			static int Start(int i) { return (i-Degree+1)>>1; }
			Real values[Size];
		};
		RestrictionStencil caseTable[2*Degree+2];	// Need +2 to handle even/odd parity in the interior
	};
	static bool RestrictionStencil(int highD,FullRestrictionStencil &s,int& lowD);
#endif // NEW_LAPLACIAN_CODE
};

template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil);
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree>::FullProlongationStencil& pStencil,int lowDim,int highDim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil);
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree-1>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil);
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree+1>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil);

template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil);
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree>::FullProlongationStencil& pStencil,int lowDim,int highDim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil);
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree-1>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil);
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree+1>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil);
#include "LaplacianMatrix/LaplacianMatrix1D.inl"
#endif // LAPLACIAN_MATRIX_INCLUDED