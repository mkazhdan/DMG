#include "FunctionBasis/PPolynomial.h"

//////////////////////
// FiniteElements2D //
//////////////////////
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::TensorMatrix(int inDim1,int outDim1,const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,SparseMatrix<Real>& m)
{
	int stride;
	if(m1.rowMajor!=m2.rowMajor)	return  false;
	if(m1.rowMajor)	stride=inDim1;
	else			stride=outDim1;
	m.Resize(m1.groups*m2.groups,m1.rowMajor);
	for(int i=0;i<m1.groups;i++)
		for(int j=0;j<m2.groups;j++)
		{
			int idx=i+j*m1.groups;
			m.SetGroupSize(idx,m1.groupSizes[i]*m2.groupSizes[j]);
			for(int ii=0;ii<m1.groupSizes[i];ii++)
				for(int jj=0;jj<m2.groupSizes[j];jj++)
				{
					int iidx=ii+jj*m1.groupSizes[i];
					m.m_ppElements[idx][iidx].N=m1.m_ppElements[i][ii].N+m2.m_ppElements[j][jj].N*stride;
					m.m_ppElements[idx][iidx].Value=m1.m_ppElements[i][ii].Value*m2.m_ppElements[j][jj].Value;
				}
		}

	return true;
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::TensorMatrixMultiply(int inDim1,int outDim1,int inDim2,int outDim2,
																			  const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,
																			  const Vector<Data>& in,Vector<Data>& out)
{
	int stride;
	if(m1.rowMajor!=m2.rowMajor || in.Dimensions()!=inDim1*inDim2)	return  false;
	if(m1.rowMajor)	stride= inDim1;
	else			stride=outDim1;
	out.Resize(outDim1*outDim2);
	for(int i=0;i<m1.groups;i++)
		for(int j=0;j<m2.groups;j++)
		{
			int idx=i+j*m1.groups;
			for(int ii=0;ii<m1.groupSizes[i];ii++)
			{
				if(m1.m_ppElements[i][ii].N<0)	continue;
				for(int jj=0;jj<m2.groupSizes[j];jj++)
				{
					if(m2.m_ppElements[j][jj].N<0)	continue;
					int N=m1.m_ppElements[i][ii].N+m2.m_ppElements[j][jj].N*stride;
					Real Value=m1.m_ppElements[i][ii].Value*m2.m_ppElements[j][jj].Value;
					if(m1.rowMajor)	out[idx]+=in[N]*Value;
					else			out[N]+=in[idx]*Value;
				}
			}
		}
	return true;
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::TensorMatrixMultiply(int inDim1,int outDim1,int inDim2,int outDim2,
																			  const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,
																			  const Data* in,Data* out)
{
	int stride;
	if(m1.rowMajor!=m2.rowMajor)	return  false;
	if(m1.rowMajor)	stride= inDim1;
	else			stride=outDim1;
	for(int i=0;i<m1.groups;i++)
		for(int j=0;j<m2.groups;j++)
		{
			int idx=i+j*m1.groups;
			for(int ii=0;ii<m1.groupSizes[i];ii++)
			{
				if(m1.m_ppElements[i][ii].N<0)	continue;
				for(int jj=0;jj<m2.groupSizes[j];jj++)
				{
					if(m2.m_ppElements[j][jj].N<0)	continue;
					int N=m1.m_ppElements[i][ii].N+m2.m_ppElements[j][jj].N*stride;
					Real Value=m1.m_ppElements[i][ii].Value*m2.m_ppElements[j][jj].Value;
					if(m1.rowMajor)	out[idx]+=in[N]*Value;
					else			out[N]+=in[idx]*Value;
				}
			}
		}
	return true;
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::TensorMatrixMultiply(int inDim1,int outDim1,int inDim2,int outDim2,
																			  const SparseMatrix<Real>& m1,const SparseMatrix<Real>& m2,
																			  const Real* in,Real* out,int Channels)
{
	int stride;
	if(m1.rowMajor!=m2.rowMajor)	return  false;
	if(m1.rowMajor)	stride= inDim1;
	else			stride=outDim1;
	for(int c=0;c<Channels;c++)
		for(int i=0;i<m1.groups;i++)
			for(int j=0;j<m2.groups;j++)
			{
				int idx=i+j*m1.groups;
				for(int ii=0;ii<m1.groupSizes[i];ii++)
				{
					if(m1.m_ppElements[i][ii].N<0)	continue;
					for(int jj=0;jj<m2.groupSizes[j];jj++)
					{
						if(m2.m_ppElements[j][jj].N<0)	continue;
						Real Value=m1.m_ppElements[i][ii].Value*m2.m_ppElements[j][jj].Value;
#if 1
						int myIdx=i+j*m1.groups*Channels+c*m1.groups;
						int N=m1.m_ppElements[i][ii].N+m2.m_ppElements[j][jj].N*stride*Channels+c*stride;
						if(m1.rowMajor)	out[myIdx]+=in[N]*Value;
						else			out[N]+=in[myIdx]*Value;
#else
						int N=m1.m_ppElements[i][ii].N+m2.m_ppElements[j][jj].N*stride*Channels+c*stride;
						if(m1.rowMajor)	out[idx]+=in[N]*Value;
						else			out[N]+=in[idx]*Value;
#endif
					}
				}
			}
	return true;
}


template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::IsUpSamplable(int lowD1,int lowD2,int& highD1,int& highD2)
{
	return FiniteElements1D<Real,Type1,Degree1>::IsUpSamplable(lowD1,highD1) && FiniteElements1D<Real,Type2,Degree2>::IsUpSamplable(lowD2,highD2);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::IsDownSamplable(int highD1,int highD2,int& lowD1,int& lowD2)
{
	return FiniteElements1D<Real,Type1,Degree1>::IsDownSamplable(highD1,lowD1) && FiniteElements1D<Real,Type2,Degree2>::IsDownSamplable(highD2,lowD2);
}

template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DotProductMatrix(int dim1,int dim2,SparseMatrix<Real>& m)
{
	return DotProductMatrix(dim1,dim2,m,0,FiniteElements1D<Real,Type1,Degree1>::DomainSize(dim1),0,FiniteElements1D<Real,Type2,Degree2>::DomainSize(dim2));
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DotProductMatrix(int dim1,int dim2,SparseMatrix<Real>& m,Real min1,Real max1,Real min2,Real max2)
{
	SparseMatrix<Real> dot1,d2Dot1,dot2,d2Dot2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::DotProductMatrix(dim1,dot1,min1,max1) ||
		!FiniteElements1D<Real,Type2,Degree2>::DotProductMatrix(dim2,dot2,min2,max2) )
		return false;
	return TensorMatrix(dim1,dim1,dot1,dot2,m);
}


template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::LaplacianMatrix(int dim1,int dim2,SparseMatrix<Real>& m,bool weakForm,bool negate)
{
	return LaplacianMatrix(dim1,dim2,m,0,Real(FiniteElements1D<Real,Type1,Degree1>::DomainSize(dim1)),0,Real(FiniteElements1D<Real,Type2,Degree2>::DomainSize(dim2)),weakForm,negate);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::LaplacianMatrix(int dim1,int dim2,SparseMatrix<Real>& m,Real min1,Real max1,Real min2,Real max2,bool weakForm,bool negate)
{
	m.Resize(dim1*dim2);
	for(int i=0;i<dim1*dim2;i++) m.SetGroupSize(i,(2*Degree1+1)*(2*Degree2+1));
	SparseMatrix<Real> dot1,d2Dot1,dot2,d2Dot2;
	FiniteElements1D<Real,Type1,Degree1>::DotProductMatrix(dim1,dot1,min1,max1);
	FiniteElements1D<Real,Type2,Degree2>::DotProductMatrix(dim2,dot2,min2,max2);
	FiniteElements1D<Real,Type1,Degree1>::LaplacianMatrix(dim1,d2Dot1,min1,max1,weakForm,negate);
	FiniteElements1D<Real,Type2,Degree2>::LaplacianMatrix(dim2,d2Dot2,min2,max2,weakForm,negate);

	for(int y=0;y<dim2;y++)
		for(int x=0;x<dim1;x++)
			for(int yy=0;yy<(2*Degree2+1);yy++)
				for(int xx=0;xx<(2*Degree1+1);xx++)
				{
					m.m_ppElements[x+y*dim1][xx+yy*(2*Degree1+1)].N=
						dot1.m_ppElements[x][xx].N+dot2.m_ppElements[y][yy].N*dim1;
					m.m_ppElements[x+y*dim1][xx+yy*(2*Degree1+1)].Value=
						  dot1.m_ppElements[x][xx].Value*d2Dot2.m_ppElements[y][yy].Value+
						d2Dot1.m_ppElements[x][xx].Value*  dot2.m_ppElements[y][yy].Value;
				}
	return true;
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::LaplacianMatrix(int dim1,int dim2,double iWeight , double gWeight , SparseMatrix<Real>& m,bool weakForm,bool negate)
{
	return LaplacianMatrix(dim1,dim2,iWeight , gWeight , m,0,Real(FiniteElements1D<Real,Type1,Degree1>::DomainSize(dim1)),0,Real(FiniteElements1D<Real,Type2,Degree2>::DomainSize(dim2)),weakForm,negate);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::LaplacianMatrix(int dim1,int dim2 , double iWeight , double gWeight , SparseMatrix<Real>& m,Real min1,Real max1,Real min2,Real max2,bool weakForm,bool negate)
{
	m.Resize(dim1*dim2);
	for(int i=0;i<dim1*dim2;i++) m.SetGroupSize(i,(2*Degree1+1)*(2*Degree2+1));
	SparseMatrix<Real> dot1,d2Dot1,dot2,d2Dot2;
	FiniteElements1D<Real,Type1,Degree1>::DotProductMatrix(dim1,dot1,min1,max1);
	FiniteElements1D<Real,Type2,Degree2>::DotProductMatrix(dim2,dot2,min2,max2);
	FiniteElements1D<Real,Type1,Degree1>::LaplacianMatrix(dim1,d2Dot1,min1,max1,weakForm,negate);
	FiniteElements1D<Real,Type2,Degree2>::LaplacianMatrix(dim2,d2Dot2,min2,max2,weakForm,negate);

	for(int y=0;y<dim2;y++)
		for(int x=0;x<dim1;x++)
			for(int yy=0;yy<(2*Degree2+1);yy++)
				for(int xx=0;xx<(2*Degree1+1);xx++)
				{
					m.m_ppElements[x+y*dim1][xx+yy*(2*Degree1+1)].N = dot1.m_ppElements[x][xx].N+dot2.m_ppElements[y][yy].N*dim1;
					m.m_ppElements[x+y*dim1][xx+yy*(2*Degree1+1)].Value = 
						(dot1.m_ppElements[x][xx].Value*d2Dot2.m_ppElements[y][yy].Value + d2Dot1.m_ppElements[x][xx].Value*  dot2.m_ppElements[y][yy].Value ) * gWeight +
						(dot1.m_ppElements[x][xx].Value*dot2.m_ppElements[y][yy].Value) * iWeight;
				}
	return true;
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DivergenceStencil(int dim1,int dim2,FullDivergenceStencil& s)
{
	typename FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct<Type1,Degree1>::FullDotProductStencil dot1;
	typename FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct<Type2,Degree2>::FullDotProductStencil dot2;
	typename FiniteElements1D< Real , DERIVATIVE(Type1) , Degree1-1 >::template DotProduct< Type1 , Degree1 >::FullDotProductStencil d1;
	typename FiniteElements1D< Real , DERIVATIVE(Type2) , Degree2-1 >::template DotProduct< Type2 , Degree2 >::FullDotProductStencil d2;

	FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::DotProductStencil( dim1 , dot1 , 0 , 0 );
	FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::DotProductStencil( dim2 , dot2 , 0 , 0 );
	FiniteElements1D< Real , DERIVATIVE(Type1) , Degree1-1 >::template DotProduct< Type1 , Degree1 >::DotProductStencil(FiniteElements1D< Real , DERIVATIVE(Type1) , Degree1-1 >::Dimension( FiniteElements1D< Real , Type1 , Degree1 >::DomainSize( dim1 ) ) , d1 , 0 , 1 );
	FiniteElements1D< Real , DERIVATIVE(Type2) , Degree2-1 >::template DotProduct< Type2 , Degree2 >::DotProductStencil(FiniteElements1D< Real , DERIVATIVE(Type2) , Degree2-1 >::Dimension( FiniteElements1D< Real , Type2 , Degree2 >::DomainSize( dim2 ) ) , d2 , 0 , 1 );

	for(int i=0;i<=2*Degree1;i++)
		for(int j=0;j<=2*Degree2;j++)
		{
			for( int k=0 ; k<2*Degree1 ; k++ ) for( int l=0 ; l<=2*Degree2 ; l++ )
				s.caseTable[i][j].values1[k][l]=d1.caseTable[i].values[k]*dot2.caseTable[j].values[l];
			for( int k=0 ; k<=2*Degree1 ; k++ ) for( int l=0 ; l<2*Degree2 ; l++ )
				s.caseTable[i][j].values2[k][l]=dot1.caseTable[i].values[k]*d2.caseTable[j].values[l];
		}
	return true;
}

template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::BiLaplacianStencil( int dim1 , int dim2 , FullMatrixStencil& s , bool weakForm , bool negate )
{
	typename FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::FullDotProductStencil dot1 , biLap1;
	typename FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::FullDotProductStencil dot2 , biLap2;

	FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::DotProductStencil( dim1 , dot1 , 0 , 0 );
	if( weakForm ) FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::DotProductStencil( dim1 , biLap1 , 2 , 2 , negate );
	else           FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::DotProductStencil( dim1 , biLap1 , 4 , 0 , negate );
	FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::DotProductStencil( dim2 , dot2 , 0 , 0 );
	if( weakForm ) FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::DotProductStencil( dim2 , biLap2 , 2 , 2 , negate );
	else           FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::DotProductStencil( dim2 , biLap2 , 4 , 0 , negate );

	for( int i=0 ; i<=2*Degree1 ; i++ )
		for( int j=0 ; j<=2*Degree2 ; j++ )
			for( int k=0 ; k<=2*Degree1 ; k++ )
				for( int l=0 ; l<=2*Degree2 ; l++ )
					s.caseTable[i][j].values[k][l] =
					dot1.caseTable[i].values[k]*biLap2.caseTable[j].values[l]+
					biLap1.caseTable[i].values[k]*dot2.caseTable[j].values[l];
	return true;
}

template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::LaplacianStencil(int dim1,int dim2,FullMatrixStencil& s,bool weakForm,bool negate)
{
	typename FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::FullDotProductStencil dot1,lap1;
	typename FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::FullDotProductStencil dot2,lap2;

	FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::DotProductStencil(dim1,dot1,0,0);
	if(weakForm)	FiniteElements1D< Real , Type1 , Degree1>::template DotProduct< Type1 , Degree1 >::DotProductStencil(dim1,lap1,1,1,negate);
	else			FiniteElements1D< Real , Type1 , Degree1>::template DotProduct< Type1 , Degree1 >::DotProductStencil(dim1,lap1,2,0,negate);
	FiniteElements1D<Real,Type2,Degree2>::template DotProduct<Type2,Degree2>::DotProductStencil(dim2,dot2,0,0);
	if(weakForm)	FiniteElements1D<Real,Type2,Degree2>::template DotProduct<Type2,Degree2>::DotProductStencil(dim2,lap2,1,1,negate);
	else			FiniteElements1D<Real,Type2,Degree2>::template DotProduct<Type2,Degree2>::DotProductStencil(dim2,lap2,2,0,negate);

	for(int i=0;i<=2*Degree1;i++)
		for(int j=0;j<=2*Degree2;j++)
			for(int k=0;k<=2*Degree1;k++)
				for(int l=0;l<=2*Degree2;l++)
					s.caseTable[i][j].values[k][l]=
					dot1.caseTable[i].values[k]*lap2.caseTable[j].values[l]+
					lap1.caseTable[i].values[k]*dot2.caseTable[j].values[l];
	return true;
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DotProductStencil( int dim1 , int dim2 , FullMatrixStencil& s )
{
	typename FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::FullDotProductStencil dot1 , dot2;

	FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::DotProductStencil( dim1 , dot1 , 0 , 0 );
	FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::DotProductStencil( dim2 , dot2 , 0 , 0 );

	for( int i=0 ; i<=2*Degree1 ; i++ ) for( int j=0 ; j<=2*Degree2 ; j++ )
		for( int k=0 ; k<=2*Degree1 ; k++ ) for( int l=0 ; l<=2*Degree2 ; l++ )
			s.caseTable[i][j].values[k][l] = dot1.caseTable[i].values[k] * dot2.caseTable[j].values[l];
	return true;
}

template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::ProlongationStencil(int lowD1,int lowD2,FullProlongationStencil& s,int& highD1,int& highD2)
{
	typename FiniteElements1D<Real,Type1,Degree1>::FullProlongationStencil stencil1;
	typename FiniteElements1D<Real,Type2,Degree2>::FullProlongationStencil stencil2;

	FiniteElements1D<Real,Type1,Degree1>::ProlongationStencil(lowD1,stencil1,highD1);
	FiniteElements1D<Real,Type1,Degree2>::ProlongationStencil(lowD2,stencil2,highD2);

	for(int i=0;i<=2*Degree1;i++)
		for(int j=0;j<=2*Degree2;j++)
			for(int k=0;k<Degree1+2;k++)
				for(int l=0;l<Degree2+2;l++)
					s.caseTable[i][j].values[k][l]=stencil1.caseTable[i].values[k]*stencil2.caseTable[j].values[l];
	return true;
}
#if NEW_LAPLACIAN_CODE
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::RestrictionStencil(int highD1,int highD2,FullRestrictionStencil& s,int& lowD1,int& lowD2)
{
	typename FiniteElements1D<Real,Type1,Degree1>::FullRestrictionStencil stencil1;
	typename FiniteElements1D<Real,Type2,Degree2>::FullRestrictionStencil stencil2;

	FiniteElements1D<Real,Type1,Degree1>::RestrictionStencil(highD1,stencil1,lowD1);
	FiniteElements1D<Real,Type1,Degree2>::RestrictionStencil(highD2,stencil2,lowD2);

	for(int i=0;i<2*Degree1+2;i++)
		for(int j=0;j<2*Degree2+2;j++)
			for(int k=0;k<(Degree1+3)>>1;k++)
				for(int l=0;l<(Degree2+3)>>1;l++)
					s.caseTable[i][j].values[k][l]=stencil1.caseTable[i].values[k]*stencil2.caseTable[j].values[l];
	return true;
}
#endif // NEW_LAPLACIAN_CODE
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
void FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::StripDiagonal(const SparseMatrix<Real>& M,Vector<MatrixEntry<Real> >& D,bool clearDiagonal)
{
	int idx=((2*Degree1+1)*(2*Degree2+1))/2;
	D.Resize(M.groups);
	for(int i=0;i<M.groups;i++)
	{
		D(i)=M.m_ppElements[i][idx];
		if(clearDiagonal)	M.m_ppElements[i][idx].Value=0;
	}
}

template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::UpSample(Vector<Data>& in,Vector<Data>& out,int lowD1,int lowD2,int& highD1,int& highD2)
{
	SparseMatrix<Real> up1,up2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::UpSampleMatrix(up1,lowD1,highD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::UpSampleMatrix(up2,lowD2,highD2) )
		return false;
	return TensorMatrixMultiply(lowD1,highD1,lowD2,highD2,up1,up2,in,out);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DownSample(Vector<Data>& in,Vector<Data>& out,
																	int highD1,int highD2,int& lowD1,int& lowD2)
{
	SparseMatrix<Real> down1,down2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::DownSampleMatrix(down1,highD1,lowD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::DownSampleMatrix(down2,highD2,lowD2) )
		return false;
	return TensorMatrixMultiply(highD1,lowD1,highD2,lowD2,down1,down2,in,out);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::UpSample(const Data* in,Data* out,int lowD1,int lowD2,int& highD1,int& highD2)
{
	SparseMatrix<Real> up1,up2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::UpSampleMatrix(up1,lowD1,highD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::UpSampleMatrix(up2,lowD2,highD2) )
		return false;
	return TensorMatrixMultiply(lowD1,highD1,lowD2,highD2,up1,up2,in,out);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::UpSample(const Real* in,Real* out,int lowD1,int lowD2,int& highD1,int& highD2,int Channels)
{
	SparseMatrix<Real> up1,up2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::UpSampleMatrix(up1,lowD1,highD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::UpSampleMatrix(up2,lowD2,highD2) )
		return false;
	return TensorMatrixMultiply(lowD1,highD1,lowD2,highD2,up1,up2,in,out,Channels);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DownSample(const Data* in,Data* out,int highD1,int highD2,int& lowD1,int& lowD2)
{
	SparseMatrix<Real> down1,down2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::DownSampleMatrix(down1,highD1,lowD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::DownSampleMatrix(down2,highD2,lowD2) )
		return false;
	return TensorMatrixMultiply(highD1,lowD1,highD2,lowD2,down1,down2,in,out);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DownSample(const Real* in,Real* out,int highD1,int highD2,int& lowD1,int& lowD2,int Channels)
{
	SparseMatrix<Real> down1,down2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::DownSampleMatrix(down1,highD1,lowD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::DownSampleMatrix(down2,highD2,lowD2) )
		return false;
	return TensorMatrixMultiply(highD1,lowD1,highD2,lowD2,down1,down2,in,out,Channels);
}


template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::UpSampleMatrix(SparseMatrix<Real>& M,
																		int lowD1,int lowD2,int& highD1,int& highD2)
{
	SparseMatrix<Real> matrix1,matrix2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::UpSampleMatrix(matrix1,lowD1,highD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::UpSampleMatrix(matrix2,lowD2,highD2) )
		return false;
	return TensorMatrix(lowD1,highD1,matrix1,matrix2,M);
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::DownSampleMatrix(SparseMatrix<Real>& M,
																		  int highD1,int highD2,int& lowD1,int& lowD2)
{
	SparseMatrix<Real> matrix1,matrix2;
	if(	!FiniteElements1D<Real,Type1,Degree1>::DownSampleMatrix(matrix1,highD1,lowD1) ||
		!FiniteElements1D<Real,Type2,Degree2>::DownSampleMatrix(matrix2,highD2,lowD2) )
		return false;
	return TensorMatrix(highD1,lowD1,matrix1,matrix2,M);
}


template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::Gradient(const Vector<Data>& in,Vector<Data>& dX,Vector<Data>& dY,
																  int w,int h,int& ww,int& hh)
{
	SparseMatrix< Real > identity1,identity2,derivative1,derivative2;
	ww = FiniteElements1D< Real , Type1 , Degree1 >::DerivativeMatrix( w , derivative1 );
	hh = FiniteElements1D< Real , Type2 , Degree2 >::DerivativeMatrix( h , derivative2 );
	FiniteElements1D< Real , Type1 , Degree1 >::IdentityMatrix( w , identity1 , derivative1.rowMajor );
	FiniteElements1D< Real , Type2 , Degree2 >::IdentityMatrix( h , identity2 , derivative2.rowMajor );

	return
		FiniteElements2D< Real , Type1 , Degree1 >::TensorMatrixMultiply( w , ww , h ,h , derivative1 , identity2 , in , dX ) &&
		FiniteElements2D< Real , Type2 , Degree2 >::TensorMatrixMultiply( w , w , h , hh , identity1 , derivative2 , in , dY );
}
template<class Real,int Type1,int Degree1,int Type2,int Degree2>
template<class Data>
bool FiniteElements2D<Real,Type1,Degree1,Type2,Degree2>::Laplacian(const Vector<Data>& dX,const Vector<Data>& dY,Vector<Data>& out,
																   int w,int h,int ww,int hh)
{
#if 1
	FullDivergenceStencil stencil;
	DivergenceStencil(w,h,stencil);

	out.Resize(w*h);

	for(int i=0;i<h;i++)
	{
		int ii;
		if(i<Degree2)			ii=i;
		else if(i>=h-Degree2)	ii=2*Degree2+(i-(h-1));
		else					ii=Degree2;
		int dI=i+FiniteElements1D< Real , DERIVATIVE(Type2) , Degree2-1 >::template DotProduct< Type2 , Degree2 >::Helper::StartOffset();
		int I =i+FiniteElements1D< Real , Type2 , Degree2 >::template DotProduct< Type2 , Degree2 >::Helper::StartOffset();
		for(int j=0;j<w;j++)
		{
			Data temp=0;
			int jj;
			if(j<Degree1)			jj=j;
			else if(j>=w-Degree2)	jj=2*Degree1+(j-(w-1));
			else					jj=Degree1;
			int dJ=j+FiniteElements1D< Real , DERIVATIVE(Type1) , Degree1-1 >::template DotProduct< Type1 , Degree1 >::Helper::StartOffset();
			int J =j+FiniteElements1D< Real , Type1 , Degree1 >::template DotProduct< Type1 , Degree1 >::Helper::StartOffset();

			// Partial w.r.t minor index
			for(int di=0;di<2*Degree2;di++)
			{
				if(dI+di<0 || dI+di>=hh)	continue;
				for(int dj=0;dj<=2*Degree1;dj++)
				{
					if(J+dj<0 || J+dj>=w)	continue;
					temp+=dY[(dI+di)*w+(J+dj)]*stencil.caseTable[jj][ii].values2[dj][di];
				}
			}
			// Partial w.r.t major index
			for(int di=0;di<=2*Degree2;di++)
			{
				if(I+di<0 || I+di>=h)	continue;
				for(int dj=0;dj<2*Degree1;dj++)
				{
					if(dJ+dj<0 || dJ+dj>=ww)	continue;
					temp+=dX[(I+di)*ww+(dJ+dj)]*stencil.caseTable[jj][ii].values1[dj][di];
				}
			}
			out[i*w+j]=temp;
		}
	}
#else
	SparseMatrix<Real> dotX,dotY,dDotX,dDotY;
	Vector<Data> tempX,tempY;

	FiniteElements1D<Real,DERIVATIVE(Type1),Degree1-1>::DotProduct<Type1,Degree1>::DerivativeMatrix(ww,dDotX,0,1);
	FiniteElements1D<Real,DERIVATIVE(Type2),Degree2-1>::DotProduct<Type2,Degree2>::DerivativeMatrix(hh,dDotY,0,1);
	FiniteElements1D<Real,Type1,Degree1>::DotProductMatrix(w,dotX);
	FiniteElements1D<Real,Type2,Degree2>::DotProductMatrix(h,dotY);
	if(	!FiniteElements2D<Real,Type1,Degree1>::TensorMatrixMultiply(ww,w,h,h,dDotX,dotY,dX,tempX) ||
		!FiniteElements2D<Real,Type2,Degree2>::TensorMatrixMultiply(w,w,hh,h,dotX,dDotY,dY,tempY) )
		return false;
	out.Resize(tempX.Dimensions());
	for(size_t i=0;i<tempX.Dimensions();i++)	out[i]=tempX[i]+tempY[i];
#endif
	return true;
}
