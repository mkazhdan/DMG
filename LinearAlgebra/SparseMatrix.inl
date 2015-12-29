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

#include <float.h>


///////////////////
//  SparseMatrix //
///////////////////
///////////////////////////////////////
// SparseMatrix Methods and Memebers //
///////////////////////////////////////

template<class T>
SparseMatrix<T>::SparseMatrix()
{
	rowMajor=true;
	groups=0;
	groupSizes=NULL;
	m_ppElements=NULL;
}
template< class T >       MatrixEntry< T >* SparseMatrix< T >::operator[]( int idx )       { return m_ppElements[idx]; }
template< class T > const MatrixEntry< T >* SparseMatrix< T >::operator[]( int idx ) const { return m_ppElements[idx]; }

template<class T>
bool SparseMatrix<T>::copyColumnMajor(SparseMatrix<T>& out,int inDim,int outDim)
{
	if(rowMajor && outDim!=groups)	return false;
	if(!rowMajor && inDim!=groups)	return false;
	out.Resize(inDim,false);
	if(rowMajor)
	{
		for(int i=0;i<groups;i++)
			for(int j=0;j<groupSizes[i];j++)
				if(0<=m_ppElements[i][j].N && m_ppElements[i][j].N<inDim)
					out.groupSizes[m_ppElements[i][j].N]++;

		for(int i=0;i<out.groups;i++)
		{
			int sz=out.groupSizes[i];
			out.groupSizes[i]=0;
			out.SetGroupSize(i,sz);
			out.groupSizes[i]=0;
		}
		for(int i=0;i<groups;i++)
			for(int j=0;j<groupSizes[i];j++)
			{
				int idx=m_ppElements[i][j].N;
				if(0<=idx && idx<inDim)
				{
					out.m_ppElements[idx][out.groupSizes[idx]].N=i;
					out.m_ppElements[idx][out.groupSizes[idx]].Value=m_ppElements[i][j].Value;
					out.groupSizes[idx]++;
				}
			}
	}
	else
		for(int i=0;i<groups;i++)
		{
			out.SetGroupSize(i,groupSizes[i]);
			memcpy(out.m_ppElements[i],m_ppElements[i],groupSizes[i]*sizeof(MatrixEntry<T>));
		}
	return true;
}
template<class T>
bool SparseMatrix<T>::copyRowMajor(SparseMatrix<T>& out,int inDim,int outDim)
{
	if(rowMajor && outDim!=groups)	return false;
	if(!rowMajor && inDim!=groups)	return false;

	out.Resize(outDim,true);
	if(rowMajor)
		for(int i=0;i<groups;i++)
		{
			out.SetGroupSize(i,groupSizes[i]);
			memcpy(out.m_ppElements[i],m_ppElements[i],groupSizes[i]*sizeof(MatrixEntry<T>));
		}
	else
	{
		for(int i=0;i<groups;i++)
			for(int j=0;j<groupSizes[i];j++)
				if(0<=m_ppElements[i][j].N && m_ppElements[i][j].N<outDim)
					out.groupSizes[m_ppElements[i][j].N]++;
		for(int i=0;i<out.groups;i++)
		{
			int sz=out.groupSizes[i];
			out.groupSizes[i]=0;
			out.SetGroupSize(i,sz);
			out.groupSizes[i]=0;
		}
		for(int i=0;i<groups;i++)
			for(int j=0;j<groupSizes[i];j++)
			{
				int idx=m_ppElements[i][j].N;
				if(0<=idx && idx<outDim)
				{
					out.m_ppElements[idx][out.groupSizes[idx]].N=i;
					out.m_ppElements[idx][out.groupSizes[idx]].Value=m_ppElements[i][j].Value;
					out.groupSizes[idx]++;
				}
			}
	}
	return true;
}

template<class T>
SparseMatrix<T>::SparseMatrix( int groups , bool rowMajor )
{
	rowMajor=true;
	this->groups=0;
	groupSizes=NULL;
	m_ppElements=NULL;
	Resize(groups,rowMajor);
}

template<class T>
SparseMatrix<T>::SparseMatrix( const SparseMatrix& M )
{
	rowMajor=true;
	groups=0;
	groupSizes=NULL;
	m_ppElements=NULL;
	Resize(M.groups,M.rowMajor);
	for (int i=0; i<groups; i++)
	{
		SetGroupSize(i,M.groupSizes[i]);
		for(int j=0;j<groupSizes[i];j++){m_ppElements[i][j]=M.m_ppElements[i][j];}
	}
}
template<class T>
int SparseMatrix<T>::Entries(void){
	int e=0;
	for(int i=0;i<groups;i++){e+=(int)(groupSizes[i]);}
	return e;
}
template<class T>
SparseMatrix<T>& SparseMatrix<T>::operator = (const SparseMatrix<T>& M)
{
	Resize(M.groups,M.rowMajor);
	for (int i=0; i<groups; i++){
		SetGroupSize(i,M.groupSizes[i]);
		for(int j=0;j<groupSizes[i];j++){m_ppElements[i][j]=M.m_ppElements[i][j];}
	}
	return *this;
}

template<class T>
void SparseMatrix<T>::Transpose(void){rowMajor=!rowMajor;}

template<class T>
SparseMatrix<T>::~SparseMatrix()
{
	Resize( 0 , true );
}

template<class T>
void SparseMatrix<T>::Resize( int g, bool rm )
{
	rowMajor = rm;
	if( groups>0 )
	{
		for( int i=0 ; i<groups ; i++ ) if( groupSizes[i] ) free( m_ppElements[i] );
		free( m_ppElements );
		free( groupSizes );
		groupSizes = NULL;
		m_ppElements = NULL;
	}
	groups = g;
	if( g )
	{
		groupSizes = ( int* )malloc( sizeof( int ) * g );
		if( !groupSizes ) fprintf( stderr , "Failed to allocate (1) in SparseMatrix< >::Resize( )\n" );
		memset( groupSizes , 0 , sizeof( int ) * g );
		m_ppElements = ( MatrixEntry<T>** )malloc( sizeof( MatrixEntry<T>* ) * g );
		if( !m_ppElements ) fprintf( stderr , "Failed to allocate (2) in SparseMatrix< >::Resize( )\n" );
	}
}

template<class T>
void SparseMatrix<T>::SetGroupSize(int group,int count)
{
	if(group>=0 && group<groups)
	{
		if( groupSizes[group] ) free(m_ppElements[group]);
		if( count>0 )
		{
			m_ppElements[group] = ( MatrixEntry<T>* )malloc( sizeof( MatrixEntry<T> ) * count );
			if( !m_ppElements[group] ) fprintf( stderr , "Failed to allocated in SparseMatrix< >::SetGroupSize( )\n");
		}
		groupSizes[group] = count;
	}
}
template<class T>
template<class T2>
void SparseMatrix<T>::Multiply	( const Vector<T2>& In,Vector<T2>& Out) const
{
	Out.SetZero();
	if(rowMajor)
		for (int i=0; i<groups; i++)
		{
			T2 temp=T2();
			for(int j=0;j<groupSizes[i];j++)
				temp+=In.m_pV[m_ppElements[i][j].N] * m_ppElements[i][j].Value;
//			Out.m_pV[i]=temp;
			Out.m_pV[i] += temp;
		}
	else
		for ( int i=0 ; i<groups ; i++ )
			for( int j=0 ; j<groupSizes[i] ; j++ )
				Out.m_pV[m_ppElements[i][j].N] += In.m_pV[i] * m_ppElements[i][j].Value;
}
template<class T>
template<class T2>
bool SparseMatrix<T>::Multiply ( const Vector<T2>& in,Vector<T2>& out,int startRow,int stopRow) const
{
	if(!rowMajor)	return false;
	for (int i=startRow; i<stopRow; i++)
	{
		T2 temp=T2();
		for(int j=0;j<groupSizes[i];j++)	temp+=in.m_pV[m_ppElements[i][j].N] * m_ppElements[i][j].Value;
		out.m_pV[i-startRow]=temp;
	}
	return true;
}
template<class T>
template<class T2>
void SparseMatrix<T>::MultiplyTranspose( const Vector<T2>& In,Vector<T2>& Out) const
{
	if(!rowMajor)
		for (int i=0; i<groups; i++)
		{
			T2 temp=T2();
			for(int j=0;j<groupSizes[i];j++)
				if(m_ppElements[i][j].N>=0)
					temp+=In.m_pV[m_ppElements[i][j].N] * m_ppElements[i][j].Value;
			Out.m_pV[i]=temp;
		}
	else
		for (int i=0; i<groups; i++)
			for(int j=0;j<groupSizes[i];j++)
				if(m_ppElements[i][j].N>=0)
					Out.m_pV[m_ppElements[i][j].N]+=In.m_pV[i] * m_ppElements[i][j].Value;
}
template<class T>
template<class T2>
bool SparseMatrix<T>::Multiply ( const SparseMatrix<T>& M, const Vector<MatrixEntry<T> >& D,const Vector<T2>& In,Vector<T2>& Out,int startRow,int stopRow)
{
	if(!M.rowMajor)	return false;
	for (int i=startRow; i<stopRow; i++)
	{
		T2 temp=T2();
		for(int j=0;j<M.groupSizes[i];j++)	temp+=In.m_pV[M.m_ppElements[i][j].N] * M.m_ppElements[i][j].Value;
		Out.m_pV[i-startRow]=temp+In.m_pV[D[i].N]*D[i].Value;
	}
	return true;
}
template<class T>
template<class T2>
bool SparseMatrix<T>::Multiply ( const SparseMatrix<T>& M, const Vector<MatrixEntry2<T> >& D,const Vector<T2>& In,Vector<T2>& Out,int startRow,int stopRow)
{
	if(!M.rowMajor)	return false;
	for (int i=startRow; i<stopRow; i++)
	{
		T2 temp=T2();
		for(int j=0;j<M.groupSizes[i];j++)
			temp+=In.m_pV[M.m_ppElements[i][j].N] * M.m_ppElements[i][j].Value;
		Out.m_pV[i-startRow]=temp+In.m_pV[D[i].inN]*D[i].value;
	}
	return true;
}

template<class T>
template<class T2>
Vector<T2> SparseMatrix<T>::operator * (const Vector<T2>& V) const
{
	if(!rowMajor)
	{
		fprintf(stderr,"Unknown output matrix size for column-major matrix\n");
		exit(0);
	}
	Vector<T2> R( groups );
	for (int i=0; i<groups; i++)
	{
		T2 temp=T2();
		for(int ii=0;ii<groupSizes[i];ii++)	
			if(0<=m_ppElements[i][ii].N && m_ppElements[i][ii].N<(int)(V.Dimensions() ) )
				temp+=V.m_pV[m_ppElements[i][ii].N] * m_ppElements[i][ii].Value;
//		R(i)=temp;
		R(i) += temp;
	}
	return R;
}
template<class T>
template<class T2>
int SparseMatrix<T>::SolveJacobi(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& diagonal, const Vector<T2>& b,int iters,Vector<T2>& solution,bool reset)
{
	Vector<T2> Md;

	if(reset)
	{
		solution.Resize(b.Dimensions());
		solution.SetZero();
	}

	Md.Resize(M.groups);
	for(int i=0;i<iters;i++)
	{
		M.Multiply(solution,Md);
		for(int j=0;j<(int)(diagonal.Dimensions()); j++)
		{
			int idx = diagonal[j].N;
//			solution[idx] = (b[idx]-Md[j]) / diagonal[j].Value;
			solution[idx] = (b[idx]-Md[idx]) / diagonal[j].Value;
		}
	}
	return iters;
}

template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel( const SparseMatrix<T>& M , const Vector<MatrixEntry<T> >& diagonal , const Vector<T2>& b , int iters , Vector<T2>& solution , bool reverseIndex , bool zeroDiagonal )
{
	return SolveGaussSeidel( M , diagonal , b , iters , solution , 0 , (int)(diagonal.Dimensions()) , reverseIndex , zeroDiagonal );
}
template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel( const SparseMatrix<T>& M , const Vector<MatrixEntry<T> >& diagonal , const Vector<T2>& b , int iters , Vector<T2>& solution , int start , int stop , bool reverseIndex , bool zeroDiagonal )
{
	if( !M.rowMajor ) return 0;
	for( int i=0 ; i<iters ; i++ )
	{
		if( reverseIndex )
			if( zeroDiagonal )
				for(int j=stop-1; j>=start ; j--)
				{
					T2 temp;
					temp *= 0;
					int idx = diagonal[j].N;
					if(idx<0)	continue;
#if 1
					for(int k=0;k<M.groupSizes[j];k++)	temp+=solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
					solution[idx] = (b[j]-temp) / diagonal[j].Value;
#else
					for(int k=0;k<M.groupSizes[idx];k++)	temp+=solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
					solution[idx] = (b[idx]-temp) / diagonal[j].Value;
#endif
				}
			else
				for(int j=stop-1; j>=start ; j--)
				{
					T2 temp;
					temp *= 0;
					int idx = diagonal[j].N;
					if(idx<0)	continue;
					for(int k=0;k<M.groupSizes[j];k++)	temp += solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
					temp -= solution[idx] * diagonal[j].Value;
					solution[idx] = (b[j]-temp) / diagonal[j].Value;
				}
		else
			if( zeroDiagonal )
				for(int j=start ; j<stop ; j++)
				{
					T2 temp;
					temp *= 0;
					int idx = diagonal[j].N;
					if( idx<0 ) continue;
#if 1
					for( int k=0 ; k<M.groupSizes[j] ; k++) temp += solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
					solution[idx] = (b[j]-temp) / diagonal[j].Value;
#else
					for(int k=0;k<M.groupSizes[idx];k++)	temp+=solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
					solution[idx] = (b[idx]-temp) / diagonal[j].Value;
#endif
				}
			else
				for(int j=start ; j<stop ; j++)
				{
					T2 temp;
					temp *= 0;
					int idx = diagonal[j].N;
					if( idx<0 ) continue;
					for( int k=0 ; k<M.groupSizes[j] ; k++) temp += solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
					temp -= solution[ idx ] * diagonal[j].Value;
					solution[idx] = (b[j]-temp) / diagonal[j].Value;
				}
	}
	return iters;
}
template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& diagonal, const Vector<T2>& b , int iters , const Vector<T2>& sorWeights , Vector<T2>& solution,bool reverseIndex)
{
	return SolveGaussSeidel(M,diagonal,b,iters,sorWeights,solution,0,(int)(diagonal.Dimensions()),reverseIndex);
}

template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& diagonal, const Vector<T2>& b , int iters , const Vector<T2>& sorWeights , Vector<T2>& solution,int start,int stop,bool reverseIndex)
{
	if(!M.rowMajor)	return 0;
	for(int i=0;i<iters;i++)
	{
		if(reverseIndex)
			for(int j=stop-1; j>=start ; j--)
			{
				T2 temp=0;
				int idx = diagonal[j].N;
				if(idx<0)	continue;
				for(int k=0;k<M.groupSizes[j];k++)	temp += solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
				solution[idx] = solution[idx] +  ( (b[j]-temp) / diagonal[j].Value - solution[idx] ) * sorWeights[idx];
			}
		else
			for(int j=start ; j<stop ; j++)
			{
				T2 temp=0;
				int idx = diagonal[j].N;
				if(idx<0)	continue;
				for(int k=0;k<M.groupSizes[j];k++)	temp+=solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
				solution[idx] = solution[idx] +  ( (b[j]-temp) / diagonal[j].Value - solution[idx] ) * sorWeights[idx];
			}
//printf( "%d %f %f %f\n" , i , solution[0] , solution[1]  , solution[65] );
	}
	return iters;
}


template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel( const SparseMatrix<T>& M , const Vector<MatrixEntry2<T> >& diagonal , const Vector<T2>& b , int iters , Vector<T2>& solution , bool reverseIndex )
{
	return SolveGaussSeidel(M,diagonal,b,iters,solution,0,diagonal.Dimensions(),reverseIndex);
}
template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel( const SparseMatrix<T>& M , const Vector<MatrixEntry2<T> >& diagonal , const Vector<T2>& b , int iters , Vector<T2>& solution , int start , int stop , bool reverseIndex )
{
	if( !M.rowMajor )	return 0;

	for( int i=0 ; i<iters ; i++ )
		if( reverseIndex )
			for(int j=stop-1; j>=start ; j--)
			{
				T2 temp=0;
				if(diagonal[j].outN<0)	continue;
				for(int k=0;k<M.groupSizes[j];k++)	temp+=solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
				solution[diagonal[j].inN] = (b[diagonal[j].outN]-temp) / diagonal[j].value;
			}
		else
			for( int j=start ; j<stop ; j++ )
			{
				T2 temp = 0;
				if( diagonal[j].outN<0 )	continue;
				for( int k=0;k<M.groupSizes[j];k++)	temp+=solution[M.m_ppElements[j][k].N]*M.m_ppElements[j][k].Value;
				solution[diagonal[j].inN] = (b[diagonal[j].outN]-temp) / diagonal[j].value;
			}
	return iters;
}

template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel2(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& diagonal, const Vector<T2>& b,int iters,Vector<T2>& solution,bool reverseIndex)
{
	return SolveGaussSeidel2(M,diagonal,b,iters,solution,0,diagonal.Dimensions(),reverseIndex);
}
template<class T>
template<class T2>
int SparseMatrix<T>::SolveGaussSeidel2(const SparseMatrix<T>& M,const Vector<MatrixEntry<T> >& diagonal, const Vector<T2>& b,int iters,Vector<T2>& solution,int start,int stop,bool reverseIndex)
{
	int count=0;
	if(!M.rowMajor)	return 0;
	if(start>=stop)	return 0;

	for(int i=0;i<iters;i++)
		if(reverseIndex)
			for(int j=stop-1; j>=start ; j--)
			{
				T2 temp=0;
				int idx = diagonal[j].N;
				if(idx<0)	continue;
				for(int k=0;k<M.groupSizes[idx];k++)	temp+=solution[M.m_ppElements[idx][k].N]*M.m_ppElements[idx][k].Value;
				solution[idx] = (b[idx]-temp) / diagonal[j].Value;
				count++;
			}
		else
			for(int j=start;j<stop; j++)
			{
				T2 temp=0;
				int idx = diagonal[j].N;
				if(idx<0)	continue;
				for(int k=0;k<M.groupSizes[idx];k++)	temp+=solution[M.m_ppElements[idx][k].N]*M.m_ppElements[idx][k].Value;
				solution[idx] = (b[idx]-temp) / diagonal[j].Value;
				count++;
			}
	return count;
}

template <class Matrix,class Data>
static int SolveConjugateGradient( const Matrix& SPD , const Vector<Data>& b , const int& iters , Vector<Data>& solution , const double eps , bool printError )
{
	Vector<Data> d,r,Md,temp;
	double alpha,beta,rDotR,oldRDotR;
	Md.Resize(b.Dimensions());

	temp.Resize(b.Dimensions());
	SPD.Multiply(solution,temp);
	d=r=b-temp;
	oldRDotR=rDotR=r.Dot(r);
	if(b.Dot(b)<=eps)
	{
//		printf("Badness1: %g %g\n",r.Dot(r),eps);
		solution.SetZero();
		return 0;
	}
	int i;
	for(i=0;i<iters;i++)
	{
		double temp;
		SPD.Multiply(d,Md);
		temp=d.Dot(Md);
		if(temp<=eps)
		{
//			printf("Badness1: %g %g\n",temp,eps);
			break;
		}
		alpha=rDotR/temp;
		r.SubtractScaled(Md,alpha);
		temp=r.Dot(r);
		if( printError ) printf( "Error[%d] = %f\n" , i , temp );
		// BADNESS!!! How can the size of the residual increase?
		if(temp>2*oldRDotR)
		{
//			printf("Badness2.5: %g %g\n",temp,oldRDotR);
			break;
		}
		oldRDotR=rDotR;
		if(temp/b.Dot(b)<=eps)
		{
//			printf("Badness2: %g %g\n",temp,eps);
			break;
		}
		beta=temp/rDotR;
		solution.AddScaled(d,alpha);
		if(beta<=eps)
		{
//			printf("Badness3: %g %g\n",beta,eps);
			break;
		}
		rDotR=temp;
		Vector<Data>::Add(d,beta,r,d);
	}
	return i;
}

//////////////////////////
// PivotSymmetricMatrix //
//////////////////////////
template <class Matrix>
template <class Data>
Vector<Data> PivotSymmetricMatrix<Matrix>::operator * (const Vector<Data>& in) const
{
	Vector<Data> out;
	Multiply(in,out);
	return out;
}
template <class Matrix>
template <class Data>
void PivotSymmetricMatrix<Matrix>::Multiply(const Vector<Data>& in,Vector<Data>& out) const
{
	Vector<Data> tempIn,tempOut;
	tempIn=in;
	for(int i=rightMatrices.size()-1;i>=0;i--)
	{
		tempOut.Resize(rightMatrices[i].second.second);
		rightMatrices[i].first.Multiply(tempIn,tempOut);
		tempIn=tempOut;
	}
	tempOut.Resize(pivot.second.second);
	pivot.first.Multiply(tempIn,tempOut);
	tempIn=tempOut;

	for(int i=0;i<rightMatrices.size();i++)
	{
		tempOut.Resize(rightMatrices[i].second.first);
		rightMatrices[i].first.MultiplyTranspose(tempIn,tempOut);
		tempIn=tempOut;
	}
	out=tempIn;
}
template <class Matrix>
void PivotSymmetricMatrix<Matrix>::setPivot(const Matrix& M,int inDim,int outDim)
{
	pivot=std::pair<Matrix,std::pair<int,int> >(M,std::pair<int,int>(inDim,outDim));
}
template <class Matrix>
void PivotSymmetricMatrix<Matrix>::push(const Matrix& M,int inDim,int outDim)
{
	rightMatrices.push_back(std::pair<Matrix,std::pair<int,int> >(M,std::pair<int,int>(inDim,outDim)));
}
template <class Matrix>
void PivotSymmetricMatrix<Matrix>::pop(void)
{
	rightMatrices.pop_back();
}
