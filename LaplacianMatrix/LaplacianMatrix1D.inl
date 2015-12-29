#include "FunctionBasis/PPolynomial.h"

inline int ModIndex(int idx,int mod)
{
	if( idx<0 ) return ( mod - ( (-idx) % mod ) ) % mod;
	else        return             idx            % mod;
}
template<class Real>
void StripDiagonal(const SparseMatrix<Real>& M,Vector<MatrixEntry<Real> >& D,bool clearDiagonal,int idx)
{
	D.Resize(M.groups);
	for(int i=0;i<M.groups;i++)
	{
		D(i)=M.m_ppElements[i][idx];
		if(clearDiagonal)	M.m_ppElements[i][idx].Value=0;
	}
}


template<>
inline BinomialCoefficients<0>::BinomialCoefficients(void){coeffs[0]=1;}
template<int Degree>
BinomialCoefficients<Degree>::BinomialCoefficients(void)
{
	BinomialCoefficients<Degree-1> bc;
	memset(coeffs,0,sizeof(int)*(Degree+1));
	for(int i=0;i<Degree;i++)
	{
		coeffs[i  ]+=bc.coeffs[i];
		coeffs[i+1]+=bc.coeffs[i];
	}
}
/////////////////////////
// FiniteDifferences1D //
/////////////////////////
template<class Real,int Type>
void FiniteDifferences1D<Real,Type>::SetValues(int dim,int i,MatrixEntry<Real>* values,bool makePositive)
{
	Real value=Real(1.0);
	if(makePositive)	value=-value;

	// Set the default values
	for(int j=0;j<3;j++)
	{
		values[j].N=i;
		values[j].Value=0;
	}

	switch(Type)
	{
	case PERIODIC_BOUNDARY:
		values[1].N=i;
		values[1].Value=-2*value;
		values[0].N=ModIndex(i-1,dim);
		values[0].Value=value;
		values[2].N=ModIndex(i+1,dim);
		values[2].Value=value;
		break;
	case ZERO_DERIVATIVE:
		values[1].N=i;
		values[1].Value=-2*value;
		if(i==0)
		{
			values[2].N=i+1;
			values[2].Value=2*value;
		}
		else if(i==dim-1)
		{
			values[0].N=i-1;
			values[0].Value=2*value;
		}
		else
		{
			values[0].N=i-1;
			values[0].Value=value;
			values[2].N=i+1;
			values[2].Value=value;
		}
		break;
	case ZERO_VALUE:
		values[1].N=i;
		values[1].Value=-2*value;
		if(i!=0)
		{
			values[0].N=i-1;
			values[0].Value=value;
		}
		if(i!=dim-1)
		{
			values[2].N=i+1;
			values[2].Value=value;
		}
		break;
	}
}

template<class Real,int Type>
void FiniteDifferences1D<Real,Type>::GetMatrix(const int& dim,SparseMatrix<Real>& m,bool makePositive)
{
	m.Resize(dim);
	for(int i=0;i<dim;i++) m.SetGroupSize(i,3);

	MatrixEntry<Real> values[3];
	for(int i=0;i<dim;i++)
	{
		SetValues(dim,i,values,makePositive);

		// The diagonal element
		m.m_ppElements[i][0].N=values[1].N;
		m.m_ppElements[i][0].Value=values[1].Value;

		// The rest
		m.m_ppElements[i][1].N=values[0].N;
		m.m_ppElements[i][1].Value=values[0].Value;
		m.m_ppElements[i][2].N=values[2].N;
		m.m_ppElements[i][2].Value=values[2].Value;
	}
}
template<class Real,int Type>
void FiniteDifferences1D<Real,Type>::StripDiagonal(const SparseMatrix<Real>& M,Vector<MatrixEntry<Real> >& D,bool clearDiagonal)
{
	StripDiagonal(M,D,clearDiagonal,0);
}
//////////////////////
// FiniteElements1D //
//////////////////////
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::IsDownSamplable(int highD,int& lowD)
{
	int domain=DomainSize(highD);
	if(!(domain&1))
	{
		lowD=Dimension(domain>>1);
		return true;
	}
	else return false;
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::IsUpSamplable(int lowD,int& highD)
{
	highD=Dimension(DomainSize(lowD)<<1);
	return true;
}

template<class Real,int Type,int Degree>
int FiniteElements1D<Real,Type,Degree>::DomainSize(int dim)
{
	switch(Type)
	{
	case ZERO_VALUE:
		if(Degree&1)	return dim+1;
		else			return dim  ;
	case ZERO_DERIVATIVE:
		if(Degree&1)	return dim-1;
		else			return dim  ;
	case PERIODIC_BOUNDARY:
		return dim;
	}
	return dim;
}
template<class Real,int Type,int Degree>
int FiniteElements1D<Real,Type,Degree>::Dimension(int dom)
{
	switch(Type)
	{
	case ZERO_VALUE:
		if(Degree&1)	return dom-1;
		else			return dom; 
	case ZERO_DERIVATIVE:
		if(Degree&1)	return dom+1;
		else			return dom  ;
	case PERIODIC_BOUNDARY:
		return dom;
	}
	return dom;
}
template<class Real,int Type,int Degree>
int FiniteElements1D<Real,Type,Degree>::DimensionOffset(void)
{
	switch(Type)
	{
	case ZERO_VALUE:
		if(Degree&1)	return  1;
		else			return  0;
	case ZERO_DERIVATIVE:
		if(Degree&1)	return -1;
		else			return  0;
	case PERIODIC_BOUNDARY:
		return 0;
	}
	return 0;
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::DownSampleMatrix(SparseMatrix<Real>& M,int highD,int& lowD)
{
	Real weights[Degree+2];
	BinomialCoefficients<Degree+1> bc;
	int c=0;
	for(int i=0;i<Degree+2;i++)	c+=bc.coeffs[i];
//	for(int i=0;i<Degree+2;i++)	weights[i]=float(bc.coeffs[i])*2.f/c*sqrt(2.0);
	for(int i=0;i<Degree+2;i++)	weights[i]=float(bc.coeffs[i])*2.f/c;
	return DownSampleMatrix(M,highD,weights,Degree+2,lowD);
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::UpSampleMatrix(SparseMatrix<Real>& M,int highD,int& lowD)
{
	Real weights[Degree+2];
	BinomialCoefficients<Degree+1> bc;
	int c=0;
	for(int i=0;i<Degree+2;i++)	c+=bc.coeffs[i];
//	for(int i=0;i<Degree+2;i++)	weights[i]=float(bc.coeffs[i])*2.f/c*sqrt(2.0);
	for(int i=0;i<Degree+2;i++)	weights[i]=float(bc.coeffs[i])*2.f/c;
	return UpSampleMatrix(M,highD,weights,Degree+2,lowD);
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::UpSampleMatrix(SparseMatrix<Real>& M,int lowD,const Real* weights,int weightCount,int& highD)
{
	if(!IsUpSamplable(lowD,highD))	return false;
	SparseMatrix<Real> MM;
	int start,groupSize=(weightCount+1)>>1;
	if(weightCount&1)	start=-weightCount/2;
	else				start=-weightCount/2+1;

	if(!DownSampleMatrix(MM,highD,weights,weightCount,lowD))	return false;
	M.Resize(highD,true);
	for(int i=0;i<highD;i++)
	{
		M.SetGroupSize(i,groupSize);
		for(int j=0;j<groupSize;j++)
		{
			M.m_ppElements[i][j].N=-1;
			M.m_ppElements[i][j].Value=0;
		}
	}
	for(int i=0;i<lowD;i++)
	{
		for(int j=0;j<weightCount;j++)
			if(MM.m_ppElements[i][j].N>=0)
			{
				int off=MM.m_ppElements[i][j].N-start-weightCount+2;
				if(off<0)	off=(off-1)/2;
				else		off>>=1;
				off=i-off;

				M.m_ppElements[MM.m_ppElements[i][j].N][off].N=i;
				M.m_ppElements[MM.m_ppElements[i][j].N][off].Value=MM.m_ppElements[i][j].Value;
			}
	}
	return true;
}
#if NEW_LAPLACIAN_CODE
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::RestrictionStencil(int highD,FullRestrictionStencil& s,int & lowD)
{
#if 1
	if(!IsDownSamplable(highD,lowD))	return false;
	FullProlongationStencil temp;
	if(!ProlongationStencil(lowD,temp,highD))	return false;
	for(int i=0;i<2*Degree+2;i++)	memset(s.caseTable[i].values,0,sizeof(Real)*FullRestrictionStencil::RestrictionStencil::Size);

	for(int i=0;i<Degree;i++)
	{
		int I,startI;
		I=i;							// The high-res index
		startI=FullRestrictionStencil::RestrictionStencil::Start(I);
		for(int j=0;j<FullRestrictionStencil::RestrictionStencil::Size;j++)
		{
			int J=startI+j;				// The low-res index
			if(J>=0 && J<lowD)
			{
				Real value;
				int off=I-FullProlongationStencil::ProlongationStencil::Start(J);
				if(off<0 || off>=FullProlongationStencil::ProlongationStencil::Size)	continue;
				if(J<Degree)			value=temp.caseTable[J].values[off];
				else if(J>=lowD-Degree)	value=temp.caseTable[2*Degree+(J-(lowD-1))].values[off];
				else					value=temp.caseTable[Degree].values[off];
				s.caseTable[i].values[j]=value;
			}
		}

		I=highD-1-i;					// The high-res index
		startI=FullRestrictionStencil::RestrictionStencil::Start(I);
		for(int j=0;j<FullRestrictionStencil::RestrictionStencil::Size;j++)
		{
			int J=startI+j;				// The low-res index
			if(J>=0 && J<lowD)
			{
				Real value;
				int off=I-FullProlongationStencil::ProlongationStencil::Start(J);
				if(off<0 || off>=FullProlongationStencil::ProlongationStencil::Size)	continue;
				if(J<Degree)			value=temp.caseTable[J].values[off];
				else if(J>=lowD-Degree)	value=temp.caseTable[2*Degree+(J-(lowD-1))].values[off];
				else					value=temp.caseTable[Degree].values[off];
				s.caseTable[2*Degree+1-i].values[j]=value;
			}
		}
	}
	{
		int I,startI;
		I=Degree;							// The high-res index
		startI=FullRestrictionStencil::RestrictionStencil::Start(I);
		for(int j=0;j<FullRestrictionStencil::RestrictionStencil::Size;j++)
		{
			int J=startI+j;					// The low-res index
			{
				int off=I-FullProlongationStencil::ProlongationStencil::Start(J);
				if(off<0 || off>=FullProlongationStencil::ProlongationStencil::Size)	continue;
				s.caseTable[Degree].values[j]=temp.caseTable[Degree].values[off];
			}
		}
		I=Degree+1;							// The high-res index
		startI=FullRestrictionStencil::RestrictionStencil::Start(I);
		for(int j=0;j<FullRestrictionStencil::RestrictionStencil::Size;j++)
		{
			int J=startI+j;					// The low-res index
			{
				int off=I-FullProlongationStencil::ProlongationStencil::Start(J);
				if(off<0 || off>=FullProlongationStencil::ProlongationStencil::Size)	continue;
				s.caseTable[Degree+1].values[j]=temp.caseTable[Degree].values[off];
			}
		}
	}

#else
	int halfD=(Degree+1)>>1;
	int dual= (Degree&1)?0:1;
	int parityBit=((halfD+dual)&1);
	int weightCount=Degree+2;
	if(!IsDownSamplable(highD,lowD))	return false;
	int start,groupSize=(weightCount+1)>>1;
	if(weightCount&1)	start=-weightCount/2;
	else				start=-weightCount/2+1;

	FullProlongationStencil temp;
	if(!ProlongationStencil(lowD,temp,highD))	return false;
	for(int i=0;i<2*Degree+2;i++)	memset(s.caseTable[i].values,0,sizeof(Real)*((Degree+3)>>1));

	for(int i=0;i<Degree;i++)
	{
		int I,startI;
		I=i;							// The high-res index
		startI=(I-halfD-dual+1)>>1;
		for(int j=0;j<groupSize;j++)
		{
			int J=startI+j;				// The low-res index
			if(J>=0 && J<lowD)
			{
				Real value;
				int off=I-(2*J-halfD);
				if(off<0 || off>=Degree+2)	continue;
				if(J<Degree)			value=temp.caseTable[J].values[off];
				else if(J>=lowD-Degree)	value=temp.caseTable[2*Degree+(J-(lowD-1))].values[off];
				else					value=temp.caseTable[Degree].values[off];
				s.caseTable[i].values[j]=value;
			}
		}

		I=highD-1-i;					// The high-res index
		startI=(I-halfD-dual+1)>>1;
		for(int j=0;j<groupSize;j++)
		{
			int J=startI+j;				// The low-res index
			if(J>=0 && J<lowD)
			{
				Real value;
				int off=I-(2*J-halfD);
				if(off<0 || off>=Degree+2)	continue;
				if(J<Degree)			value=temp.caseTable[J].values[off];
				else if(J>=lowD-Degree)	value=temp.caseTable[2*Degree+(J-(lowD-1))].values[off];
				else					value=temp.caseTable[Degree].values[off];
				s.caseTable[2*Degree+1-i].values[j]=value;
			}
		}
	}
	{
		int I,startI;
		I=0;								// The high-res index
		startI=(I-halfD-dual+1)>>1;
		for(int j=0;j<groupSize;j++)
		{
			int J=startI+j;					// The low-res index
			{
				int off=I-(2*J-halfD);
				if(off<0 || off>=Degree+2)	continue;
				s.caseTable[Degree].values[j]=temp.caseTable[Degree].values[off];
			}
		}
		I=1;								// The high-res index
		startI=(I-halfD-dual+1)>>1;
		for(int j=0;j<groupSize;j++)
		{
			int J=startI+j;					// The low-res index
			{
				int off=I-(2*J-halfD);
				if(off<0 || off>=Degree+2)	continue;
				s.caseTable[Degree+1].values[j]=temp.caseTable[Degree].values[off];
			}
		}
	}
#endif
	return true;
}
#endif // NEW_LAPLACIAN_CODE
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::ProlongationStencil(int lowD,FullProlongationStencil& s,int& highD)
{
	Real weights[Degree+2];
	int weightCount=Degree+2;
	BinomialCoefficients<Degree+1> bc;
	int c=0;
	for(int i=0;i<weightCount;i++)	c+=bc.coeffs[i];
	for(int i=0;i<weightCount;i++)	weights[i]=Real(bc.coeffs[i])*2.0/c;

	int start,stop;
	int halfD=(Degree+1)>>1;
	if(!IsUpSamplable(lowD,highD))	return false;

	if(weightCount&1)
	{
		start=-weightCount/2;
		stop = weightCount/2;
	}
	else
	{
		start=-weightCount/2+1;
		stop = weightCount/2;
	}
	for(int i=0;i<=2*Degree;i++)	for(int j=0;j<weightCount;j++)	s.caseTable[i].values[j]=0;
	int NN;
	for(int i=0;i<Degree && i<lowD;i++)
	{
		for(int k=start;k<=stop;k++)
			if(Type==ZERO_DERIVATIVE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				bool refX;

				x1=i;
				if(Degree&1)	stepX=2*(lowD-1);
				else			stepX=2*lowD;
				if(!(Degree&1))	x2=-x1-1;
				else			x2=-x1  ;


#if MISHA_CODE
				refX=true;
#else // !MISHA_CODE
				refX=!(Degree&1) || (ModIndex(i,stepX)!=0 && ModIndex(i,stepX)!=lowD-1);
#endif // MISHA_CODE

				startX1=(int)ceil (Real(       -k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1-k-2*x1)/(2*stepX))*stepX;

				if(refX)
				{
					startX2=(int)ceil (Real(       -k-2*x2)/(2*stepX))*stepX;
					stopX2 =(int)floor(Real(highD-1-k-2*x2)/(2*stepX))*stepX;
				}

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[i].values[2*xx+k-(2*i-halfD)]+=weights[k-start];
				}
				if(refX)
					for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
					{
						NN=ModIndex(2*xx+k-2*i+1,weightCount);
						s.caseTable[i].values[2*xx+k-(2*i-halfD)]+=weights[k-start];
					}
			}
			else if(Type==ZERO_VALUE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				int offset=Degree&1;

				x1=i+offset;
				stepX=2*(lowD+offset);
				x2=-x1-(1-offset);

				startX1=(int)ceil (Real(        offset-k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1+offset-k-2*x1)/(2*stepX))*stepX;
				startX2=(int)ceil (Real(        offset-k-2*x2)/(2*stepX))*stepX;
				stopX2 =(int)floor(Real(highD-1+offset-k-2*x2)/(2*stepX))*stepX;

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[i].values[2*xx+k-offset-(2*i-halfD)]+=weights[k-start];
				}
				for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[i].values[2*xx+k-offset-(2*i-halfD)]-=weights[k-start];
				}
			}
	}
	for(int i=lowD-1;i>=lowD-Degree && i>=0;i--)
	{
		for(int k=start;k<=stop;k++)
			if(Type==ZERO_DERIVATIVE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				bool refX;

				x1=i;
				if(Degree&1)	stepX=2*(lowD-1);
				else			stepX=2*lowD;
				if(!(Degree&1))	x2=-x1-1;
				else			x2=-x1  ;

#if MISHA_CODE
				refX=true;
#else // !MISHA_CODE
				refX=!(Degree&1) || (ModIndex(i,stepX)!=0 && ModIndex(i,stepX)!=lowD-1);
#endif // MISHA_CODE

				startX1=(int)ceil (Real(       -k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1-k-2*x1)/(2*stepX))*stepX;

				if(refX)
				{
					startX2=(int)ceil (Real(       -k-2*x2)/(2*stepX))*stepX;
					stopX2 =(int)floor(Real(highD-1-k-2*x2)/(2*stepX))*stepX;
				}

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[2*Degree-((lowD-1)-i)].values[2*xx+k-(2*i-halfD)]+=weights[k-start];
				}
				if(refX)
					for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
					{
						NN=ModIndex(2*xx+k-2*i+1,weightCount);
						s.caseTable[2*Degree-((lowD-1)-i)].values[2*xx+k-(2*i-halfD)]+=weights[k-start];
					}
			}
			else if(Type==ZERO_VALUE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				int offset=Degree&1;

				x1=i+offset;
				stepX=2*(lowD+offset);
				x2=-x1-(1-offset);

				startX1=(int)ceil (Real(        offset-k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1+offset-k-2*x1)/(2*stepX))*stepX;
				startX2=(int)ceil (Real(        offset-k-2*x2)/(2*stepX))*stepX;
				stopX2 =(int)floor(Real(highD-1+offset-k-2*x2)/(2*stepX))*stepX;

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[2*Degree-((lowD-1)-i)].values[2*xx+k-offset-(2*i-halfD)]+=weights[k-start];
				}
				for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[2*Degree-((lowD-1)-i)].values[2*xx+k-offset-(2*i-halfD)]-=weights[k-start];
				}
			}
	}
	if(lowD>2*Degree)
	{
		int i=Degree;
		for(int k=start;k<=stop;k++)
			if(Type==ZERO_DERIVATIVE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				bool refX;

				x1=i;
				if(Degree&1)	stepX=2*(lowD-1);
				else			stepX=2*lowD;
				if(!(Degree&1))	x2=-x1-1;
				else			x2=-x1  ;

#if MISHA_CODE
				refX=true;
#else // !MISHA_CODE
				refX=!(Degree&1) || (ModIndex(i,stepX)!=0 && ModIndex(i,stepX)!=lowD-1);
#endif // MISHA_CODE

				startX1=(int)ceil (Real(       -k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1-k-2*x1)/(2*stepX))*stepX;

				if(refX)
				{
					startX2=(int)ceil (Real(       -k-2*x2)/(2*stepX))*stepX;
					stopX2 =(int)floor(Real(highD-1-k-2*x2)/(2*stepX))*stepX;
				}

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[i].values[2*xx+k-(2*i-halfD)]+=weights[k-start];
				}
				if(refX)
					for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
					{
						NN=ModIndex(2*xx+k-2*i+1,weightCount);
						s.caseTable[i].values[2*xx+k-(2*i-halfD)]+=weights[k-start];
					}
			}
			else if(Type==ZERO_VALUE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				int offset=Degree&1;

				x1=i+offset;
				stepX=2*(lowD+offset);
				x2=-x1-(1-offset);

				startX1=(int)ceil (Real(        offset-k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1+offset-k-2*x1)/(2*stepX))*stepX;
				startX2=(int)ceil (Real(        offset-k-2*x2)/(2*stepX))*stepX;
				stopX2 =(int)floor(Real(highD-1+offset-k-2*x2)/(2*stepX))*stepX;

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[i].values[2*xx+k-offset-(2*i-halfD)]+=weights[k-start];
				}
				for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					s.caseTable[i].values[2*xx+k-offset-(2*i-halfD)]-=weights[k-start];
				}
			}
	}
	return true;
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::DownSampleMatrix(SparseMatrix<Real>& M,int highD,const Real* weights,int weightCount,int& lowD)
{
	int start,stop;
	if(!IsDownSamplable(highD,lowD))	return false;
	M.Resize(lowD,true);

	if(weightCount&1)
	{
		start=-weightCount/2;
		stop = weightCount/2;
	}
	else
	{
		start=-weightCount/2+1;
		stop = weightCount/2;
	}
	int NN;
	for(int i=0;i<lowD;i++)
	{
		M.SetGroupSize(i,weightCount);
		for(int j=0;j<weightCount;j++)
		{
//			M.m_ppElements[i][j].N=0;
			M.m_ppElements[i][j].N=-1;
			M.m_ppElements[i][j].Value=0;
		}

		for(int k=start;k<=stop;k++)
			if(Type==PERIODIC_BOUNDARY)
			{
				int startX,stopX;
				startX=(int)ceil (Real(       -k-2*i)/(2*lowD))*lowD;
				stopX =(int)floor(Real(highD-1-k-2*i)/(2*lowD))*lowD;
				for(int xx=i+startX;xx<=i+stopX;xx+=lowD)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					M.m_ppElements[i][NN].N=2*xx+k;
					M.m_ppElements[i][NN].Value+=weights[k-start];
				}
			}
			else if(Type==ZERO_DERIVATIVE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				bool refX;

				x1=i;
				if(Degree&1)	stepX=2*(lowD-1);
				else			stepX=2*lowD;
				if(!(Degree&1))	x2=-x1-1;
				else			x2=-x1  ;

#if MISHA_CODE
				refX=true;
#else // !MISHA_CODE
				refX=!(Degree&1) || (ModIndex(i,stepX)!=0 && ModIndex(i,stepX)!=lowD-1);
#endif // MISHA_CODE

				startX1=(int)ceil (Real(       -k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1-k-2*x1)/(2*stepX))*stepX;

				if(refX)
				{
					startX2=(int)ceil (Real(       -k-2*x2)/(2*stepX))*stepX;
					stopX2 =(int)floor(Real(highD-1-k-2*x2)/(2*stepX))*stepX;
				}

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					M.m_ppElements[i][NN].N=2*xx+k;
					M.m_ppElements[i][NN].Value+=weights[k-start];
				}
				if(refX)
					for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
					{
						NN=ModIndex(2*xx+k-2*i+1,weightCount);
						M.m_ppElements[i][NN].N=2*xx+k;
						M.m_ppElements[i][NN].Value+=weights[k-start];
					}
			}
			else if(Type==ZERO_VALUE)
			{
				int stepX,startX1,stopX1,x1,startX2,stopX2,x2;
				int offset=Degree&1;

				x1=i+offset;
				stepX=2*(lowD+offset);
				x2=-x1-(1-offset);

				startX1=(int)ceil (Real(        offset-k-2*x1)/(2*stepX))*stepX;
				stopX1 =(int)floor(Real(highD-1+offset-k-2*x1)/(2*stepX))*stepX;
				startX2=(int)ceil (Real(        offset-k-2*x2)/(2*stepX))*stepX;
				stopX2 =(int)floor(Real(highD-1+offset-k-2*x2)/(2*stepX))*stepX;

				for(int xx=x1+startX1;xx<=x1+stopX1;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					M.m_ppElements[i][NN].N=2*xx+k-offset;
					M.m_ppElements[i][NN].Value+=weights[k-start];
				}
				for(int xx=x2+startX2;xx<=x2+stopX2;xx+=stepX)
				{
					NN=ModIndex(2*xx+k-2*i+1,weightCount);
					M.m_ppElements[i][NN].N=2*xx+k-offset;
					M.m_ppElements[i][NN].Value-=weights[k-start];
				}
			}
	}

	return true;
}
template<class Real,int Type,int Degree>
int FiniteElements1D<Real,Type,Degree>::DerivativeMatrix(int dim,SparseMatrix<Real>& m)
{
	int outDim;
	outDim=FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(DomainSize(dim));
	int offset=Degree&1;

	m.Resize(dim,false);
	for(int i=0;i<dim;i++)
	{
		m.SetGroupSize(i,2);
		for(int j=0;j<2;j++)
		{
			m.m_ppElements[i][j].N=0;
			m.m_ppElements[i][j].Value=0;
		}
	}
	switch(Type)
	{
	case PERIODIC_BOUNDARY:
		for(int i=0;i<dim;i++)
		{
			m.m_ppElements[i][0].N=ModIndex(i-offset,outDim);
			m.m_ppElements[i][0].Value+=1;
			m.m_ppElements[i][1].N=ModIndex(i+1-offset,outDim);
			m.m_ppElements[i][1].Value-=1;
		}
		break;
	case ZERO_DERIVATIVE:
		for(int i=1;i<dim-1;i++)
		{
			m.m_ppElements[i][0].N=i-1;
			m.m_ppElements[i][0].Value+=1;
			m.m_ppElements[i][1].N=i;
			m.m_ppElements[i][1].Value-=1;
		}
		m.m_ppElements[0][1].N=0;
		m.m_ppElements[0][1].Value-=1;
		m.m_ppElements[dim-1][0].N=dim-2;
		m.m_ppElements[dim-1][0].Value+=1;
		break;
	case ZERO_VALUE:
		for(int i=0;i<dim;i++)
		{
			m.m_ppElements[i][0].N=i;
			m.m_ppElements[i][0].Value+=1;
			m.m_ppElements[i][1].N=i+1;
			m.m_ppElements[i][1].Value-=1;
		}
		if(!(Degree&1))
		{
			m.m_ppElements[0][0].Value+=1;
			m.m_ppElements[dim-1][1].Value-=1;
		}
		break;
	}
	return outDim;
}

template<class Real,int Type,int Degree>
template<class Data>
bool FiniteElements1D<Real,Type,Degree>::UpSample(Vector<Data>& in,Vector<Data>& out)
{
	Real weights[Degree+2];
	BinomialCoefficients<Degree+1> bc;
	int c=0;
	for(int i=0;i<Degree+2;i++)	c+=bc.coeffs[i];
	for(int i=0;i<Degree+2;i++)	weights[i]=float(bc.coeffs[i])*2.f/c;
	return UpSample(in,out,weights,Degree+2);
}
template<class Real,int Type,int Degree>
template<class Data>
bool FiniteElements1D<Real,Type,Degree>::DownSample(Vector<Data>& in,Vector<Data>& out)
{
	Real weights[Degree+2];
	BinomialCoefficients<Degree+1> bc;
	int c=0;
	for(int i=0;i<Degree+2;i++)	c+=bc.coeffs[i];
	for(int i=0;i<Degree+2;i++)	weights[i]=float(bc.coeffs[i])*2.f/c;
	return DownSample(in,out,weights,Degree+2);
}
template<class Real,int Type,int Degree>
template<class Data>
bool FiniteElements1D<Real,Type,Degree>::DownSample(Vector<Data>& in,Vector<Data>& out,
											   const Real* weights,int weightCount)
{
	int outDim;
	if(!IsDownSamplable(in.Dimensions(),outDim))	return false;
	SparseMatrix<Real> downMatrix;
	if(!DownSampleMatrix(downMatrix,in.Dimensions(),weights,weightCount,outDim) )	return false;
	out=downMatrix*in;
	return true;
}
template<class Real,int Type,int Degree>
template<class Data>
bool FiniteElements1D<Real,Type,Degree>::UpSample(Vector<Data>& in,Vector<Data>& out,																						
												  const Real* weights,int weightCount)
{
	int outDim;
	if(!IsUpSamplable(in.Dimensions(),outDim))	return false;

	SparseMatrix<Real> upMatrix;
	if(!UpSampleMatrix(upMatrix,in.Dimensions(),weights,weightCount,outDim) )	return false;
	out.Resize(outDim);
	upMatrix.Multiply(in,out);
	return true;
}

template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::LaplacianMatrix(int dim,SparseMatrix<Real>& m,bool weakForm,bool negate)
{
	return LaplacianMatrix(dim,m,0,Real(DomainSize(dim)),weakForm,negate);
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::LaplacianMatrix(int dim,SparseMatrix<Real>& m,Real min,Real max,bool weakForm,bool negate)
{
	if(weakForm)	return FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::DerivativeMatrix(dim,m,min,max,1,1,negate);
	else			return FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::DerivativeMatrix(dim,m,min,max,2,0,negate);
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::IdentityMatrix(int dim,SparseMatrix<Real>& m,bool rowMajor)
{
	m.Resize(dim,rowMajor);
	for(int i=0;i<dim;i++)
	{
		m.SetGroupSize(i,1);
		m.m_ppElements[i][0].N=i;
		m.m_ppElements[i][0].Value=1;
	}
	return true;
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::DotProductMatrix(int dim,SparseMatrix<Real>& m)
{
	return DotProductMatrix(dim,m,0,Real(DomainSize(dim)));
}
template<class Real,int Type,int Degree>
bool FiniteElements1D<Real,Type,Degree>::DotProductMatrix(int dim,SparseMatrix<Real>& m,Real min,Real max)
{
	return FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::DerivativeMatrix(dim,m,min,max,0,0);
}


template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
bool FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::DotProductStencil( int dim , FullDotProductStencil& s , int d1 , int d2 , bool negate )
{
	int dotDim=FiniteElements1D<Real,DotType,DotDegree>::Dimension(DomainSize(dim));

	Helper helper;
	helper.setDerivatives( d1 , d2 );
	MatrixEntry<Real> values[Degree+DotDegree+1];
	for( int i=0 ; i<=2*DotDegree ; i++ ) memset( s.caseTable[i].values , 0 , sizeof(Real)*(Degree+DotDegree+1) );

	for( int i=0 ; i<dotDim && i<DotDegree ; i++ )
	{
		helper.SetValues( dotDim , i , values , Real(0) , Real( DomainSize(dim) ) );
		for( int k=0 ; k<=Degree+DotDegree ; k++ )
		{
			if(negate)	s.caseTable[i].values[k] = -values[k].Value;
			else		s.caseTable[i].values[k] =  values[k].Value;
		}
	}
	for(int i=dotDim-1;i>=dotDim-DotDegree && i>=0;i--)
	{
		helper.SetValues( dotDim , i , values , Real(0) , Real( DomainSize(dim) ) );
		for(int k=0;k<=Degree+DotDegree;k++)
		{
			if(negate)	s.caseTable[2*DotDegree-(dotDim-1-i)].values[k]=-values[k].Value;
			else		s.caseTable[2*DotDegree-(dotDim-1-i)].values[k]= values[k].Value;
		}
	}
	if(dotDim>2*DotDegree)
	{
		helper.SetValues( dotDim , DotDegree , values , Real(0) , Real( DomainSize(dim) ) );
		for(int k=0;k<=Degree+DotDegree;k++)
		{
			if(negate)	s.caseTable[DotDegree].values[k]=-values[k].Value;
			else		s.caseTable[DotDegree].values[k]= values[k].Value;
		}
	}
	return true;
}
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
bool FiniteElements1D< Real , Type , Degree >::DotProduct<DotType,DotDegree>::PaddedDotProductStencil( int dim , FullPaddedDotProductStencil& s , int d1 , int d2 , bool negate )
{
	FullDotProductStencil _s;
	if( !DotProductStencil( dim , _s , d1 , d2 , negate ) ) return false;
	for( int i=0 ; i<DotDegree ; i++ ) for(int j=0 ; j<Degree+DotDegree+1 ; j++ ) s.caseTable[i].values[j] = _s.caseTable[i].values[j];
	for( int i=0 ; i<DotDegree ; i++ ) for(int j=0 ; j<Degree+DotDegree+1 ; j++ ) s.caseTable[2*DotDegree+2-i].values[j] = _s.caseTable[2*DotDegree-i].values[j];
	for( int i=0 ; i<3 ; i++ ) for(int j=0 ; j<Degree+DotDegree+1 ; j++ ) s.caseTable[DotDegree+i].values[j] = _s.caseTable[DotDegree].values[j];
	return true;
}

template<class Real,int Type,int Degree>
void FiniteElements1D< Real , Type , Degree >::StripDiagonal(const SparseMatrix<Real>& M,Vector<MatrixEntry<Real> >& D,bool clearDiagonal)
{
	int idx=(2*Degree+1)/2;
	D.Resize(M.groups);
	for(int i=0;i<M.groups;i++)
	{
		D(i)=M.m_ppElements[i][idx];
		if(clearDiagonal)	M.m_ppElements[i][idx].Value=0;
	}
}

//////////////////////////////////////////
// FiniteElements1D::DotProduct::Helper //
//////////////////////////////////////////

// Set the helper up so that it facilitates the computation of the
// matrix entries for a row-major representation.
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
const int FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::Start=
((Degree&1)==(DotDegree&1))?-(Degree+DotDegree)/2:((Degree&1)?-(Degree+DotDegree)/2:-(Degree+DotDegree+1)/2);
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
const int FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::Stop=
((Degree&1)==(DotDegree&1))? (Degree+DotDegree)/2:((Degree&1)? (Degree+DotDegree+1)/2: (Degree+DotDegree)/2);

template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::Helper(void)
{
	F1=PPolynomial<Degree   >::GaussianApproximation();
	F2=PPolynomial<DotDegree>::GaussianApproximation();
}
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
void FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::setDerivatives(int d1,int d2)
{
	double min1,min2,max1,max2;
	shift=dotShift=0;

	if(!(Degree&1) )	shift=0.5;
	if(!(DotDegree&1) )	dotShift=0.5;
	min1=-(Degree+1)*0.5+shift;
	max1= (Degree+1)*0.5+shift;
	min2=-(DotDegree+1)*0.5+dotShift;
	max2= (DotDegree+1)*0.5+dotShift;

	F1=F1.shift(shift);
	F2=F2.shift(dotShift);
	for(int i=0;i<d1;i++)	F1=((PPolynomial<Degree   +1>)F1).derivative();
	for(int i=0;i<d2;i++)	F2=((PPolynomial<DotDegree+1>)F2).derivative();

	for(int d=Start;d<=Stop;d++)
	{
		F1F2[d-Start]=(F1.shift(d)*F2);
		min[d-Start]=(min1+d)>min2?(min1+d):min2;
		max[d-Start]=(max1+d)<max2?(max1+d):max2;
		fullValues[d-Start]=F1F2[d-Start].integral(min[d-Start],max[d-Start]);
	}
}
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
double FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::GetValue(int i,int j,double min,double max) const
{
	double lower,upper;
	int d=i-j;

	if(d<Start || d>Stop)	return 0;

	if(	(i+shift   -(Degree   +1)*0.5>min && i+shift   +(Degree   +1)*0.5<max) ||
		(j+dotShift-(DotDegree+1)*0.5>min && j+dotShift+(DotDegree+1)*0.5<max))
		return fullValues[d-Start];
	else
	{
		lower = min-j > this->min[d-Start] ? min-j : this->min[d-Start];
		upper = max-j < this->max[d-Start] ? max-j : this->max[d-Start];
		if(lower<upper)	return F1F2[d-Start].integral(lower,upper);
	}
	return 0;
}

template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
int FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::StartOffset(void)
{
	int i=Start;
	if		(Type==ZERO_VALUE && (Degree&1) && DotType==ZERO_DERIVATIVE)	i--;
	else if	(DotType==ZERO_VALUE && (DotDegree&1) && Type==ZERO_DERIVATIVE)	i++;
	if(Type==ZERO_VALUE || Type==ZERO_DERIVATIVE)	return i;
	else
	{
		fprintf(stderr,"Boundary type is not supported\n");
		return -1;
	}
}
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
int FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::StopOffset(void)
{
	return Degree+DotDegree+StartOffset();
}
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
void FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::Helper::SetValues(int dotDim,int j,MatrixEntry<Real>* values,Real min,Real max)
{
	bool center[]={!(Degree&1),!(DotDegree&1)};
	int dim[]={Dimension(FiniteElements1D<Real,DotType,DotDegree>::DomainSize(dotDim)),dotDim};
	int deg[]={Degree,DotDegree};
	int type[]={Type,DotType};
	Real _shift[]={Real(shift),Real(dotShift)};
	for(int d=Start;d<=Stop;d++)
	{
		values[d-Start].N=0;
		values[d-Start].Value=0;
	}
	for(int d=Start;d<=Stop;d++)
	{
		int i=j+d;
		if		(Type==ZERO_VALUE && (Degree&1) && DotType==ZERO_DERIVATIVE)	i--;
		else if	(DotType==ZERO_VALUE && (DotDegree&1) && Type==ZERO_DERIVATIVE)	i++;
		if(Type==ZERO_VALUE || Type==ZERO_DERIVATIVE)
			if(i<0 || i>=dim[0])	continue;
			else					values[d-Start].N=i;
		else
			if(d-Start>=dim[0])		continue;
			else					values[d-Start].N=ModIndex(i,dim[0]);

		Real _min[2],_max[2];
		int _step[2],_start[2],_stop[2],_x[2],x[2],_refStart[2],_refStop[2],_refX[2];
		int _ref[2];

		_x[0]=i;
		_x[1]=j;

		for(int k=0;k<2;k++)
		{
			switch(type[k])
			{
			case ZERO_VALUE:
				_step[k]=2*DomainSize(dim[0]);
				_ref[k]=-1;
				if(!center[k])	_x[k]++;
				break;
			case ZERO_DERIVATIVE:
				_step[k]=2*DomainSize(dim[0]);
#if MISHA_CODE
				_ref[k]=1;
#else // !MISHA_CODE
				if(center[k] || (ModIndex(_x[k],_step[k])!=0 && ModIndex(_x[k],_step[k])!=dim[k]-1))	_ref[k]=1;
				else																					_ref[k]=0;
#endif // MISHA_CODE
				break;
			case PERIODIC_BOUNDARY:
				_step[k]=  DomainSize(dim[0]);
				_ref[k]=0;
				break;
			};
			if(center[k])	_refX[k]=-_x[k]-1;
			else			_refX[k]=-_x[k];
			_min[k]=Real( min-_shift[k]-0.5*(deg[k]+1) );
			_max[k]=Real( max-_shift[k]+0.5*(deg[k]+1) );
			_start[k]=(int)ceil (Real(_min[k]-_x[k])/_step[k])*_step[k];
			_stop[k] =(int)floor(Real(_max[k]-_x[k])/_step[k])*_step[k];
			_refStart[k]=(int)ceil (Real(_min[k]-_refX[k])/_step[k])*_step[k];
			_refStop[k] =(int)floor(Real(_max[k]-_refX[k])/_step[k])*_step[k];
		}

		for(x[0]=_x[0]+_start[0];x[0]<=_x[0]+_stop[0];x[0]+=_step[0])
		{
			for(x[1]=_x[1]+_start[1];x[1]<=_x[1]+_stop[1];x[1]+=_step[1])
				values[d-Start].Value += Real( GetValue(x[0],x[1],min,max) );
			if(_ref[1])
				for(x[1]=_refX[1]+_refStart[1];x[1]<=_refX[1]+_refStop[1];x[1]+=_step[1])
					values[d-Start].Value += Real( GetValue(x[0],x[1],min,max)*_ref[1] );
		}
		if(_ref[0])
			for(x[0]=_refX[0]+_refStart[0];x[0]<=_refX[0]+_refStop[0];x[0]+=_step[0])
			{
				for(x[1]=_x[1]+_start[1];x[1]<=_x[1]+_stop[1];x[1]+=_step[1])
					values[d-Start].Value += Real( GetValue(x[0],x[1],min,max)*_ref[0] );
				if(_ref[1])
					for(x[1]=_refX[1]+_refStart[1];x[1]<=_refX[1]+_refStop[1];x[1]+=_step[1])
						values[d-Start].Value += Real( GetValue(x[0],x[1],min,max)*_ref[0]*_ref[1] );
			}
	}
}
///////////////////////////////////////////////
// FiniteElements1D::FullProlongationStencil //
///////////////////////////////////////////////
template<class Real,int Type,int Degree>
void FiniteElements1D<Real,Type,Degree>::FullProlongationStencil::Normalize(int lowD)
{
	int highD;
	int halfD=(Degree+1)>>1;
	IsUpSamplable(lowD,highD);
	for(int ii=0;ii<=2*Degree;ii++)
	{
		int i;
		if(ii<=Degree)	i=ii;
		else			i=lowD-1-(2*Degree-ii);
		int startI=2*i-halfD;
		for(int j=0;j<Degree+2;j++)	if(startI+j<0 || startI+j>=highD)	caseTable[ii].values[j]=0;
	}
}
//////////////////////////////////
// FiniteElements1D::DotProduct //
//////////////////////////////////
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
bool FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::DerivativeMatrix(int dim,SparseMatrix<Real>& m,int d1,int d2,bool negate)
{
	return DerivativeMatrix(dim,m,0,Real(DomainSize(dim)),d1,d2,negate);
}
template<class Real,int Type,int Degree>
template<int DotType,int DotDegree>
bool FiniteElements1D<Real,Type,Degree>::DotProduct<DotType,DotDegree>::DerivativeMatrix(int dim,SparseMatrix<Real>& m,Real min,Real max,int d1,int d2,bool negate)
{
	int dotDim=FiniteElements1D<Real,DotType,DotDegree>::Dimension(DomainSize(dim));
	Helper helper;
	helper.setDerivatives(d1,d2);

	m.Resize(dotDim);
	for(int i=0;i<dotDim;i++) m.SetGroupSize(i,Degree+DotDegree+1);

	MatrixEntry<Real> values[Degree+DotDegree+1];
	for(int i=0;i<dotDim;i++)
	{
		helper.SetValues(dotDim,i,values,min,max);
		for(int k=0;k<=Degree+DotDegree;k++)
		{
			m.m_ppElements[i][k].N=values[k].N;
			if(negate)	m.m_ppElements[i][k].Value=-values[k].Value;
			else		m.m_ppElements[i][k].Value= values[k].Value;
		}
	}
	return true;
}
template<class Real,int Type,int Degree>
void CombineStencils( const typename FiniteElements1D< Real , Type , Degree >::template DotProduct< Type , Degree >::FullDotProductStencil& dStencil ,
					  const typename FiniteElements1D< Real , Type , Degree >::FullProlongationStencil& pStencil , int dim ,
					  typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil)
{
	int lowDim;
	FiniteElements1D< Real , Type , Degree >::IsDownSamplable(dim,lowDim);
	CombineStencils< Real , Type , Degree >(dStencil,pStencil,lowDim,dim,newDStencil);
}
template<class Real,int Type,int Degree>
void CombineStencils( const typename FiniteElements1D< Real , Type , Degree >::template DotProduct< Type , Degree >::FullDotProductStencil& dStencil ,
					  const typename FiniteElements1D< Real , Type , Degree >::FullProlongationStencil& pStencil ,
					  int lowDim , int dim ,
					  typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil )
{
	int halfD=(Degree+1)>>1;

	// Clear the output stencil
	for( int ii=0 ; ii<=2*Degree ; ii++ ) for( int jj=0 ; jj<=2*Degree ; jj++ ) newDStencil.caseTable[ii].values[jj] = 0;

	for( int ii=0 ; ii<=2*Degree ; ii++ )
	{
		int i;											// Index in the low resolution
		if( ii<=Degree)	i=ii;
		else			i=(lowDim-1)-(2*Degree-ii);

		int startI=2*i-halfD;
		int endI=startI+Degree+1;
		for(int j=startI;j<=endI;j++)					// Index in the high resolution
			if(j>=0 && j<dim)
			{
				int jj;
				if		(j<Degree)			jj=j;
				else if	(j>dim-1-Degree)	jj=2*Degree+(j-(dim-1));
				else						jj=Degree;
				for(int k=j-Degree;k<=j+Degree;k++)		// Anything that overlaps in the high resolution
					if(k>=0 && k<dim)
					{
						int startK=(k+halfD-Degree)>>1;
						int endK  =(k+halfD       )>>1;

						for(int l=startK;l<=endK;l++)	// Anything that overlaps in the lower resolution
							if(l>=0 && l<lowDim)
							{
								int startL=2*l-halfD;
								int ll;
								if		(l<Degree)			ll=l;
								else if	(l>lowDim-1-Degree)	ll=2*Degree+(l-(lowDim-1));
								else						ll=Degree;

								// Sanity Check
								if(k<startL || k>startL+Degree+1)	fprintf(stderr,"Badness 1\n");

								// Sanity Check
								if(i-l<-Degree || i-l>Degree)	fprintf(stderr,"Badness 2\n");

								newDStencil.caseTable[ii].values[l-i+Degree]+=
									pStencil.caseTable[ll].values[k-startL]*
									dStencil.caseTable[jj].values[k-j+Degree]*
									pStencil.caseTable[ii].values[j-startI];
							}
					}
			}
	}
}
template<class Real,int Type,int Degree>
void CombineStencils( const typename FiniteElements1D< Real , Type , Degree >::template DotProduct< Type , Degree >::FullDotProductStencil& dStencil,
					  const typename FiniteElements1D< Real , Type , Degree-1 >::FullProlongationStencil& pStencil , int dim ,
					  typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil)
{
	int LDegree=Degree-1;
	int halfD=(LDegree+1)>>1;
	int lowDim;
	FiniteElements1D<Real,Type,Degree-1>::IsDownSamplable(dim,lowDim);

	for(int ii=0;ii<=2*Degree;ii++)	for(int jj=0;jj<=2*Degree;jj++)	newDStencil.caseTable[ii].values[jj]=0;

	for(int c=0;c<=2*Degree;c++)
	{
		int i;											// Index in the low resolution
		if(c<=Degree)	i=c;
		else			i=(lowDim-1)-(2*Degree-c);

		int ii;
		if		(i<LDegree)			ii=i;
		else if	(i>lowDim-1-LDegree)ii=2*LDegree+(i-(lowDim-1));
		else						ii=LDegree;

		int startI=2*i-halfD;
		int endI=startI+LDegree+1;
		for(int j=startI;j<=endI;j++)					// Index in the high resolution
			if(j>=0 && j<dim)
			{
				int jj;
				if		(j<Degree)			jj=j;
				else if	(j>dim-1-Degree)	jj=2*Degree+(j-(dim-1));
				else						jj=Degree;
				for(int k=j-Degree;k<=j+Degree;k++)		// Anything that overlaps in the high resolution
					if(k>=0 && k<dim)
					{
						int startK=(k+halfD-LDegree)>>1;
						int endK=(k+halfD)>>1;

						for(int l=startK;l<=endK;l++)	// Anything that overlaps in the lower resolution
							if(l>=0 && l<lowDim)
							{
								int startL=2*l-halfD;
								int ll;
								if		(l<LDegree)			ll=l;
								else if	(l>lowDim-1-LDegree)ll=2*LDegree+(l-(lowDim-1));
								else						ll=LDegree;

								// Sanity Check
								if(k<startL || k>startL+LDegree+1)	fprintf(stderr,"Badness 1\n");

								// Sanity Check
								if(i-l<-Degree || i-l>Degree)	fprintf(stderr,"Badness 2\n");

								newDStencil.caseTable[c].values[l-i+Degree]+=
									pStencil.caseTable[ll].values[k-startL]*
									dStencil.caseTable[jj].values[k-j+Degree]*
									pStencil.caseTable[ii].values[j-startI];
							}
					}
			}
	}
}
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree+1>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullDotProductStencil& newDStencil)
{
	int LDegree=Degree+1;
	int halfD=(LDegree+1)>>1;
	int lowDim;
	FiniteElements1D<Real,Type,Degree+1>::IsDownSamplable(dim,lowDim);

	for(int ii=0;ii<=2*Degree;ii++)	for(int jj=0;jj<=2*Degree;jj++)	newDStencil.caseTable[ii].values[jj]=0;

	for(int c=0;c<=2*Degree;c++)
	{
		int i;											// Index in the low resolution
		if(c<=Degree)	i=c;
		else			i=(lowDim-1)-(2*Degree-c);

		int ii;
		if		(i<LDegree)			ii=i;
		else if	(i>lowDim-1-LDegree)ii=2*LDegree+(i-(lowDim-1));
		else						ii=LDegree;

		int startI=2*i-halfD;
		int endI=startI+LDegree+1;
		for(int j=startI;j<=endI;j++)					// Index in the high resolution
			if(j>=0 && j<dim)
			{
				int jj;
				if		(j<Degree)			jj=j;
				else if	(j>dim-1-Degree)	jj=2*Degree+(j-(dim-1));
				else						jj=Degree;
				for(int k=j-Degree;k<=j+Degree;k++)		// Anything that overlaps in the high resolution
					if(k>=0 && k<dim)
					{
						int startK=(k+halfD-LDegree)>>1;
						int endK=(k+halfD)>>1;

						for(int l=startK;l<=endK;l++)	// Anything that overlaps in the lower resolution
							if(l>=0 && l<lowDim)
							{
								int startL=2*l-halfD;
								int ll;
								if		(l<LDegree)			ll=l;
								else if	(l>lowDim-1-LDegree)ll=2*LDegree+(l-(lowDim-1));
								else						ll=LDegree;

								// Sanity Check
								if		(k<startL || k>startL+LDegree+1)	fprintf(stderr,"Badness 1\n");
								else if	(i-l<-Degree || i-l>Degree)			fprintf(stderr,"Badness 2\n");
								else
									newDStencil.caseTable[c].values[l-i+Degree]+=
										pStencil.caseTable[ll].values[k-startL]*
										dStencil.caseTable[jj].values[k-j+Degree]*
										pStencil.caseTable[ii].values[j-startI];
							}
					}
			}
	}
}
///////////////////////
template<class Real,int Type,int Degree>
void CombineStencils( const typename FiniteElements1D< Real , Type , Degree >::template DotProduct< Type , Degree >::FullPaddedDotProductStencil& dStencil,
					  const typename FiniteElements1D< Real , Type , Degree >::FullProlongationStencil& pStencil,int dim,
					  typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil)
{
	int lowDim;
	FiniteElements1D<Real,Type,Degree>::IsDownSamplable(dim,lowDim);
	CombineStencils<Real,Type,Degree>(dStencil,pStencil,lowDim,dim,newDStencil);
}
template<class Real,int Type,int Degree>
void CombineStencils( const typename FiniteElements1D< Real , Type , Degree >::template DotProduct< Type , Degree >::FullPaddedDotProductStencil& dStencil ,
					  const typename FiniteElements1D< Real , Type , Degree >::FullProlongationStencil& pStencil ,
					  int lowDim , int dim ,
					  typename FiniteElements1D<Real,Type,Degree>::template DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil )
{
	int halfD=(Degree+1)>>1;

	// Clear the output stencil
	for( int ii=0 ; ii<=2*Degree+2 ; ii++ ) for( int jj=0 ; jj<=2*Degree ; jj++ ) newDStencil.caseTable[ii].values[jj] = 0;

	for( int dI=0 ; dI<=2*Degree+2 ; dI++ )
	{
		int i;											// Index in the low resolution
		if( dI<=Degree+1 )	i = dI;
		else				i = (lowDim-1)-(2*Degree+2-dI);
		int pI;
		if		( i<Degree )			pI = i;
		else if	( i>lowDim-1-Degree )	pI = 2*Degree-((lowDim-1)-i);
		else							pI = Degree;

		int startI = 2*i-halfD;
		int endI = startI+Degree+1;
		for( int j=startI ; j<=endI ; j++ )				// Index in the high resolution
			if( j>=0 && j<dim )
			{
				int dJ;
				if		( j<Degree+1 )			dJ = j;
				else if	( j>dim-1-Degree-1 )	dJ = 2*Degree+2-((dim-1)-j);
				else							dJ = Degree+1;
				for( int k=j-Degree ; k<=j+Degree ; k++ )	// Anything that overlaps in the high resolution
					if( k>=0 && k<dim )
					{
						int startK = (k+halfD-Degree)>>1;
						int endK   = (k+halfD       )>>1;

						for( int l=startK ; l<=endK ; l++ )	// Anything that overlaps in the lower resolution
							if( l>=0 && l<lowDim )
							{
								int startL = 2*l-halfD;
								int pL;
								if		( l<Degree )			pL = l;
								else if	( l>lowDim-1-Degree )	pL = 2*Degree-((lowDim-1)-l);
								else							pL = Degree;

								// Sanity Check
								if(k<startL || k>startL+Degree+1)	fprintf(stderr,"Badness 1\n");

								// Sanity Check
								if(i-l<-Degree || i-l>Degree)	fprintf(stderr,"Badness 2\n");

								newDStencil.caseTable[dI].values[l-i+Degree]+=
									pStencil.caseTable[pL].values[k-startL]*
									dStencil.caseTable[dJ].values[k-j+Degree]*
									pStencil.caseTable[pI].values[j-startI];
							}
					}
			}
	}
}
#if 0
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::FullPaddedDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree-1>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil)
{
	int LDegree=Degree-1;
	int halfD=(LDegree+1)>>1;
	int lowDim;
	FiniteElements1D<Real,Type,Degree-1>::IsDownSamplable(dim,lowDim);

	for(int ii=0;ii<=2*Degree+2;ii++)	for(int jj=0;jj<=2*Degree+2;jj++)	newDStencil.caseTable[ii].values[jj]=0;

	for(int c=0;c<=2*Degree+2;c++)
	{
		int i;											// Index in the low resolution
		if(c<=Degree+1)	i=c;
		else			i=(lowDim-1)-(2*Degree+2-c);

		int ii;
		if		(i<LDegree+1)			ii = i;
		else if	(i>lowDim-1-LDegree-1)	ii = 2*LDegree+2+(i-(lowDim-1));
		else							ii = LDegree+1;

		int startI=2*i-halfD;
		int endI=startI+LDegree+1;
		for(int j=startI;j<=endI;j++)					// Index in the high resolution
			if(j>=0 && j<dim)
			{
				int jj;
				if		(j<Degree+1)		jj = j;
				else if	(j>dim-1-Degree-1)	jj = 2*Degree+2+(j-(dim-1));
				else						jj = Degree+1;
				for(int k=j-Degree;k<=j+Degree;k++)		// Anything that overlaps in the high resolution
					if(k>=0 && k<dim)
					{
						int startK=(k+halfD-LDegree)>>1;
						int endK=(k+halfD)>>1;

						for(int l=startK;l<=endK;l++)	// Anything that overlaps in the lower resolution
							if(l>=0 && l<lowDim)
							{
								int startL=2*l-halfD;
								int ll;
								if		(l<LDegree+1)			ll = l;
								else if	(l>lowDim-1-LDegree-1)	ll = 2*LDegree+2+(l-(lowDim-1));
								else							ll = LDegree+1;

								// Sanity Check
								if(k<startL || k>startL+LDegree+1)	fprintf(stderr,"Badness 1\n");

								// Sanity Check
								if(i-l<-Degree || i-l>Degree)	fprintf(stderr,"Badness 2\n");

								newDStencil.caseTable[c].values[l-i+Degree]+=
									pStencil.caseTable[ll].values[k-startL]*
									dStencil.caseTable[jj].values[k-j+Degree]*
									pStencil.caseTable[ii].values[j-startI];
							}
					}
			}
	}
}
template<class Real,int Type,int Degree>
void CombineStencils(const typename FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::FullPaddedDotProductStencil& dStencil,
					 const typename FiniteElements1D<Real,Type,Degree+1>::FullProlongationStencil& pStencil,int dim,
					 typename FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::FullPaddedDotProductStencil& newDStencil)
{
	int LDegree=Degree+1;
	int halfD=(LDegree+1)>>1;
	int lowDim;
	FiniteElements1D<Real,Type,Degree+1>::IsDownSamplable(dim,lowDim);

	for(int ii=0;ii<=2*Degree+2;ii++)	for(int jj=0;jj<=2*Degree+2;jj++)	newDStencil.caseTable[ii].values[jj]=0;

	for(int c=0;c<=2*Degree+2;c++)
	{
		int i;											// Index in the low resolution
		if(c<=Degree+1)	i=c;
		else			i=(lowDim-1)-(2*Degree+2-c);

		int ii;
		if		(i<LDegree+1)			ii = i;
		else if	(i>lowDim-1-LDegree-1)	ii = 2*LDegree+2+(i-(lowDim-1));
		else							ii = LDegree+1;

		int startI=2*i-halfD;
		int endI=startI+LDegree+1;
		for(int j=startI;j<=endI;j++)					// Index in the high resolution
			if(j>=0 && j<dim)
			{
				int jj;
				if		(j<Degree+1)		jj = j;
				else if	(j>dim-1-Degree-1)	jj = 2*Degree+2+(j-(dim-1));
				else						jj = Degree+1;
				for(int k=j-Degree;k<=j+Degree;k++)		// Anything that overlaps in the high resolution
					if(k>=0 && k<dim)
					{
						int startK=(k+halfD-LDegree)>>1;
						int endK=(k+halfD)>>1;

						for(int l=startK;l<=endK;l++)	// Anything that overlaps in the lower resolution
							if(l>=0 && l<lowDim)
							{
								int startL=2*l-halfD;
								int ll;
								if		(l<LDegree+1)			ll = l;
								else if	(l>lowDim-1-LDegree-1)	ll = 2*LDegree+2+(l-(lowDim-1));
								else							ll = LDegree+1;

								// Sanity Check
								if		(k<startL || k>startL+LDegree+1)	fprintf(stderr,"Badness 1\n");
								else if	(i-l<-Degree || i-l>Degree)			fprintf(stderr,"Badness 2\n");
								else
									newDStencil.caseTable[c].values[l-i+Degree]+=
										pStencil.caseTable[ll].values[k-startL]*
										dStencil.caseTable[jj].values[k-j+Degree]*
										pStencil.caseTable[ii].values[j-startI];
							}
					}
			}
	}
}
#endif