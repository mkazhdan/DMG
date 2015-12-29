#define ScaledSum128d(x1,x2,x3,x4) (_mm_add_pd(_mm_mul_pd(x1,x2),_mm_mul_pd(x3,x4)))

template<class Real,int Type,int Degree,class Data>
inline double StreamingSolver<Real,Type,Degree,Data>::GetLaplacianValue0(const __m128d matrixValues[][TemplateCount],int j)
{
	ALIGN( const __m128d *mValues , SSE_ALIGNMENT );
	ALIGN( const __m128d* xPtr , SSE_ALIGNMENT );
	ALIGN( double scratch[2] , SSE_ALIGNMENT );
	register __m128d temp1,temp2,sum;

	double temp=0;
	switch(Degree)
	{
#if TEST_LAPLACIAN_DEFAULT
		case 1000:
#else // !TEST_LAPLACIAN_DEFAULT
		case 2:
#endif // TEST_LAPLACIAN_DEFAULT
			{
				int jj=(-Degree+j  )/2;
				// xx=0
				mValues=matrixValues[0];
				xPtr=&localXPtr[0][jj];
				temp1=ScaledSum128d(xPtr[0],mValues[0],xPtr[1],mValues[1]);
				temp2=_mm_mul_pd(xPtr[2],mValues[2]);
				sum=_mm_add_pd(temp1,temp2);
				// xx=1
				mValues=matrixValues[1];
				xPtr=&localXPtr[1][jj];
				temp1=ScaledSum128d(xPtr[0],mValues[0],xPtr[1],mValues[1]);
				temp2=_mm_mul_pd(xPtr[2],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
				// xx=2
				mValues=matrixValues[2];
				xPtr=&localXPtr[2][jj];
				temp1=ScaledSum128d(xPtr[0],mValues[0],xPtr[1],mValues[1]);
				temp2=_mm_mul_pd(xPtr[2],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
				// xx=3
				mValues=matrixValues[3];
				xPtr=&localXPtr[3][jj];
				temp1=ScaledSum128d(xPtr[0],mValues[0],xPtr[1],mValues[1]);
				temp2=_mm_mul_pd(xPtr[2],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
				// xx=4
				mValues=matrixValues[4];
				xPtr=&localXPtr[4][jj];
				temp1=ScaledSum128d(xPtr[0],mValues[0],xPtr[1],mValues[1]);
				temp2=_mm_mul_pd(xPtr[2],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
			}
			_mm_store_pd(scratch,sum);
			return scratch[0]+scratch[1];
#if TEST_LAPLACIAN_DEFAULT
		case 1001:
#else // !TEST_LAPLACIAN_DEFAULT
		case 1:
#endif // TEST_LAPLACIAN_DEFAULT
			{
				int jj=(-Degree+j+1)/2;
				// xx=0
				mValues=matrixValues[0];
				xPtr=&localXPtr[0][jj];
				sum=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				// xx=1
				mValues=matrixValues[1];
				xPtr=&localXPtr[1][jj];
				temp1=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				sum=_mm_add_pd(sum,temp1);
				// xx=2
				mValues=matrixValues[2];
				xPtr=&localXPtr[2][jj];
				temp1=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				sum=_mm_add_pd(sum,temp1);
			}
			_mm_store_pd(scratch,sum);
			return scratch[0]+scratch[1];
		default:
			for(int xx=0;xx<=2*Degree;xx++)
			{
				if(Degree&1)
				{
					mValues=&matrixValues[xx][1];
					xPtr=&localXPtr[xx][(-Degree+j+1)/2];
					temp2=_mm_mul_pd(xPtr[-1],mValues[-1]);
					for(int i=0;i<Degree;i++)
					{
						temp1=_mm_mul_pd(xPtr[i],mValues[i]);
						temp2=_mm_add_pd(temp1,temp2);
					}
					_mm_store_pd(scratch,temp2);
					temp+=scratch[0]+scratch[1];
				}
				else
				{
					mValues=matrixValues[xx];
					xPtr=&localXPtr[xx][(-Degree+j)/2];
					temp2=_mm_mul_pd(xPtr[Degree],mValues[Degree]);
					for(int i=0;i<Degree;i++)
					{
						temp1=_mm_mul_pd(xPtr[i],mValues[i]);
						temp2=_mm_add_pd(temp1,temp2);
					}
					_mm_store_pd(scratch,temp2);
					temp+=scratch[0]+scratch[1];
				}
			}
			return temp;
	}
}
template<class Real,int Type,int Degree,class Data>
inline double StreamingSolver<Real,Type,Degree,Data>::GetLaplacianValue1(const __m128d matrixValues[][TemplateCount],int j)
{
	ALIGN( const __m128d *mValues , SSE_ALIGNMENT );
	ALIGN( const __m128d* xPtr , SSE_ALIGNMENT );
	ALIGN( double scratch[2] , SSE_ALIGNMENT );
	register __m128d temp1,temp2,sum;

	double temp=0;
	switch(Degree)
	{
#if TEST_LAPLACIAN_DEFAULT
		case 1000:
#else // !TEST_LAPLACIAN_DEFAULT
		case 2:
#endif // TEST_LAPLACIAN_DEFAULT
			{
				int jj=(-Degree+j+1)/2;
				// xx=0
				mValues=matrixValues[0];
				xPtr=&localXPtr[0][jj];
				temp1=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				temp2=_mm_mul_pd(xPtr[1],mValues[2]);
				sum=_mm_add_pd(temp1,temp2);
				// xx=1
				mValues=matrixValues[1];
				xPtr=&localXPtr[1][jj];
				temp1=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				temp2=_mm_mul_pd(xPtr[1],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
				// xx=2
				mValues=matrixValues[2];
				xPtr=&localXPtr[2][jj];
				temp1=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				temp2=_mm_mul_pd(xPtr[1],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
				// xx=3
				mValues=matrixValues[3];
				xPtr=&localXPtr[3][jj];
				temp1=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				temp2=_mm_mul_pd(xPtr[1],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
				// xx=4
				mValues=matrixValues[4];
				xPtr=&localXPtr[4][jj];
				temp1=ScaledSum128d(xPtr[-1],mValues[0],xPtr[0],mValues[1]);
				temp2=_mm_mul_pd(xPtr[1],mValues[2]);
				temp1=_mm_add_pd(temp1,temp2);
				sum=_mm_add_pd(sum,temp1);
			}
			_mm_store_pd(scratch,sum);
			return scratch[0]+scratch[1];
#if TEST_LAPLACIAN_DEFAULT
		case 1001:
#else // !TEST_LAPLACIAN_DEFAULT
		case 1:
#endif // TEST_LAPLACIAN_DEFAULT
			{
				int jj=(-Degree+j)/2;
				// xx=0
				mValues=matrixValues[0];
				xPtr=&localXPtr[0][jj];
				sum=ScaledSum128d(xPtr[1],mValues[1],xPtr[0],mValues[0]);
				// xx=1
				mValues=matrixValues[1];
				xPtr=&localXPtr[1][jj];
				temp1=ScaledSum128d(xPtr[1],mValues[1],xPtr[0],mValues[0]);
				sum=_mm_add_pd(sum,temp1);
				// xx=2
				mValues=matrixValues[2];
				xPtr=&localXPtr[2][jj];
				temp1=ScaledSum128d(xPtr[1],mValues[1],xPtr[0],mValues[0]);
				sum=_mm_add_pd(sum,temp1);
			}
			_mm_store_pd(scratch,sum);
			return scratch[0]+scratch[1];
		default:
			for(int xx=0;xx<=2*Degree;xx++)
			{
				if(Degree&1)
				{
					mValues=matrixValues[xx];
					xPtr=&localXPtr[xx][(-Degree+j)/2];
					temp2=_mm_mul_pd(xPtr[Degree],mValues[Degree]);
					for(int i=0;i<Degree;i++)
					{
						temp1=_mm_mul_pd(xPtr[i],mValues[i]);
						temp2=_mm_add_pd(temp1,temp2);
					}
					_mm_store_pd(scratch,temp2);
					temp+=scratch[0]+scratch[1];
				}
				else
				{
					mValues=&matrixValues[xx][1];
					xPtr=&localXPtr[xx][(-Degree+j+1)/2];
					temp2=_mm_mul_pd(xPtr[-1],mValues[-1]);
					for(int i=0;i<Degree;i++)
					{
						temp1=_mm_mul_pd(xPtr[i],mValues[i]);
						temp2=_mm_add_pd(temp1,temp2);
					}
					_mm_store_pd(scratch,temp2);
					temp+=scratch[0]+scratch[1];
				}
			}
			return temp;
	}
}
template<class Real,int Type,int Degree,class Data>
inline double StreamingSolver<Real,Type,Degree,Data>::GaussSeidelUpdate0(const __m128d mValues[][TemplateCount],double diagonal,int j,int iB)
{
	return (localB[iB+j]-GetLaplacianValue0(mValues,j))/diagonal;
}
template<class Real,int Type,int Degree,class Data>
inline double StreamingSolver<Real,Type,Degree,Data>::GaussSeidelUpdate1(const __m128d mValues[][TemplateCount],double diagonal,int j,int iB)
{
	return (localB[iB+j]-GetLaplacianValue1(mValues,j))/diagonal;
}
