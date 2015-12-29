#ifndef STREAMING_SOLVER_128_INCLUDED
#define STREAMING_SOLVER_128_INCLUDED
#include <Util/XPlatform.h>

#define ScaledSum128( x1 , x2 , x3 , x4 ) ( _mm_add_ps( _mm_mul_ps( x1 , x2 ) , _mm_mul_ps( x3 , x4 ) ) )

// The dot-product functions
inline void SetInteriorDotSum( const __m128 matrixValues[] , ConstPointer( __m128 ) xPtrs[] , size_t j , __m128& dotSum )
{
	dotSum =                      _mm_mul_ps( xPtrs[0][j] , matrixValues[0] )  ;
	dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[1][j] , matrixValues[1] ) );
	dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[2][j] , matrixValues[2] ) );
}
inline void SetFullInteriorDotSum( const __m128 matrixValues[] , ConstPointer( __m128 ) xPtrs[] , size_t j , __m128& dotSum )
{
	dotSum =				      _mm_mul_ps( xPtrs[2][j] , matrixValues[2] ) ;
	dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[0][j] , matrixValues[0] ) );
	dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[1][j] , matrixValues[1] ) );
	dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[3][j] , matrixValues[3] ) );
	dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[4][j] , matrixValues[4] ) );
}

inline void SetDotSum( const __m128 matrixValues[] , ConstPointer( __m128 ) xPtrs[] , size_t j , __m128& dotSum )
{
	dotSum =									 _mm_mul_ps( xPtrs[2][j] , matrixValues[2] ) ;
	if( xPtrs[0] ) dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[0][j] , matrixValues[0] ) );
	if( xPtrs[1] ) dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[1][j] , matrixValues[1] ) );
	if( xPtrs[3] ) dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[3][j] , matrixValues[3] ) );
	if( xPtrs[4] ) dotSum = _mm_add_ps( dotSum , _mm_mul_ps( xPtrs[4][j] , matrixValues[4] ) );
}
inline void SetInteriorRestrictionDotSum( const __m128 matrixValues[] , ConstPointer( __m128 ) rPtrs[] , size_t j , __m128& dotSum )
{
	dotSum=						_mm_mul_ps(rPtrs[0][j],matrixValues[0]) ;
	dotSum=_mm_add_ps(dotSum,	_mm_mul_ps(rPtrs[1][j],matrixValues[1]));
}
inline void SetRestrictionDotSum( const __m128& matrixValues , ConstPointer( __m128 ) rPtrs , size_t j , __m128& dotSum )
{
	dotSum=_mm_mul_ps(rPtrs[j],matrixValues) ;
}
inline void SetRestrictionDotSum( const __m128 matrixValues[] , ConstPointer( __m128 ) rPtrs[] , size_t j , __m128& dotSum )
{
	dotSum=										_mm_mul_ps(rPtrs[2][j],matrixValues[2]) ;
	if(rPtrs[0])	dotSum=_mm_add_ps(dotSum,	_mm_mul_ps(rPtrs[0][j],matrixValues[0]));
	if(rPtrs[1])	dotSum=_mm_add_ps(dotSum,	_mm_mul_ps(rPtrs[1][j],matrixValues[1]));
	if(rPtrs[3])	dotSum=_mm_add_ps(dotSum,	_mm_mul_ps(rPtrs[3][j],matrixValues[3]));
}
/////////////////////////////////
// The base template functions //
/////////////////////////////////
template< typename DotFunction >
inline float GaussSeidelUpdate0( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F( mValues , xPtrs , j>>2 , dotSum );
	scratch[0] = scratch[0] + scratch[1] + scratch[2] + previousDotSum;
	previousDotSum = scratch[2] + scratch[3];
	return (localBPtr[j]-scratch[0]);
}
template< typename DotFunction >
inline float GaussSeidelUpdate1( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F(mValues,xPtrs,j>>2,dotSum);
	scratch[0]=(scratch[0]+scratch[1])+(scratch[2]+scratch[3])+previousDotSum;
	previousDotSum=scratch[3];
	return (localBPtr[j]-scratch[0]);
}
template< typename DotFunction >
inline float GaussSeidelUpdate2( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F(mValues,xPtrs,(j>>2)+1,dotSum);
	float temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return (localBPtr[j]-temp);
}
template< typename DotFunction >
inline float GaussSeidelUpdate3( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F(mValues,xPtrs,(j>>2)+1,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return (localBPtr[j]-scratch[0]);
}

template< typename DotFunction >
inline float GetLaplacianValue0( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j)
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F(mValues,xPtrs,j>>2,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[2]+scratch[3];
	return scratch[0];
}
template< typename DotFunction >
inline float GetLaplacianValue1( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F(mValues,xPtrs,j>>2,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+scratch[3]+previousDotSum;
	previousDotSum=scratch[3];
	return scratch[0];
}
template< typename DotFunction >
inline float GetLaplacianValue2( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	float temp;
	F(mValues,xPtrs,(j>>2)+1,dotSum);
	temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return temp;
}
template< typename DotFunction >
inline float GetLaplacianValue3( DotFunction F , const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F(mValues,xPtrs,(j>>2)+1,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return scratch[0];
}

template< typename DotFunction >
inline float ReverseGaussSeidelUpdate0( DotFunction F , const __m128 mValues[] , const __m128* xPtrs[] , const float* localBPtr , float& nextDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F( mValues , xPtrs , (j>>2)-1 , dotSum );
	scratch[3] = ( scratch[2] + scratch[3] ) + nextDotSum;
	nextDotSum = scratch[0] + scratch[1] + scratch[2];
	return (localBPtr[j]-scratch[3]);
}
template< typename DotFunction >
inline float ReverseGaussSeidelUpdate1( DotFunction F , const __m128 mValues[] , const __m128* xPtrs[] , const float* localBPtr , float& nextDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F( mValues , xPtrs , (j>>2)-1 , dotSum );
	float temp = nextDotSum + scratch[3];
	nextDotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
	return (localBPtr[j]-temp);
}
template< typename DotFunction >
inline float ReverseGaussSeidelUpdate2( DotFunction F , const __m128 mValues[] , const __m128* xPtrs[] , const float* localBPtr , float& nextDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F( mValues , xPtrs , (j>>2) , dotSum );
	scratch[3] = nextDotSum + scratch[0] + scratch[1] + scratch[2] + scratch[3];
	nextDotSum = scratch[0];
	return (localBPtr[j]-scratch[3]);
}
template< typename DotFunction >
inline float ReverseGaussSeidelUpdate3( DotFunction F , const __m128 mValues[] , const __m128* xPtrs[] , const float* localBPtr , float& nextDotSum , int j )
{
	__m128 dotSum;
	float* scratch = (float*)&dotSum;
	F( mValues , xPtrs , (j>>2) , dotSum );
	scratch[3] = nextDotSum + scratch[1] + scratch[2] + scratch[3];
	nextDotSum = scratch[0] + scratch[1];
	return (localBPtr[j]-scratch[3]);
}

inline float FullInteriorGaussSeidelUpdate0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetFullInteriorDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[2]+scratch[3];
	return (localBPtr[j]-scratch[0]);
}
inline float FullInteriorGaussSeidelUpdate1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetFullInteriorDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=(scratch[0]+scratch[1])+(scratch[2]+scratch[3])+previousDotSum;
	previousDotSum=scratch[3];
	return (localBPtr[j]-scratch[0]);
}
inline float FullInteriorGaussSeidelUpdate2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetFullInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	float temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return (localBPtr[j]-temp);
}
inline float FullInteriorGaussSeidelUpdate3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetFullInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return (localBPtr[j]-scratch[0]);
}

inline float GetFullInteriorLaplacianValue0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetFullInteriorDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[2]+scratch[3];
	return scratch[0];
}
inline float GetFullInteriorLaplacianValue1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetFullInteriorDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+scratch[3]+previousDotSum;
	previousDotSum=scratch[3];
	return scratch[0];
}
inline float GetFullInteriorLaplacianValue2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	float temp;
	SetFullInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return temp;
}
inline float GetFullInteriorLaplacianValue3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetFullInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return scratch[0];
}
inline float GaussSeidelUpdate0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum( mValues , xPtrs , j>>2 , dotSum );
	_mm_store_ps(scratch,dotSum);
	scratch[0] = ( scratch[0] + scratch[1] ) + ( scratch[2] + previousDotSum );
	previousDotSum = scratch[2] + scratch[3];
	return (localBPtr[j]-scratch[0]);
}
inline float GaussSeidelUpdate1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+scratch[3]+previousDotSum;
	previousDotSum=scratch[3];
	return (localBPtr[j]-scratch[0]);
}
inline float GaussSeidelUpdate2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	float temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return (localBPtr[j]-temp);
}
inline float GaussSeidelUpdate3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return (localBPtr[j]-scratch[0]);
}
inline float ReverseGaussSeidelUpdate0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum( mValues , xPtrs , (j>>2)-1 , dotSum );
	_mm_store_ps( scratch , dotSum );
	scratch[3] = ( scratch[2] + scratch[3] ) + nextDotSum;
	nextDotSum = scratch[0] + scratch[1] + scratch[2];
	return (localBPtr[j]-scratch[3]);
}
inline float ReverseGaussSeidelUpdate1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum( mValues , xPtrs , (j>>2)-1 , dotSum );
	_mm_store_ps( scratch , dotSum );
	float temp = nextDotSum + scratch[3];
	nextDotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
	return (localBPtr[j]-temp);
}
inline float ReverseGaussSeidelUpdate2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum( mValues , xPtrs , (j>>2) , dotSum );
	_mm_store_ps(scratch,dotSum);
	scratch[3] = nextDotSum + scratch[0] + scratch[1] + scratch[2] + scratch[3];
	nextDotSum = scratch[0];
	return (localBPtr[j]-scratch[3]);
}
inline float ReverseGaussSeidelUpdate3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum( mValues , xPtrs , (j>>2) , dotSum );
	_mm_store_ps( scratch , dotSum );
	scratch[3] = nextDotSum + scratch[1] + scratch[2] + scratch[3];
	nextDotSum = scratch[0] + scratch[1];
	return (localBPtr[j]-scratch[3]);
}
inline float InteriorGaussSeidelUpdate0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum( mValues , xPtrs , j>>2 , dotSum );
	_mm_store_ps( scratch , dotSum );
	scratch[0] = scratch[0] + scratch[1] + scratch[2] + previousDotSum;
	previousDotSum = scratch[2] + scratch[3];
	return ( localBPtr[j] - scratch[0] );
}
inline float InteriorGaussSeidelUpdate1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=(scratch[0]+scratch[1])+(scratch[2]+scratch[3])+previousDotSum;
	previousDotSum=scratch[3];
	return (localBPtr[j]-scratch[0]);
}
inline float InteriorGaussSeidelUpdate2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	float temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return (localBPtr[j]-temp);
}
inline float InteriorGaussSeidelUpdate3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return (localBPtr[j]-scratch[0]);
}
inline float ReverseInteriorGaussSeidelUpdate0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum( mValues , xPtrs , (j>>2)-1 , dotSum );
	_mm_store_ps( scratch , dotSum );
	scratch[3] = ( scratch[2] + scratch[3] ) + nextDotSum;
	nextDotSum = scratch[0] + scratch[1] + scratch[2];
	return (localBPtr[j]-scratch[3]);
}
inline float ReverseInteriorGaussSeidelUpdate1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum( mValues , xPtrs , (j>>2)-1 , dotSum );
	_mm_store_ps( scratch , dotSum );
	float temp = nextDotSum + scratch[3];
	nextDotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
	return (localBPtr[j]-temp);
}
inline float ReverseInteriorGaussSeidelUpdate2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum( mValues , xPtrs , (j>>2) , dotSum );
	_mm_store_ps(scratch,dotSum);
	scratch[3] = nextDotSum + scratch[0] + scratch[1] + scratch[2] + scratch[3];
	nextDotSum = scratch[0];
	return (localBPtr[j]-scratch[3]);
}
inline float ReverseInteriorGaussSeidelUpdate3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , ConstPointer( float ) localBPtr , float& nextDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum( mValues , xPtrs , (j>>2) , dotSum );
	_mm_store_ps( scratch , dotSum );
	scratch[3] = nextDotSum + scratch[1] + scratch[2] + scratch[3];
	nextDotSum = scratch[0] + scratch[1];
	return (localBPtr[j]-scratch[3]);
}
inline float GetLaplacianValue0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[2]+scratch[3];
	return scratch[0];
}
inline float GetLaplacianValue1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+scratch[3]+previousDotSum;
	previousDotSum=scratch[3];
	return scratch[0];
}
inline float GetLaplacianValue2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	float temp;
	SetDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return temp;
}
inline float GetLaplacianValue3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return scratch[0];
}
inline float GetInteriorLaplacianValue0( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[2]+scratch[3];
	return scratch[0];
}
inline float GetInteriorLaplacianValue1( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum(mValues,xPtrs,j>>2,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+scratch[3]+previousDotSum;
	previousDotSum=scratch[3];
	return scratch[0];
}
inline float GetInteriorLaplacianValue2( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	float temp;
	SetInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	temp=(previousDotSum+scratch[0]);
	previousDotSum=(scratch[0]+scratch[1])+(scratch[2]+scratch[3]);
	return temp;
}
inline float GetInteriorLaplacianValue3( const __m128 mValues[] , ConstPointer( __m128 ) xPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorDotSum(mValues,xPtrs,(j>>2)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+previousDotSum;
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return scratch[0];
}
inline float InteriorRestrictionUpdate0( const __m128 mValues[] , ConstPointer( __m128 ) rPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorRestrictionDotSum(mValues,rPtrs,j>>1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[3];
	return scratch[0];
}
inline float InteriorRestrictionUpdate1( const __m128 mValues[] , ConstPointer( __m128 ) rPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetInteriorRestrictionDotSum(mValues,rPtrs,(j>>1)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	float temp=(previousDotSum+scratch[0]);
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return temp;
}
inline float RestrictionUpdate0( const __m128& mValues , ConstPointer( __m128 ) rPtrs , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetRestrictionDotSum(mValues,rPtrs,j>>1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[3];
	return scratch[0];
}
inline float RestrictionUpdate1( const __m128& mValues , ConstPointer( __m128 ) rPtrs , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetRestrictionDotSum(mValues,rPtrs,(j>>1)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	float temp=(previousDotSum+scratch[0]);
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return temp;
}
inline float RestrictionUpdate0( const __m128 mValues[] , ConstPointer( __m128 ) rPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetRestrictionDotSum(mValues,rPtrs,j>>1,dotSum);
	_mm_store_ps(scratch,dotSum);
	scratch[0]=scratch[0]+scratch[1]+scratch[2]+previousDotSum;
	previousDotSum=scratch[3];
	return scratch[0];
}
inline float RestrictionUpdate1( const __m128 mValues[] , ConstPointer( __m128 ) rPtrs[] , float& previousDotSum , int j )
{
	ALIGN( float scratch[4] , 16 );
	__m128 dotSum;
	SetRestrictionDotSum(mValues,rPtrs,(j>>1)+1,dotSum);
	_mm_store_ps(scratch,dotSum);
	float temp=(previousDotSum+scratch[0]);
	previousDotSum=scratch[1]+scratch[2]+scratch[3];
	return temp;
}
#endif // STREAMING_SOLVER_128_INCLUDED
