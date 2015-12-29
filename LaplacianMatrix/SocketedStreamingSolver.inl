#define FIX_DERAMP 1

#include <stdint.h>
#define MyFree( buffer ) if( buffer ) free( buffer ) , buffer = NULL

void SetDownSampledStencil( int dim , int iters , DotProductStencil& outDot , DotProductStencil& outD2Dot , bool lump )
{
	DotProductStencil dot , d2Dot;
	FiniteElements1D< double , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( dim ,   dot , 0 , 0 );
	FiniteElements1D< double , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( dim , d2Dot , 1 , 1 , false );

	if( lump )
	{
		for( int i=0 ; i<=2*Degree ; i++ )
		{
			double sum = 0;
			for( int j=0 ; j<=2*Degree ; j++ ) sum += dot.caseTable[i].values[j] , dot.caseTable[i].values[j] = 0;
			dot.caseTable[i].values[Degree] = sum;
		}
	}

	SetDownSampledStencil( dot , d2Dot , dim , iters , outDot , outD2Dot );
}
void SetDownSampledStencil( const DotProductStencil& inDot , const DotProductStencil& inD2Dot , int dim , int iters , DotProductStencil& outDot , DotProductStencil& outD2Dot )
{
	if( !iters ) outDot = inDot , outD2Dot = inD2Dot;
	else
	{
		DotProductStencil newDot , newD2Dot;
		FiniteElements1D< double , Type , Degree >::FullProlongationStencil prolongationStencil;
		FiniteElements1D< double , Type , Degree >::ProlongationStencil( dim >>1 , prolongationStencil , dim );
		CombineStencils< double , Type , Degree >(    inDot , prolongationStencil , dim ,   newDot );
		CombineStencils< double , Type , Degree >(  inD2Dot , prolongationStencil , dim , newD2Dot );
		SetDownSampledStencil( newDot , newD2Dot , dim>>1 , iters-1 , outDot , outD2Dot );
	}
}
////////////////////////
// SocketedStreamData //
////////////////////////

bool SocketedStreamData::_set( int width , int height , int start , int end , DataStream* leftStream , Socket syncSocket , DataStream* rightStream , int periodicType )
{
	if(width&3)														{ fprintf(stderr,"width must be a multiple of four: %d\n",width) ; return false; }
	if(start&3)														{ fprintf(stderr,"start must be a multiple of four: %d\n",start) ; return false; }
	if(end&3)														{ fprintf(stderr,"end must be a multiple of four: %d\n",end)	 ; return false; }
	if( (start || periodicType!=NO_PERIODIC ) && !leftStream )			{ fprintf(stderr,"left socket invalid for stream\t%dx%d\n",width,height) ; return false; }
	if( (end!=width || periodicType!=NO_PERIODIC ) && !rightStream )	{ fprintf(stderr,"right socket invalid for stream\t%dx%d\n",width,height)  ; return false; }
	if( periodicType==SPHERICAL_PERIODIC && syncSocket==_INVALID_SOCKET_ ){ fprintf(stderr,"sync socket invalid for spherical stream\n") ; return false; }

	this->leftStream = leftStream;
	this->rightStream =rightStream;
	syncXSocket = syncRSocket = syncSocket;

	_start128=start>>2;
	_end128=end	>>2;
	_size128=_end128-_start128;

	_start=_start128<<2;
	_end=_end128<<2;
	_size=_size128<<2;

	_paddedSize128=_size128+2*_padSize128;
	_padSize=_padSize128<<2;
	_paddedSize=_paddedSize128<<2;

	return true;
}

// This guys gets called by the streaming solver

bool SocketedStreamData::SetSocketedStreamData(int width,int height,int start,int end,int iters,
											   DataStream* leftStream,Socket syncSocket,DataStream* rightStream,
											   int periodicType )
{
	// It's iters+2 because we need:
	//	1] iters-1 for offsetting the solver
	//	2] 1 for reading in the solver
	//	3] 1 for computing the residual without having to re-sync.
	//	4] 1 for zig-zag since you start one in.
	_padSize128=(iters+2)*WordPerDegree;

	if( !_set( width , height , start , end , leftStream , syncSocket , rightStream , periodicType ) )	return false;

	if(end-start<4*_padSize)	{ fprintf(stderr,"stream too narrow: %d < 4*%d\n",end-start,_padSize) ; return false; }
	if(height<2*Degree)			{ fprintf(stderr,"stream too short: %d < 2*%d\n",height,Degree)		  ; return false; }
	return true;
}
// This guys gets called by the streaming laplacian
bool SocketedStreamData::SetSocketedStreamData(int width,int height,int start,int end,
											   DataStream* leftStream,Socket syncSocket,DataStream* rightStream,
											   int periodicType )
{
	_padSize128=WordPerDegree;

	if( !_set( width , height , start , end , leftStream , syncSocket , rightStream , periodicType ) )	return false;

	if(end-start<2*_padSize)	{ fprintf(stderr,"stream too narrow: %d < 2*%d\n",end-start,_padSize) ; return false; }
	return true;
}
int SocketedStreamData::start	(void)	const	{ return _start; }
int SocketedStreamData::end		(void)	const	{ return _end; }
int SocketedStreamData::size	(void)	const	{ return _size; }

/////////////////////////////
// SocketedStreamingSolver //
/////////////////////////////
template< int Channels , class SyncType > int SocketedStreamingSolver< Channels , SyncType >::OffsetX(int iters)	{	return (iters+1)*Degree;	}
template< int Channels , class SyncType > int SocketedStreamingSolver< Channels , SyncType >::OffsetB(int iters)	{	return iters*Degree;	}
//template<int Channels>	int SocketedStreamingSolver< Channels , SyncType >::OffsetR(void)	{	return -1;	}
template< int Channels , class SyncType > int SocketedStreamingSolver< Channels , SyncType >::OffsetR(void)		{	return -3;	}

template< int Channels , class SyncType >
SocketedStreamingSolver< Channels , SyncType >::SocketedStreamingSolver(void)
{
	_server = NULL;
	_deleteServer = true;
	showProgress=false;
	progressCount=0;

#if CLEAN_RESIDUAL
	setResidual = false;
#else // !CLEAN_RESIDUAL
	setResidual = true;
#endif // CLEAN_RESIDUAL
	xSize = bSize = rSize = 0;
	laneNum = 1;

	RStream = NullPointer< Pointer( __m128 ) >( );
	XStream = NullPointer< Pointer( __m128 ) >( );
	BStream = NullPointer< Pointer( __m128 ) >( );
	for( int d=0 ; d<Degree*Channels ; d++ ) RBuffer[d] = XBuffer[d] = NullPointer< __m128 >( );

	lapTemplates = AlignedAllocPointer< TemplateSSE >( 3 * (2*Degree+1) , ALIGNMENT );
	zeroLapTemplates = AlignedAllocPointer< TemplateSSE >( 3 * (2*Degree+1)  , ALIGNMENT );

	syncBuffer = NullPointer< SyncType >( );
	for( int i=0 ; i<Degree ; i++ ) localXAccum[i] = NullPointer<  __m128 >( );

}
template< int Channels , class SyncType >
SocketedStreamingSolver< Channels , SyncType >::~SocketedStreamingSolver(void)
{
	freeStreams();
	FreePointer( syncBuffer );
	AlignedFreePointer( lapTemplates );
	AlignedFreePointer( zeroLapTemplates );
	for( int i=0 ; i<Degree ; i++ )
		if( localXAccum[i] )
		{
			localXAccum[i] -= _padSize128;
			FreePointer( localXAccum[i] );
		}
	if( _server && _deleteServer ) delete _server , _server=NULL;
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::freeStreams( void )
{
	if( XStream )
	{
		for( int i=0 ; i<xSize*Channels ; i++ )
			if( XStream[i] )
			{
				XStream[i] -= _padSize128;
				AlignedFreePointer( XStream[i] );
			}
		FreePointer( XStream );
	}
	if( BStream )
	{
		for( int i=0 ; i<bSize*Channels ; i++ )
			if( BStream[i] )
			{
				BStream[i] -= _padSize128;
				AlignedFreePointer( BStream[i] );
			}
		FreePointer( BStream );
	}
	if( RStream )
	{
		for( int i=0 ; i<rSize*Channels ; i++ )
			if( RStream[i] )
			{
				RStream[i] -= _padSize128;
				AlignedFreePointer( RStream[i] );
			}
		FreePointer( RStream );
	}
	if( RBuffer ) for( int i=0 ; i<Degree*Channels ; i++ ) AlignedFreePointer( RBuffer[i] );
	if( XBuffer ) for( int i=0 ; i<Degree*Channels ; i++ ) AlignedFreePointer( XBuffer[i] );

	RStream = NullPointer< Pointer( __m128 ) >( );
	XStream = NullPointer< Pointer( __m128 ) >( );
	BStream = NullPointer< Pointer( __m128 ) >( );

	rSize=xSize=bSize=0;
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Init( int start , int end , int major , int minor , int iters , 
											 DataStream* leftStream , Socket syncSocket , DataStream* rightStream ,
											 int periodicType , MultiStreamIOServer* server
											 )
{
	MatrixStencil lStencil;
	FiniteElements2D< double , Type , Degree >::LaplacianStencil( major , minor , lStencil , true );
	Init( lStencil , start , end , major , minor , iters ,
		leftStream , syncSocket , rightStream ,
		periodicType , server
		);
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Init( const MatrixStencil& lStencil , int start , int end , int major , int minor , int iters ,
											 DataStream* leftStream , Socket syncSocket , DataStream* rightStream ,
											 int periodicType , MultiStreamIOServer* server
											 )
{

	if( !SetSocketedStreamData( major , minor , start , end , iters , leftStream , syncSocket , rightStream , periodicType ) )	exit(0);
	if( _server && _deleteServer ) delete _server , _server = NULL;
	if( server ) _server = server , _deleteServer = false;
	else		 _server = new MultiStreamIOServer() , _deleteServer = true;

	this->periodicType = periodicType;
	this->major=major;
	this->minor=minor;
	this->iters=iters;
	if( periodicType==SPHERICAL_PERIODIC )
	{
		FreePointer( syncBuffer );
		syncBuffer = AllocPointer< SyncType >( _size * Channels * Degree );
	}

	MatrixStencil laplacianStencil;
	laplacianStencil = lStencil;
	laplacianScale   = float( 1./laplacianStencil.caseTable[Degree][Degree].values[Degree][Degree] );
	laplacianScaleR  = float(    laplacianStencil.caseTable[Degree][Degree].values[Degree][Degree] );

	for( int i=0 ; i<=2*Degree ; i++ ) for( int j=0 ; j<=2*Degree ; j++ )
		for( int k=0 ; k<=2*Degree ; k++ ) for( int l=0 ; l<=2*Degree ; l++ )
			laplacianStencil.caseTable[i][j].values[k][l] *= laplacianScale;
	// BADNESS!!! should fix the indexing of lapTemplates so that major index iterates faster
	ALIGN( float scratch[4] , 16 );
	for( int i=0 ; i<=2*Degree ; i++ )	// Iterate over the minor index in the mask
	{
		for( int k=0 ; k<4 ; k++ ) lapTemplates    [3*i+1].diagonalR[k] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		for( int k=0 ; k<4 ; k++ ) zeroLapTemplates[3*i+1].diagonalR[k] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		lapTemplates    [3*i  ].diagonalR[0] = float( 1./laplacianStencil.caseTable[i][0].values[Degree][Degree] );
		lapTemplates    [3*i  ].diagonalR[1] = float( 1./laplacianStencil.caseTable[i][1].values[Degree][Degree] );
		lapTemplates    [3*i  ].diagonalR[2] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		lapTemplates    [3*i  ].diagonalR[3] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		lapTemplates    [3*i+2].diagonalR[0] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		lapTemplates    [3*i+2].diagonalR[1] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		lapTemplates    [3*i+2].diagonalR[2] = float( 1./laplacianStencil.caseTable[i][2*Degree-1].values[Degree][Degree] );
		lapTemplates    [3*i+2].diagonalR[3] = float( 1./laplacianStencil.caseTable[i][2*Degree  ].values[Degree][Degree] );
		zeroLapTemplates[3*i  ].diagonalR[0] = float( 1./laplacianStencil.caseTable[i][0].values[Degree][Degree] );
		zeroLapTemplates[3*i  ].diagonalR[1] = float( 1./laplacianStencil.caseTable[i][1].values[Degree][Degree] );
		zeroLapTemplates[3*i  ].diagonalR[2] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		zeroLapTemplates[3*i  ].diagonalR[3] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		zeroLapTemplates[3*i+2].diagonalR[0] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		zeroLapTemplates[3*i+2].diagonalR[1] = float( 1./laplacianStencil.caseTable[i][Degree].values[Degree][Degree] );
		zeroLapTemplates[3*i+2].diagonalR[2] = float( 1./laplacianStencil.caseTable[i][2*Degree-1].values[Degree][Degree] );
		zeroLapTemplates[3*i+2].diagonalR[3] = float( 1./laplacianStencil.caseTable[i][2*Degree  ].values[Degree][Degree] );
	}
	for(int i=0;i<=2*Degree;i++)	// Iterate over the minor index in the mask
		for(int j=0;j<3;j++)
			for(int k=0;k<=2*Degree;k++)
			{
				int jj;

				if(j==0)	jj=0;
				else		jj=Degree;
				scratch[0] = float( laplacianStencil.caseTable[i][jj].values[k][2] );
				scratch[1] = float( laplacianStencil.caseTable[i][jj].values[k][3] );
				scratch[2] = float( laplacianStencil.caseTable[i][jj].values[k][4] );
				if(j!=2)	scratch[3] = float( laplacianStencil.caseTable[i][Degree].values[k][1] );
				else		scratch[3] = float( 0 );
				lapTemplates[3*i+j].matrixValues[0][k]=_mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l] = float( laplacianStencil.caseTable[i][jj].values[k][l+1] );
				lapTemplates[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l] = float( laplacianStencil.caseTable[i][jj].values[k][l] );
				lapTemplates[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0] = float( laplacianStencil.caseTable[i][Degree].values[k][3] );
				else		scratch[0] = float( 0 );
				scratch[1] = float( laplacianStencil.caseTable[i][jj].values[k][0] );
				scratch[2] = float( laplacianStencil.caseTable[i][jj].values[k][1] );
				scratch[3] = float( laplacianStencil.caseTable[i][jj].values[k][2] );
				lapTemplates[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
			}
	// Zero out the diagonal entry
	for( int i=0 ; i<=2*Degree ; i++ ) for( int j=0 ; j<=2*Degree ; j++ ) laplacianStencil.caseTable[i][j].values[Degree][Degree] = 0;
	for( int i=0 ; i<=2*Degree ; i++ )	// Iterate over the minor index in the mask
		for(int j=0;j<3;j++)
			for(int k=0;k<=2*Degree;k++)
			{
				int jj;

				if(j==0)	jj=0;
				else		jj=Degree;
				scratch[0] = float( laplacianStencil.caseTable[i][jj].values[k][2] );
				scratch[1] = float( laplacianStencil.caseTable[i][jj].values[k][3] );
				scratch[2] = float( laplacianStencil.caseTable[i][jj].values[k][4] );
				if(j!=2)	scratch[3] = float( laplacianStencil.caseTable[i][Degree].values[k][1] );
				else		scratch[3] = float( 0 );
				zeroLapTemplates[3*i+j].matrixValues[0][k]=_mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l] = float( laplacianStencil.caseTable[i][jj].values[k][l+1] );
				zeroLapTemplates[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l] = float( laplacianStencil.caseTable[i][jj].values[k][l] );
				zeroLapTemplates[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0] = float( laplacianStencil.caseTable[i][Degree].values[k][3] );
				else		scratch[0] = float( 0 );
				scratch[1] = float( laplacianStencil.caseTable[i][jj].values[k][0] );
				scratch[2] = float( laplacianStencil.caseTable[i][jj].values[k][1] );
				scratch[3] = float( laplacianStencil.caseTable[i][jj].values[k][2] );
				zeroLapTemplates[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
			}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Set( int rSize ){ Set( -OffsetX , 0 , 0 , 0 , 0 , rSize ); }
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Set( int start , int bStart , int xStart , int bEnd , int xEnd , int rSize )
{
	// At index idx, we solve for the X coefficients
	// [idx+Degree-1,idx+2*Degree-1,...,idx+iters*Degree-1]
	// Which requires read access to the X coefficients:
	// [idx-1,idx+(iters+1)*Degree-1]
	// and the B coefficients
	// [idx+Degree-1,idx+iters*Degree-1]

	// After solving at idx, the values of X at coefficient idx+Degree-1 (and below) are finalized.
	// Which means that we can solve for the value of the residual in row idx-1, requiring read access to the X coefficients:
	// [idx-Degree-1,idx+Degree-1]
	// and the B coefficients
	// [idx-1]

	clearX = clearB = true;
	index = start;
	freeStreams();

	xSize += ( OffsetX(iters) > xEnd ? OffsetX(iters) : xEnd ) + 1;
	bSize += ( OffsetB(iters) > bEnd ? OffsetB(iters) : bEnd ) + 1;
	if( setResidual )
	{
		xSize+=-Degree+OffsetR()<xStart?Degree-OffsetR():-xStart;
		bSize+=        OffsetR()<bStart?      -OffsetR():-bStart;
	}
	else
	{
		xSize += -Degree<xStart ? Degree : - xStart;
		bSize +=       0<bStart ?      0 : - bStart;
	}
#if SAME_SIZE_BUFFERS
	bSize = xSize;
#endif // SAME_SIZE_BUFFERS
	if( bSize>minor )
	{
		bSize=minor;
		clearB=false;
	}
	if( xSize>minor )
	{
		xSize=minor;
		clearX=false;
	}

	this->rSize=rSize;
	if(setResidual)
	{
		RStream = AllocPointer< Pointer( __m128 ) >( rSize * Channels );
		for( int i=0 ; i<rSize*Channels ; i++ )
		{
			RStream[ i ]  = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
			RStream[ i ] += _padSize128;
		}
	}
	XStream = AllocPointer< Pointer( __m128 ) >( xSize * Channels );
	for( int i=0 ; i<xSize*Channels ; i++ )
	{
		XStream[ i ]  = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
		memset( XStream[i] , 0 , sizeof( __m128 ) * _paddedSize128 );
		XStream[ i ] += _padSize128;
	}

	BStream = AllocPointer< Pointer( __m128 ) >( bSize * Channels );
	for( int i=0 ; i<bSize*Channels ; i++ )
	{
		BStream[ i ] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
		BStream[ i ] += _padSize128;
	}

	for( int i=0 ; i<Degree ; i++ )
	{
		localXAccum[i] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
		localXAccum[i] += _padSize128;
	}

	if( periodicType==SPHERICAL_PERIODIC )
	{
		for( int d=0 ; d<Degree*Channels ; d++ ) XBuffer[d] = AlignedAllocPointer< __m128 >(  _paddedSize128 , ALIGNMENT );
		if( setResidual ) for( int d=0 ; d<Degree*Channels ; d++ ) RBuffer[d] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
	}
}
template< int Channels , class SyncType >
Pointer( float ) SocketedStreamingSolver< Channels , SyncType >::GetXRow( int row , int channel )
{
	if( row<0 )
		if( periodicType==SPHERICAL_PERIODIC && row >=-Degree)		return ( Pointer( float ) )( XBuffer[(Degree+row)*Channels+channel] + _padSize128 );
		else								return NullPointer< float >( );
	else if(row>=minor)
		if( periodicType==SPHERICAL_PERIODIC && row < minor+Degree)	return ( Pointer( float ) )( XBuffer[(row-minor)*Channels+channel] + _padSize128 );
		else								return NullPointer< float >( );
	else									return ( Pointer( float ) )XStream[ MyModIndex( row*Channels+channel , xSize*Channels ) ];
}
template< int Channels , class SyncType >
Pointer( float ) SocketedStreamingSolver< Channels , SyncType >::GetBRow(int row,int channel)
{
	return ( Pointer( float ) )BStream[MyModIndex(row*Channels+channel , bSize*Channels)];
}
template< int Channels , class SyncType >
Pointer( float ) SocketedStreamingSolver< Channels , SyncType >::GetRRow(int row,int channel)
{
	if(row<0)
		if( periodicType==SPHERICAL_PERIODIC && row >=-Degree)		return ( Pointer( float ) )( RBuffer[(Degree+row)*Channels + channel] + _padSize128 );
		else								return NullPointer< float >( );
	else if(row>=minor)
		if( periodicType==SPHERICAL_PERIODIC && row < minor+Degree)	return ( Pointer( float ) )( RBuffer[(row-minor)*Channels + channel] + _padSize128 );
		else								return NullPointer< float >( );
	else									return ( Pointer( float ) )RStream[MyModIndex(row*Channels+channel,rSize*Channels)];
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::UnSet(void)
{
	clearX = clearB = true;
	index = 0;
	setResidual = false;

	freeStreams();
	for( int i=0 ; i<Degree ; i++ )
		if( localXAccum[i] )
		{
			localXAccum[i] -= _padSize128;
			AlignedFreePointer( localXAccum[i] );
		}
}
template< int Channels , class SyncType >
bool SocketedStreamingSolver< Channels , SyncType >::Increment(void)
{
	int idx;
	if( clearX )
	{
		if( setResidual ) idx=index-Degree+OffsetR()+xSize;
		else              idx=index-Degree          +xSize;
		if( idx>=0 && idx<minor ) for( int c=0 ; c<Channels ; c++ ) memset( GetXRow( idx , c ) , 0 , sizeof( __m128 ) * _size128 );
	}
	if( clearB )
	{
#if SAME_SIZE_BUFFERS
		if( setResidual ) idx=index-Degree+OffsetR()+bSize;
		else              idx=index-Degree          +bSize;
#else // !SAME_SIZE_BUFFERS
		if( setResidual ) idx=index+OffsetR()+bSize;
		else              idx=index          +bSize;
#endif // SAME_SIZE_BUFFERS
		if( idx>=0 && idx<minor ) for( int c=0 ; c<Channels ; c++ ) memset( GetBRow( idx , c ) , 0 , sizeof( __m128 ) * _size128 );
	}
	index++;
	int restrictionIndex  = FiniteElements1D< float , Type , Degree >::FullRestrictionStencil::RestrictionStencil::Start( index+OffsetR() );
	int prolongationIndex = FiniteElements1D< float , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( index-1 );
	if( restrictionIndex>=minor/2 && index>=minor && prolongationIndex>=minor*2 && (!setResidual || (index+OffsetR()>=minor) ) ) return false;
	else return true;
}
template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateXInput( StreamingGrid* X )
{
	if( !X ) return false;
	int idx = index+OffsetX(iters);
	// Read in from the X vector
	if( idx>=0 && idx<minor )
	{
		Pointer( StorageType ) xPtr = ( Pointer( StorageType ) )(*X)[idx];
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) xRow = GetXRow( idx , c );
			for( int jj=0 ; jj<_size ; jj++ ) xRow[jj] += float( xPtr[jj*Channels+c] );
		}
		return true;
	}
	return false;
}
template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateBInput( StreamingGrid* B )
{
	if( !B ) return false;
	int idx = index+OffsetB(iters);
	// Read in from the B vector
	if( idx>=0 && idx<minor )
	{
		Pointer( StorageType ) bPtr = ( Pointer( StorageType ) )(*B)[idx];
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) bRow = GetBRow( idx , c );
			for( int jj=0 ; jj<_size ; jj++ ) bRow[jj] = float( bPtr[jj*Channels+c] );
		}
		return true;
	}
	return false;
}

template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateXOutput( StreamingGrid* X )
{
	if( !X ) return false;
	// Copy the solution
	if( index>=0 && index<minor )
	{
		Pointer( StorageType ) xPtr = ( Pointer( StorageType ) )(*X)[index];
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) xRow = GetXRow( index , c );
			for( int jj=0 ; jj<_size ; jj++ ) xPtr[jj*Channels+c] = StorageType( xRow[jj] );
		}
		return true;
	}
	return false;
}
template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateBOutput( StreamingGrid* B )
{
	if( !B ) return false;
	if( index>=0 && index<minor )
	{
		Pointer( StorageType ) bPtr = ( Pointer( StorageType ) )(*B)[index];
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) bRow = GetBRow( index , c );
			for( int jj=0 ; jj<_size ; jj++ ) bPtr[jj*Channels+c] = StorageType( bRow[jj] );
		}
		return true;
	}
	return false;
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverHead( int idx , bool read )
{
	if( idx<0 || idx>=Degree )
	{
		fprintf( stderr , "Synchronizing solver in non-head row: %d\n" , idx );
		return;
	}
	int r = -idx-1;
	if( syncXSocket!=_INVALID_SOCKET_ )
		if( read )
		{
			ReceiveOnSocket( syncXSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverHead" );
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( r , c );
				Pointer( SyncType ) sRow = syncBuffer + _size*c;
				for( int i=0 ; i<_size ; i++ ) xRow[i] = float( sRow[i] );
			}
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( idx , c );
				Pointer( SyncType ) sRow = syncBuffer+_size*c;
				for( int i=0 ; i<_size ; i++ ) sRow[i] = SyncType( xRow[i] );
			}
			SendOnSocket( syncXSocket , ( ConstPointer( SyncType ) )syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverHead" );
		}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverTail( int idx , bool read )
{
	if( idx<minor-Degree || idx>=minor )
	{
		fprintf( stderr , "Synchronizing solver in non-tail row: %d\n" , idx );
		return;
	}
	int r = 2 * minor - 1 - idx;
	if( syncXSocket!=_INVALID_SOCKET_ )
		if( read )
		{
			ReceiveOnSocket( syncXSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverTail" );
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( r , c );
				Pointer( SyncType ) sRow = syncBuffer+_size*c;
				for( int i=0 ; i<_size ; i++ ) xRow[i] = float( sRow[i] );
			}
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( idx , c );
				Pointer( SyncType ) sRow = syncBuffer+_size*c;
				for( int i=0 ; i<_size ; i++ ) sRow[i] = SyncType( xRow[i] );
			}
			SendOnSocket( syncXSocket , ( ConstPointer( SyncType ) )syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverTail" );
		}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncResidualHead( int idx , bool read , bool syncPeriodic )
{
	if( idx<0 || idx>=Degree )
	{
		fprintf( stderr , "Synchronizing residual in non-head row: %d\n" , idx );
		return;
	}
	int r = -idx-1;
	if( syncRSocket!=_INVALID_SOCKET_ )
	{
		if( !syncPeriodic )
			if( read )
			{
				ReceiveOnSocket( syncRSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualHead" );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( idx , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				SendOnSocket( syncRSocket , ( ConstPointer( SyncType ) )syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualHead" );
			}
		else
		{
			int size = WordPerDegree*RealPerWord;
			if( read )
			{
				leftStream->read( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
				rightStream->read( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r  , c ) + _size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				leftStream->write( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) + _size - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				rightStream->write( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
			}
		}
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncResidualTail( int idx , bool read , bool syncPeriodic )
{
	if( idx<minor-Degree || idx>=minor )
	{
		fprintf( stderr , "Synchronizing residual in non-tail row: %d\n" , idx );
		return;
	}
	int r = 2 * minor - 1 - idx;
	if( syncRSocket!=_INVALID_SOCKET_ )
	{
		if( !syncPeriodic )
			if( read )
			{
				ReceiveOnSocket( syncRSocket , syncBuffer , sizeof( SyncType ) * _size * Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualTail" );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( idx , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				SendOnSocket( syncRSocket , ( ConstPointer( SyncType ) )syncBuffer , sizeof(SyncType)*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualTail" );
			}
		else
		{
			int size = WordPerDegree*RealPerWord;
			if(read)
			{
				leftStream->read( ( Pointer( byte ) )syncBuffer , sizeof(SyncType) * size * Channels );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
				rightStream->read( ( Pointer( byte ) )syncBuffer , sizeof(SyncType)*size*Channels );
				for( int c=0 ; c<Channels ; c++)
				{
					Pointer( float ) rRow = GetRRow( r , c ) + _size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				leftStream->write( ( Pointer( byte ) )syncBuffer , sizeof(SyncType)*size*Channels );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) + _size - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				rightStream->write( ( Pointer( byte ) )syncBuffer , sizeof(SyncType)*size*Channels );
			}
		}
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverLeft( int j , bool read , bool overlapped )
{
	if( leftStream )
	{
		int xSize = _padSize;
		int offset = _padSize-2*WordPerDegree*RealPerWord;
		int ySize = (iters+1)*Degree+2;
		if(!iters)	ySize+=Degree;
		int size = (xSize+offset) * ySize * Channels + 2 * xSize * Channels;
		int yEnd = Degree+1;
		Pointer( SyncType ) left = AllocPointer< SyncType >( size );
		if( read )
		{
			if ( !leftStream->read( ( Pointer( byte ) )left , sizeof( SyncType )*size ) )	exit(0);
			for( int y=-(ySize-yEnd-1) ; y<=yEnd ; y++ )
			{
				int yy = y + (ySize-yEnd-1);
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ -xSize + x ] = float( left[ x * ySize * Channels + yy * Channels + c ] );
					}
			}
			if( overlapped )
				for( int k=0 ; k<iters ; k++ )
				{
					int y = -k*Degree;
					int yy = y + (ySize-yEnd-1);
					if( j+y>=0 && j+y<minor )
						for(int c=0;c<Channels;c++)
						{
							Pointer( float ) tmp = GetXRow( j+y , c );
							if( tmp ) for( int x=0 ; x<offset-WordPerDegree*RealPerWord*k ; x++ ) tmp[ x ] = float( left[ (xSize+x) * ySize * Channels + yy * Channels + c ] );
						}
				}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ -xSize+x ] = float( left[ (xSize+offset) * ySize * Channels + c * xSize + x ] );
			}
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ -xSize+x ] = float( left[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] );
				}
		}
		else
		{
			memset( left , 0 , sizeof( SyncType )*size );
			for( int y=-(ySize-yEnd-1) ; y<=yEnd ; y++ )
			{
				int yy = y + (ySize-yEnd-1);
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize+offset ; x++ ) left[ x * ySize * Channels + yy * Channels + c ] = SyncType( tmp[ x-offset ] );
					}
			}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) left[ (xSize+offset) * ySize * Channels + c * xSize + x ] = SyncType( tmp[ x ] );
			}
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) left[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] = SyncType( tmp[ x ] );
				}
			if ( !leftStream->write( ( Pointer( byte ) )left , sizeof( SyncType )*size ) )	exit(0);
		}
		FreePointer( left );
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverRight( int j , bool read , bool overlapped )
{
	if( rightStream )
	{
		int xSize = _padSize;
		int ySize = (iters+1)*Degree+2;
		if(!iters)	ySize+=Degree;
		int offset = _padSize-2*WordPerDegree*RealPerWord;
		int size = (xSize+offset) * ySize * Channels + 2 * xSize * Channels;
		int yEnd = Degree+1;
		Pointer( SyncType ) right = AllocPointer< SyncType >( size );
		if( read )
		{
			if ( !rightStream->read( ( Pointer( byte ) )right , sizeof( SyncType )*size ) )	exit(0);
			for( int y=-(ySize-yEnd-1) ; y<=yEnd ; y++ )
			{
				int yy = y + (ySize-yEnd-1);
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ _size + x ] = float( right[ (x+offset) * ySize * Channels + yy * Channels + c ] );
					}
			}
			if( overlapped )
				for( int k=0 ; k<iters ; k++ )
				{
					int y = -k*Degree;
					int yy = y + (ySize-yEnd-1);
					if( j+y>=0 && j+y<minor )
						for( int c=0 ; c<Channels ; c++ )
						{
							Pointer( float ) tmp = GetXRow( j+y , c );
							if( tmp ) for( int x=WordPerDegree*RealPerWord*k ; x<offset ; x++ ) tmp[ _size - offset + x ] = float( right[ x * ySize * Channels + yy * Channels + c ] );
						}
				}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ _size+x ] = float( right[ (xSize+offset) * ySize * Channels + c * xSize + x ] );
			}
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ _size+x ] = float( right[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] );
				}
		}
		else
		{
			memset(right,0, sizeof( SyncType )*size);
			for(int y=-(ySize-yEnd-1);y<=yEnd;y++)
			{
				int yy = y + (ySize-yEnd-1);
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize+offset ; x++ ) right[ x * ySize * Channels + yy * Channels + c ] = SyncType( tmp[ _size-xSize+x ] );
					}
			}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) right[ (xSize+offset) * ySize * Channels + c * xSize + x ] = SyncType( tmp[ _size - xSize + x ] );
			}
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) right[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] = SyncType( tmp[ _size - xSize + x ] );
				}
			if ( !rightStream->write( ( Pointer( byte ) )right , sizeof( SyncType )*size ) ) exit(0);
		}
		FreePointer( right );
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SolveInterior( int j , int c , int sSolve , int eSolve )
{
	int jj=Degree*3;
	{
		Pointer( float ) localXPtrs[2*Degree+1];
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=0 ; y<2*Degree+1 ; y++ ) localXPtrs[y] = GetXRow( j-Degree+y , c );
		{
			Pointer( __m128 ) xPtrs[] = { ( Pointer( __m128 ) )localXPtrs[0] , ( Pointer( __m128 ) )localXPtrs[1] , ( Pointer( __m128 ) )localXPtrs[3] , ( Pointer( __m128 ) )localXPtrs[4] };

			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for( int i=-WordPerDegree+((sSolve-_start)>>2) ; i<((eSolve-_start)>>2)+WordPerDegree ; i++ )
			{
				localXAccum[0][i] = _mm_add_ps( xPtrs[0][i] , xPtrs[3][i] );
				localXAccum[1][i] = _mm_add_ps( xPtrs[1][i] , xPtrs[2][i] );
			}
		}

		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) xPtrs[] = { ( Pointer( __m128 ) )localXAccum[0] , ( Pointer( __m128 ) )localXAccum[1] , ( Pointer( __m128 ) )localXPtrs[2] };
		int s,e;

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localXPtrs[2][0] = InteriorGaussSeidelUpdate0( zeroLapTemplates[jj].matrixValues[0] , xPtrs , localBPtr , dotSum , 0 ) * zeroLapTemplates[jj].diagonalR[0];
				s=4;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
				s=sSolve-_start;
			}
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localXPtrs[2][i]=InteriorGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i );
			int i=eSolve-_start-4;
			if( eSolve==major && periodicType==NO_PERIODIC ) localXPtrs[2][i] = InteriorGaussSeidelUpdate0(zeroLapTemplates[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonalR[0];
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localXPtrs[2][1]=InteriorGaussSeidelUpdate1(zeroLapTemplates[jj].matrixValues[1],xPtrs,localBPtr,dotSum,1)*zeroLapTemplates[jj].diagonalR[1];
				s=5;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
				s=sSolve-_start+1;
			}
			if( eSolve==major && periodicType==NO_PERIODIC)	e=eSolve-_start-4;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)				localXPtrs[2][i]=InteriorGaussSeidelUpdate1(mValues,xPtrs,localBPtr,dotSum,i);
			int i=eSolve-_start-3;
			if( eSolve==major && periodicType==NO_PERIODIC )		localXPtrs[2][i]=InteriorGaussSeidelUpdate1(zeroLapTemplates[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonalR[1];
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(zeroLapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localXPtrs[2][i]=InteriorGaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,i);
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-6;
				localXPtrs[2][i]=InteriorGaussSeidelUpdate2(zeroLapTemplates[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,i);
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonalR[2];
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(zeroLapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localXPtrs[2][i]=InteriorGaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,i);
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-5;
				localXPtrs[2][i]=InteriorGaussSeidelUpdate3(zeroLapTemplates[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,i);
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonalR[3];
			}
		}
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Solve( int j , int c , int sSolve , int eSolve )
{
	int jj;
	if      ( j<Degree        ) jj=j;
	else if	( j>minor-1-Degree) jj=2*Degree+(j-(minor-1));
	else                        jj=Degree;
	if( periodicType==NO_PERIODIC )
	{
		if( sSolve<0     ) sSolve = 0;
		if( eSolve>major ) eSolve = major;
	}
	if( jj==Degree || periodicType==SPHERICAL_PERIODIC ) return SolveInterior( j , c , sSolve , eSolve );
	jj*=3;

	{
		Pointer( float ) localXPtrs[2*Degree+1];
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=-Degree ; y<=Degree ; y++ )
			if(j+y>=0 && j+y<minor)	localXPtrs[y+Degree] = GetXRow( j+y , c );
			else                    localXPtrs[y+Degree] = NullPointer< float >( );
		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) xPtrs[] =
		{
			( ConstPointer( __m128 ) )localXPtrs[0],
			( ConstPointer( __m128 ) )localXPtrs[1],
			( ConstPointer( __m128 ) )localXPtrs[2],
			( ConstPointer( __m128 ) )localXPtrs[3],
			( ConstPointer( __m128 ) )localXPtrs[4]
		};
		int s,e;

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2],
				zeroLapTemplates[jj+1].matrixValues[0][3],
				zeroLapTemplates[jj+1].matrixValues[0][4]
			};
			float d=zeroLapTemplates[jj+1].diagonalR[0];
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localXPtrs[2][0] = GaussSeidelUpdate0( zeroLapTemplates[jj].matrixValues[0] , xPtrs , localBPtr , dotSum , 0 ) * zeroLapTemplates[jj].diagonalR[0];
				s=4;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
				s=sSolve-_start;
			}
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
			else							e=eSolve-_start;
			for( int i=s ; i<e ; i+=4 ) localXPtrs[2][i]=GaussSeidelUpdate0(mValues,xPtrs,localBPtr,dotSum,i)*d;
			int i=eSolve-_start-4;
			if( eSolve==major && periodicType==NO_PERIODIC )		localXPtrs[2][i]=GaussSeidelUpdate0(zeroLapTemplates[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonalR[0];
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2],
				zeroLapTemplates[jj+1].matrixValues[1][3],
				zeroLapTemplates[jj+1].matrixValues[1][4]
			};
			float d=zeroLapTemplates[jj+1].diagonalR[1];
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localXPtrs[2][1]=GaussSeidelUpdate1(zeroLapTemplates[jj].matrixValues[1],xPtrs,localBPtr,dotSum,1)*zeroLapTemplates[jj].diagonalR[1];
				s=5;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
				s=sSolve-_start+1;
			}
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
			else							e=eSolve-_start;
			for( int i=s ; i<e ; i+=4 ) localXPtrs[2][i]=GaussSeidelUpdate1(mValues,xPtrs,localBPtr,dotSum,i)*d;
			int i=eSolve-_start-3;
			if( eSolve==major && periodicType==NO_PERIODIC )		localXPtrs[2][i]=GaussSeidelUpdate1(zeroLapTemplates[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonalR[1];
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2],
				zeroLapTemplates[jj+1].matrixValues[2][3],
				zeroLapTemplates[jj+1].matrixValues[2][4]
			};
			float d=zeroLapTemplates[jj+1].diagonalR[2];
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(zeroLapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for( int i=s ; i<e ; i+=4 ) localXPtrs[2][i]=GaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,i)*d;
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-6;
				localXPtrs[2][i]=GaussSeidelUpdate2(zeroLapTemplates[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,i)*d;
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonalR[2];
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2],
				zeroLapTemplates[jj+1].matrixValues[3][3],
				zeroLapTemplates[jj+1].matrixValues[3][4]
			};
			float d=zeroLapTemplates[jj+1].diagonalR[3];
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(zeroLapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for( int i=s ; i<e ; i+=4 ) localXPtrs[2][i]=GaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,i)*d;
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-5;
				localXPtrs[2][i]=GaussSeidelUpdate3(zeroLapTemplates[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,i)*d;
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonalR[3];
			}
		}
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SolveInteriorReverse( int j , int c , int sSolve , int eSolve )
{
	int jj = Degree*3;

	{
		Pointer( float ) localXPtrs[ 2*Degree + 1 ];
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=0 ; y<2*Degree+1 ; y++ ) localXPtrs[y] = GetXRow( j-Degree+y , c );
		{
			ConstPointer( __m128 ) xPtrs[] = { ( Pointer( __m128 ) )localXPtrs[0] , ( Pointer( __m128 ) )localXPtrs[1] , ( Pointer( __m128 ) )localXPtrs[3] , ( Pointer( __m128 ) )localXPtrs[4] };

			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for(int i=-WordPerDegree+((sSolve-_start)>>2);i<((eSolve-_start)>>2)+WordPerDegree;i++)
			{
				localXAccum[0][i]=_mm_add_ps( xPtrs[0][i] , xPtrs[3][i] );
				localXAccum[1][i]=_mm_add_ps( xPtrs[1][i] , xPtrs[2][i] );
			}
		}
		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) xPtrs[]={ localXAccum[0] , localXAccum[1] , ( Pointer( __m128 ) )localXPtrs[2] };
		int s , e ;

		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2]
			};
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				dotSum = 0;
				int i = eSolve - _start - 1;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate3( zeroLapTemplates[jj+2].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[3];
				e = eSolve - _start - 5;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1];
				e = eSolve - _start - 1;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 7;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate3( mValues , xPtrs , localBPtr , dotSum , i );
			int i = 3;
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate3( zeroLapTemplates[jj].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonalR[3];
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2]
			};
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				dotSum = 0;
				int i = eSolve - _start - 2;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate2( zeroLapTemplates[jj+2].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[2];
				e = eSolve - _start - 6;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0];
				e = eSolve - _start - 2;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 6;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate2( mValues , xPtrs , localBPtr , dotSum , i );
			int i = 2;
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate2( zeroLapTemplates[jj].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonalR[2];
		}

		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2]
			};
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				SetInteriorDotSum( zeroLapTemplates[jj+2].matrixValues[1] , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				int i = eSolve - _start - 3;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[1];
				e = eSolve - _start - 7;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				e = eSolve - _start - 3;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 9;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i );
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				int i = 5;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate1( zeroLapTemplates[jj].matrixValues[1] , xPtrs , localBPtr , dotSum , i );
				i = 1;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonalR[1];			}
		}

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2]
			};
			float d = zeroLapTemplates[jj+1].diagonalR[0];

			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				SetInteriorDotSum( zeroLapTemplates[jj+2].matrixValues[0] , xPtrs , (eSolve - 1 - _start)>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				int i = eSolve - _start - 4;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[0];
				e = eSolve - _start - 8;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , ( eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				e = eSolve - _start - 4;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 8;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i );
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				int i = 4;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate0( zeroLapTemplates[jj].matrixValues[0] , xPtrs , localBPtr , dotSum , i );
				i = 0;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonalR[0];
			}
		}
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SolveReverse( int j , int c , int sSolve , int eSolve )
{
	int jj;
	if		(j<Degree)			jj=j;
	else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
	else						jj=Degree;

	if( periodicType==NO_PERIODIC )
	{
		if( sSolve<0 )		sSolve = 0;
		if( eSolve>major )	eSolve = major;
	}
	if( jj==Degree || periodicType==SPHERICAL_PERIODIC )	return SolveInteriorReverse( j , c , sSolve , eSolve );

	//	if( spherical ) jj=Degree;
	jj*=3;

	{
		Pointer( float ) localXPtrs[ 2*Degree + 1 ];
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=0 ; y<2*Degree+1 ; y++ )
			if( (j-Degree+y>=0 && j-Degree+y<minor) || periodicType==SPHERICAL_PERIODIC )	localXPtrs[y] = GetXRow( j-Degree+y , c );
			else													localXPtrs[y] = NullPointer< float >( );

		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) xPtrs[] = 
		{
			( Pointer( __m128 ) )localXPtrs[0],
			( Pointer( __m128 ) )localXPtrs[1],
			( Pointer( __m128 ) )localXPtrs[2],
			( Pointer( __m128 ) )localXPtrs[3],
			( Pointer( __m128 ) )localXPtrs[4]
		};
		int s , e ;

		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2],
				zeroLapTemplates[jj+1].matrixValues[3][3],
				zeroLapTemplates[jj+1].matrixValues[3][4]
			};
			float d = zeroLapTemplates[jj+1].diagonalR[3];

			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				dotSum = 0;
				int i = eSolve - _start - 1;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate3( zeroLapTemplates[jj+2].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[3];
				e = eSolve - _start - 5;
			}
			else
			{
				SetDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1];
				e = eSolve - _start - 1;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 7;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate3( mValues , xPtrs , localBPtr , dotSum , i ) * d;
			int i = 3;
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseGaussSeidelUpdate3( zeroLapTemplates[jj].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonalR[3];
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2],
				zeroLapTemplates[jj+1].matrixValues[2][3],
				zeroLapTemplates[jj+1].matrixValues[2][4]
			};
			float d = zeroLapTemplates[jj+1].diagonalR[2];

			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				dotSum = 0;
				int i = eSolve - _start - 2;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate2( zeroLapTemplates[jj+2].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[2];
				e = eSolve - _start - 6;
			}
			else
			{
				SetDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0];
				e = eSolve - _start - 2;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 6;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate2( mValues , xPtrs , localBPtr , dotSum , i ) * d;
			int i = 2;
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseGaussSeidelUpdate2( zeroLapTemplates[jj].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonalR[2];
		}

		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2],
				zeroLapTemplates[jj+1].matrixValues[1][3],
				zeroLapTemplates[jj+1].matrixValues[1][4]
			};
			float d = zeroLapTemplates[jj+1].diagonalR[1];

			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				SetDotSum( zeroLapTemplates[jj+2].matrixValues[1] , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				int i = eSolve - _start - 3;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[1];
				e = eSolve - _start - 7;
			}
			else
			{
				SetDotSum( mValues , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				e = eSolve - _start - 3;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 9;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i ) * d;
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				int i = 5;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate1( zeroLapTemplates[jj].matrixValues[1] , xPtrs , localBPtr , dotSum , i ) * d;
				i = 1;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonalR[1];
			}
		}

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2],
				zeroLapTemplates[jj+1].matrixValues[0][3],
				zeroLapTemplates[jj+1].matrixValues[0][4]
			};
			float d = zeroLapTemplates[jj+1].diagonalR[0];

			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				SetDotSum( zeroLapTemplates[jj+2].matrixValues[0] , xPtrs , (eSolve - 1 - _start)>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				int i = eSolve - _start - 4;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonalR[0];
				e = eSolve - _start - 8;
			}
			else
			{
				SetDotSum( mValues , xPtrs , ( eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				e = eSolve - _start - 4;
			}
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 8;
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i ) * d;
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				int i = 4;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate0( zeroLapTemplates[jj].matrixValues[0] , xPtrs , localBPtr , dotSum , i ) * d;
				i = 0;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonalR[0];
			}
		}
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SetInteriorResidual(int j,int c,int sSolve,int eSolve)
{
	int jj=Degree*3;
	{
		ConstPointer( float ) localXPtrs[2*Degree+1];
		Pointer( float )       localRPtr = GetRRow(j,c);
		ConstPointer( float ) localBPtr = GetBRow(j,c);
		for(int y=0;y<2*Degree+1;y++)	localXPtrs[y]=GetXRow(j-Degree+y,c);
		{
			ConstPointer( __m128 ) xPtrs[] =
			{
				( ConstPointer( __m128 ) )localXPtrs[0],
				( ConstPointer( __m128 ) )localXPtrs[1],
				( ConstPointer( __m128 ) )localXPtrs[3],
				( ConstPointer( __m128 ) )localXPtrs[4]
			};
			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for(int i=-WordPerDegree+((sSolve-_start)>>2);i<((eSolve-_start)>>2)+WordPerDegree;i++)
			{
				localXAccum[0][i]=_mm_add_ps(xPtrs[0][i],xPtrs[3][i]);
				localXAccum[1][i]=_mm_add_ps(xPtrs[1][i],xPtrs[2][i]);
			}
		}


		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) xPtrs[] = { localXAccum[0] , localXAccum[1] , ( ConstPointer( __m128 ) )localXPtrs[2]};
		int s,e;

		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localRPtr[0]=localBPtr[0]-GetInteriorLaplacianValue0(lapTemplates[jj].matrixValues[0],xPtrs,dotSum,0);
				s=4;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
			else							e=eSolve-_start;

			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue0(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-4;
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue0(lapTemplates[jj+2].matrixValues[0],xPtrs,dotSum,i);
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localRPtr[1]=localBPtr[1]-GetInteriorLaplacianValue1(lapTemplates[jj].matrixValues[1],xPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue1(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-3;
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue1(lapTemplates[jj+2].matrixValues[1],xPtrs,dotSum,i);
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue2(mValues,xPtrs,dotSum,i);
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-6;
				localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue2(lapTemplates[jj+2].matrixValues[2],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue3(mValues,xPtrs,dotSum,i);
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-5;
				localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue3(lapTemplates[jj+2].matrixValues[3],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SetResidual( int j , int c , int sSolve , int eSolve )
{
	int jj;
	if		(j<Degree)			jj=j;
	else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
	else						jj=Degree;

	if( periodicType==NO_PERIODIC )
	{
		if(sSolve<0)		sSolve=0;
		if(eSolve>major)	eSolve=major;
	}
	if( jj==Degree || periodicType==SPHERICAL_PERIODIC ) return SetInteriorResidual( j , c , sSolve , eSolve );
	jj*=3;

	{
		ConstPointer( float ) localXPtrs[2*Degree+1];
		Pointer( float )      localRPtr = GetRRow( j , c );
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=0 ; y<2*Degree+1 ; y++ )
			if(j-Degree+y>=0 && j-Degree+y<minor)	localXPtrs[y] = GetXRow(j-Degree+y,c);
			else									localXPtrs[y] = NullPointer< float >( );

		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) xPtrs[] =
		{
			( ConstPointer( __m128 ) )localXPtrs[0],
			( ConstPointer( __m128 ) )localXPtrs[1],
			( ConstPointer( __m128 ) )localXPtrs[2],
			( ConstPointer( __m128 ) )localXPtrs[3],
			( ConstPointer( __m128 ) )localXPtrs[4]
		};
		int s,e;
		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2],
				lapTemplates[jj+1].matrixValues[0][3],
				lapTemplates[jj+1].matrixValues[0][4]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localRPtr[0]=localBPtr[0]-GetLaplacianValue0(lapTemplates[jj].matrixValues[0],xPtrs,dotSum,0);
				s=4;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
			if( eSolve==major && periodicType==NO_PERIODIC)	e=eSolve-_start-4;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetLaplacianValue0(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-4;
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetLaplacianValue0(lapTemplates[jj+2].matrixValues[0],xPtrs,dotSum,i);
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2],
				lapTemplates[jj+1].matrixValues[1][3],
				lapTemplates[jj+1].matrixValues[1][4]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				localRPtr[1]=localBPtr[1]-GetLaplacianValue1(lapTemplates[jj].matrixValues[1],xPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetLaplacianValue1(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-3;
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetLaplacianValue1(lapTemplates[jj+2].matrixValues[1],xPtrs,dotSum,i);
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2],
				lapTemplates[jj+1].matrixValues[2][3],
				lapTemplates[jj+1].matrixValues[2][4]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue2(mValues,xPtrs,dotSum,i);
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-6;
				localRPtr[i]=localBPtr[i]-GetLaplacianValue2(lapTemplates[jj+2].matrixValues[2],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2],
				lapTemplates[jj+1].matrixValues[3][3],
				lapTemplates[jj+1].matrixValues[3][4]
			};
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue3(mValues,xPtrs,dotSum,i);
			if( eSolve==major && periodicType==NO_PERIODIC )
			{
				int i=eSolve-_start-5;
				localRPtr[i]=localBPtr[i]-GetLaplacianValue3(lapTemplates[jj+2].matrixValues[3],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Solve( void )
{
	int idx = index+iters*Degree-1;
	int start , stop;

	stop = index+iters*Degree-1;
	start = stop-(iters-1)*Degree;	// = index + Degree - 1;
	if( periodicType==SPHERICAL_PERIODIC )
	{
		if( stop==0 )
		{
			// This is necessary, because when the same resolution row is solved twice, we want to reuse the solution
			for( int d=0 ; d<Degree ; d++ ) SyncSolverHead( d , false ) , SyncSolverHead( d, true );
		}
		if( stop==minor-Degree )
		{
			for( int d=minor-Degree ; d<minor ; d++ ) SyncSolverTail( d , false ) , SyncSolverTail( d , true );
			SyncSolverRight( stop , false ) , SyncSolverLeft ( stop , true );
		}
	}
	if( stop==0 ) SyncSolverRight( -1 , false ) , SyncSolverLeft( -1 , true );
	{
		int offset;
		if( stop&1 ) // Solve right to left
		{
			offset = _padSize-2*WordPerDegree*RealPerWord;
			for( int i=stop ; i>=start ; i-=Degree )
			{
				if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) SolveReverse( i , c , _end-_padSize-offset , _end+offset );
				offset -= WordPerDegree*RealPerWord;
			}
			SyncSolverRight( stop , false );
			for(int c=0;c<Channels;c++)
			{
				int s128 = _end128-_start128-3*_padSize128;
				for(int l=0;l<laneNum;l++)
				{
					int b128 = _start128+2*_padSize128+((laneNum-l-1)*s128)/laneNum;
					int e128 = _start128+2*_padSize128+((laneNum-l  )*s128)/laneNum;
					int b = b128<<2;
					int e = e128<<2;

					offset=_padSize-2*WordPerDegree*RealPerWord;
					for(int i=stop;i>=start;i-=Degree)
					{
						if(i>=0 && i<minor)	SolveReverse( i , c , b-offset , e-offset );
						offset -= WordPerDegree*RealPerWord;
					}
				}
			}
			SyncSolverLeft( stop , true );
			offset=_padSize-2*WordPerDegree*RealPerWord;
			for( int i=stop ; i>=start ; i-=Degree )
			{
				if( leftStream )
				{
					if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) SolveReverse( i , c , _start+offset , _start+2*_padSize-offset );
				}
				else
				{
					if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) SolveReverse( i , c , _start        , _start+2*_padSize-offset );
				}
				offset -= WordPerDegree*RealPerWord;
			}
		}
		else // Solve left to right
		{
			// Solve left
			offset=_padSize-2*WordPerDegree*RealPerWord;
			for( int i=stop ; i>=start ; i-=Degree )
			{
				if(i>=0 && i<minor)	for(int c=0;c<Channels;c++)	Solve( i , c , _start-offset , _start+_padSize+offset );
				offset -= WordPerDegree*RealPerWord;
			}
			// Write out to right buffer
			SyncSolverLeft( stop , false );
			// Solve center
			for( int c=0 ; c<Channels ; c++ )
			{
				int s128 = _end128-_start128-3*_padSize128;
				for(int l=0;l<laneNum;l++)
				{
					int b128 = _start128+_padSize128+( l   *s128)/laneNum;
					int e128 = _start128+_padSize128+((l+1)*s128)/laneNum;
					int b = b128<<2;
					int e = e128<<2;

					offset=_padSize-2*WordPerDegree*RealPerWord;
					for(int i=stop;i>=start;i-=Degree)
					{
						if(i>=0 && i<minor)	Solve( i , c , b+offset , e+offset );
						offset -= WordPerDegree*RealPerWord;
					}
				}
			}
			// Read in right buffer
			SyncSolverRight( stop , true );
			// Solve Right
			offset=_padSize-2*WordPerDegree*RealPerWord;
			for( int i=stop ; i>=start ; i-=Degree )
			{
				if( rightStream )
				{
					if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) Solve( i , c , _end-2*_padSize+offset , _end-offset );
				}
				else
				{
					if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) Solve( i , c , _end-2*_padSize+offset , _end );
				}
				offset -= WordPerDegree*RealPerWord;
			}
		}
	}
	if( periodicType==SPHERICAL_PERIODIC )
	{
		if(start<Degree && stop>=0 )
		{
			int d = start;
			while( d<0 ) d += Degree;
			SyncSolverHead( d , false ) , SyncSolverHead( d , true );
			SyncSolverRight( start , false , false ) , SyncSolverLeft ( start , false , false );
			SyncSolverRight( start , true  , false ) , SyncSolverLeft ( start , true  , false );
		}
		if(stop>minor-Degree-1 && start<minor)
		{
			int d = start;
			while( d<minor-Degree ) d += Degree;
			SyncSolverTail( d , false ) , SyncSolverTail( d , true );
		}
	}
	idx=index+OffsetR();
	if( idx>=0 && idx<minor && setResidual )
	{
		for( int c=0 ; c<Channels ; c++ )
		{
			SetResidual(idx,c,_start-WordPerDegree*RealPerWord,_start+_padSize);
			int s128 = _end128-_start128-2*_padSize128;
			for(int l=0;l<laneNum;l++)
			{
				int b128 = _start128+_padSize128+((l  )*s128)/laneNum;
				int e128 = _start128+_padSize128+((l+1)*s128)/laneNum;
				int b = b128<<2;
				int e = e128<<2;
				SetResidual( idx , c , b , e );
			}
			SetResidual( idx , c , _end-_padSize , _end+WordPerDegree*RealPerWord );
		}
		if( idx<Degree )
		{
			SyncResidualHead( idx , false , false) , SyncResidualHead( idx , true ,  false);
			SyncResidualHead( idx , false , true ) , SyncResidualHead( idx , true ,  true );
		}
		if(idx<minor && idx>=minor-Degree)
		{
			SyncResidualTail( idx , false , false) , SyncResidualTail( idx , true , false );
			SyncResidualTail( idx , false , true ) , SyncResidualTail( idx , true ,  true );
		}

		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) localXPtr = GetXRow( idx , c );
			Pointer( float ) localBPtr = GetBRow( idx , c );
			Pointer( float ) localRPtr = GetRRow( idx , c );
			for( int i=0 ; i<_size ; i++ )
			{
				xSquareNorm += (double)(localXPtr[i])*localXPtr[i] * laplacianScale * laplacianScale;
				rSquareNorm += (double)(localRPtr[i])*localRPtr[i];
				bSquareNorm += (double)(localBPtr[i])*localBPtr[i];
				solutionSum[c] += (double)(localXPtr[i]);
			}
		}
	}
	if( showProgress )
	{
		if( !idx )
		{
			progressCount = 0;
			pastIndex = 0;
			startProgressTime = pastProgressTime = Time();
		}
		if( idx>=0 && idx<minor )
		{
			size_t current , peak;
			WorkingSetInfo(current,peak);
			double t = Time();
			int pCount = (idx*1000) / ( minor-1 );
			if( pCount>progressCount )
			{
				printf("[%.1f%%]  [%d/%d] [%4llu/%4llu MB] Rows/Second: %f\t(%f)         \r" , float(idx)/float(minor-1)*100 , idx , minor , (unsigned long long )(current>>20) , (unsigned long long)(peak>>20) , float(idx-pastIndex) / float( t-pastProgressTime ) , float(idx) / float( t-startProgressTime ) );
				progressCount = pCount;
				pastProgressTime = t;
				pastIndex = idx;
			}
		}
		else if( idx==minor ) printf( "\n" ) , fflush( stdout );
	}
}

//////////////////////////////////////
// SocketedMultiGridStreamingSolver //
//////////////////////////////////////
template< int Channels , class StorageType , class SyncType >
SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SocketedMultiGridStreamingSolver( void )
{
	_X = _B = NULL;
	inX = outX = inB = outB = outR = outP = NULL;
	scratchR = NullPointer< float >( );
	scratchP = NullPointer< float >( );

	parent = NULL;
	rChild = NULL;
	pChild = NULL;
	localRAccum = NullPointer< __m128 >( );
	prolongationStencil = NullPointer< ProlongationStencilSSE2 >( );

	int halfD=((Degree+1)>>1);
	int dual=(Degree&1)?0:1;

	// Restriction
	// To solve X[i] at depth d-1:
	//		Need to have B[i+iters-1] at depth d-1:
	//		Need to have R[2*(i+iters-1)-halfD,2*(i+iters-1)+halfD+dual]

	// If we update i at depth d:
	// @ depth d+1 i gets prolonged to: 	[2*i-halfD,2*i+halfD+dual]
	// @ depth d-1 i gets restricted to:	[(i-halfD-dual+1)/2,(i+halfD)/2]
	// Can start processing the parent once:
	//		i > dual+halfD-1
	// Can start processing the child once:
	//		i >= halfD/2


	// Restriction

	// Can start processing the parent once ((i+OffsetR())+1) gets restricted to a starting index bigger than zero:
	//		(i+OffsetR()+1-halfD-dual+1)/2 >= 1
	//		i+OffsetR()-halfD-dual >= 0
	//		i >= dual+halfD-OffsetR()
	startRestriction = dual + halfD - OffsetR();
	// In the coarser resolution, we only care that the B value has been set from the finer residual
	startRestriction -= 2*Degree;
	// Get the parity of the starting index so we know on what parity to pass control to the coarser resolution
	restrictionBit=startRestriction&1;


	// Prolongation
	// Can start processing the child once i gets prolonged to an index bigger than or equal to one:
	//		2*i-halfD >= 1
	//		i >= (halfD+1)/2
	startProlongation = (halfD+2)>>1;
	// If we prolong at idx=startProlongation, in the finer resolution it will map to a leading index of
	prolongationOffset = 2*startProlongation+halfD+dual;

	_B = _X = NULL;
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::Initialize( double iWeight , bool lump , double gWeight , int start , int end , int major , int minor , int iters ,
																		DataStream* leftStream,Socket* syncSockets,DataStream* rightStream,
																		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
																		)
{
#if 0
printf( "Initializing without stencils: %d\n" , minor );
#endif
	// Warning!!! A certaint amount of care needs to be taken here. It's safe to call this function when iWeight == 0 since
	// the Laplacian stencil is independent of the resolution. However, the dot-product stencil is not!!!
	// Care should be taken:
	// 1] When transitioning from out-of-core to in-core
	// 2] When transitioning from in-core to a conjugate-gradient solver
	// 3] When collapsing multiple threads into one

	DotProductStencil dotMajor , d2DotMajor , dotMinor , d2DotMinor;

#if 1
	SetDownSampledStencil( major , 0 , dotMajor , d2DotMajor , lump );
	SetDownSampledStencil( minor , 0 , dotMinor , d2DotMinor , lump );
#else
	FiniteElements1D< double , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( major , dotMajor   , 0 , 0 );
	FiniteElements1D< double , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( major , d2DotMajor , 1 , 1 , false );
	FiniteElements1D< double , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( minor , dotMinor   , 0 , 0 );
	FiniteElements1D< double , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( minor , d2DotMinor , 1 , 1 , false );

	if( lump )
	{
		for( int i=0 ; i<=2*Degree ; i++ )
		{
			double majorSum = 0 , minorSum = 0;
			for( int j=0 ; j<=2*Degree ; j++ )
			{
				majorSum += dotMajor.caseTable[i].values[j] , dotMajor.caseTable[i].values[j] = 0;
				minorSum += dotMinor.caseTable[i].values[j] , dotMinor.caseTable[i].values[j] = 0;
			}
			dotMajor.caseTable[i].values[Degree] = majorSum;
			dotMinor.caseTable[i].values[Degree] = minorSum;
		}
	}
#endif

	Initialize( dotMajor , d2DotMajor , dotMinor , d2DotMinor , iWeight , gWeight , start , end , major , minor , iters ,
		leftStream , syncSockets , rightStream ,
		memoryMappedFile , periodicType , server
		);
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::Initialize( DotProductStencil& dotMajor , DotProductStencil& d2DotMajor ,
																		DotProductStencil& dotMinor , DotProductStencil& d2DotMinor,
																		double iWeight , double gWeight ,
																		int start , int end , int major , int minor , int iters ,
																		DataStream* leftStream , Socket* syncSockets , DataStream* rightStream,
																		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
																		)
{
	this->dotMajor   = dotMajor;
	this->dotMinor   = dotMinor;
	this->d2DotMajor = d2DotMajor;
	this->d2DotMinor = d2DotMinor;

#if 0
#pragma omp critical
{
	printf( "Initializing: %d\n" , minor );
	for( int i=0 ; i<=2*Degree ; i++ ) printf( "\t[%d] %g %g %g %g\n" , i , dotMajor.caseTable[i].values[Degree] , d2DotMajor.caseTable[i].values[Degree] , dotMinor.caseTable[i].values[Degree] , d2DotMinor.caseTable[i].values[Degree] );
}
#endif

	MatrixStencil lStencil;
	for( int i=0 ; i<=2*Degree ; i++ ) for( int j=0 ; j<=2*Degree ; j++ )
		for( int k=0 ; k<=2*Degree ; k++ ) for( int l=0 ; l<=2*Degree ; l++ )
			lStencil.caseTable[j][i].values[l][k] =
			( dotMajor.caseTable[i].values[k] * d2DotMinor.caseTable[j].values[l] + d2DotMajor.caseTable[i].values[k] * dotMinor.caseTable[j].values[l] ) * gWeight +
			( dotMajor.caseTable[i].values[k] * dotMinor.caseTable[i].values[l] ) * iWeight;

	if( syncSockets )
		SocketedStreamingSolver< Channels , SyncType >::Init( lStencil , start , end , major , minor , iters ,
		leftStream , syncSockets[0] , rightStream ,
		periodicType , server
		);
	else
		SocketedStreamingSolver< Channels , SyncType >::Init( lStencil , start , end , major , minor , iters ,
		leftStream , _INVALID_SOCKET_ , rightStream ,
		periodicType , server
		);

	if( _B ) delete _B;
	if( _X ) delete _X;
	_B = _X = NULL;

	if( memoryMappedFile )
	{
		{
			_B = new MultiStreamIOClient( _size*Channels*sizeof(StorageType) , minor , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
			_X = new MultiStreamIOClient( _size*Channels*sizeof(StorageType) , minor , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
		}
		inCore = false;
	}
	else
	{
		_B = new MemoryBackedGrid( _size * Channels * sizeof( StorageType ) , minor );
		_X = new MemoryBackedGrid( _size * Channels * sizeof( StorageType ) , minor );
		inCore = true;
	}
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( major>>1 , majorProlongationStencil , major );
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( minor>>1 , minorProlongationStencil , minor );
	int foo;
	FiniteElements1D< double , Type , Degree >::RestrictionStencil( major , majorRestrictionStencil , foo );
	FiniteElements1D< double , Type , Degree >::RestrictionStencil( minor , minorRestrictionStencil , foo );
	if( parent )
	{
		DotProductStencil newDotMajor , newD2DotMajor , newDotMinor , newD2DotMinor;
		CombineStencils< double , Type , Degree >(   dotMajor , majorProlongationStencil , major ,   newDotMajor );
		CombineStencils< double , Type , Degree >( d2DotMajor , majorProlongationStencil , major , newD2DotMajor );
		CombineStencils< double , Type , Degree >(   dotMinor , minorProlongationStencil , minor ,   newDotMinor );
		CombineStencils< double , Type , Degree >( d2DotMinor , minorProlongationStencil , minor , newD2DotMinor );


		if( syncSockets )
			parent->Initialize( newDotMajor , newD2DotMajor , newDotMinor , newD2DotMinor , iWeight , gWeight , start>>1 , end>>1 , major>>1 , minor>>1 , iters ,
			leftStream , syncSockets+1 , rightStream ,
			memoryMappedFile , periodicType , _server
			);
		else
			parent->Initialize( newDotMajor , newD2DotMajor , newDotMinor , newD2DotMinor , iWeight , gWeight , start>>1 , end>>1 , major>>1 , minor>>1 , iters ,
			leftStream , NULL , rightStream ,
			memoryMappedFile , periodicType , _server
			);
	}
}
template< int Channels , class StorageType , class SyncType >
SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::~SocketedMultiGridStreamingSolver(void)
{
	if( localRAccum )
	{
		localRAccum -= _padSize128;
		AlignedFreePointer( localRAccum );
	}
	AlignedFreePointer( prolongationStencil );
	if( _X ) delete _X;
	if( _B ) delete _B;
	_X = _B = NULL;
	AlignedFreePointer( scratchR );
	AlignedFreePointer( scratchP );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::InitProlongation(void)
{
	if( parent ) parent->InitProlongation();
	if( inX )  inX->reset ( true  , 1 );
	if( inB )  inB->reset ( true  , 1 );
	if( outX ) outX->reset( false , 1 );
	if( outB ) outB->reset( false , 1 );
	if( outP ) outP->reset( false , 1 );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::InitRestriction( void )
{
	{
		localRAccum = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
		localRAccum += _padSize128;
		ALIGN( float scratch[4] , 16 );
		prolongationStencil = AlignedAllocPointer< ProlongationStencilSSE2 >( 3 , ALIGNMENT );
		for(int j=0;j<3;j++)
		{
			int jj;
			if(j==0)	jj=0;
			else		jj=Degree;
			scratch[0] = float( majorProlongationStencil.caseTable[jj].values[1] );
			scratch[1] = float( majorProlongationStencil.caseTable[jj].values[2] );
			scratch[2] = float( majorProlongationStencil.caseTable[jj].values[3] );
			if(j!=2)	scratch[3]=float( majorProlongationStencil.caseTable[Degree].values[0] );
			else		scratch[3]=float( 0 );
			prolongationStencil[j].matrixValues[0]=_mm_load_ps(scratch);

			if(j==2)	jj=2*Degree;
			else		jj=Degree;
			if(j!=0)	scratch[0]=float( majorProlongationStencil.caseTable[Degree].values[3] );
			else		scratch[0]=float( 0 );
			scratch[1] = float( majorProlongationStencil.caseTable[jj].values[0] );
			scratch[2] = float( majorProlongationStencil.caseTable[jj].values[1] );
			scratch[3] = float( majorProlongationStencil.caseTable[jj].values[2] );
			prolongationStencil[j].matrixValues[1]=_mm_load_ps(scratch);
		}
	}
	if( parent ) parent->InitRestriction();
	if( inX )  inX->reset ( true  , 1 );
	if( inB )  inB->reset ( true  , 1 );
	if( outX ) outX->reset( false , 1 );
	if( outB ) outB->reset( false , 1 );
	if( outR ) outR->reset( false , 1 );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetProlongation( void )
{
	// Zig-zagging means that we are not synchronized until the subsequent row
//	SocketedStreamingSolver< Channels , SyncType >::Set( -OffsetX(iters)   , 0 , 0 , 0 , 0 );
	SocketedStreamingSolver< Channels , SyncType >::Set( -OffsetX(iters)-1 , 0 , 0 , 0 , 0 );
	// Set the iterations for the child so that it trails accordingly
	if( pChild ) pChild->SetProlongation();
	{
		if( _B ) _B->reset( true , bSize );
		if( _X ) _X->reset( true , xSize );
		if( _B ) _B->SetServer( _server );
		if( _X ) _X->SetServer( _server );
	}
	AlignedFreePointer( scratchP );
	if( outP ) scratchP = AlignedAllocPointer< float >( _size*2 * Channels , ALIGNMENT );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetRestriction( void )
{
	// Warning: This is actually somewhat agressive as you only really need something closer Size/2
//	SocketedStreamingSolver< Channels , SyncType >::Set( -OffsetX(iters) , 0 , 0 , 0 , 0 ,   FiniteElements1D< float , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size      -OffsetR() );
	SocketedStreamingSolver< Channels , SyncType >::Set( -OffsetX(iters) , 0 , 0 , 0 , 0 , ( FiniteElements1D< float , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size+1 )/2-OffsetR() );
	// Set the iterations for the parent so that it trails accordingly
	if( parent ) parent->SetRestriction();
	{
		if( _B ) _X->reset( false , xSize );
		if( _X ) _B->reset( false , bSize );
		if( _B ) _B->SetServer( _server );
		if( _X ) _X->SetServer( _server );
	}
	AlignedFreePointer( scratchR );
	if( outR ) scratchR = AlignedAllocPointer< float >( _size / 2 * Channels , ALIGNMENT );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::UnSetProlongation(void)
{
	SocketedStreamingSolver< Channels , SyncType >::UnSet();
	{
		if( _X ) _X->unset();
		if( _B ) _B->unset();
	}
	AlignedFreePointer( prolongationStencil );
	AlignedFreePointer( scratchP );
	if( pChild ) pChild->UnSetProlongation();
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::UnSetRestriction(void)
{
	SocketedStreamingSolver< Channels , SyncType >::UnSet();
	if( _X ) _X->unset();
	if( _B ) _B->unset();
	AlignedFreePointer( scratchR );
	if( parent ) parent->UnSetRestriction();
}

template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	int highStart = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(sRestrict);
	int highEnd = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(eRestrict)+FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;
	int highStart128 = (highStart-(RealPerWord-1))/RealPerWord;
	int highEnd128   = (highEnd  +(RealPerWord-1))/RealPerWord;
	highStart128 -= _start128;
	highEnd128 -= _start128;

	if( (idx>=Degree && idx<(minor>>1)-Degree) || periodicType==SPHERICAL_PERIODIC ) SetInteriorRestriction( lB , c , idx , sRestrict , eRestrict );
	else
	{
		int startY=FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(idx);
		int jj;
		if		(idx<Degree)				jj=idx;
		else if	(idx>(minor>>1)-1-Degree)	jj=2*Degree+(idx-((minor>>1)-1));
		else								jj=Degree;

		{
			ConstPointer( __m128 ) localRPtrs[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
			for(int y=0;y<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;y++)
				if(startY+y>=0 && startY+y<minor)	localRPtrs[y] = ( Pointer( __m128 ) )SocketedStreamingSolver< Channels , SyncType >::GetRRow(startY+y,c);
				else								localRPtrs[y] = NullPointer< __m128 >( );
			ALIGN( float scratch[4] , 16 );
			__m128 res[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
			float dotSum;
			__m128 dSum;
			int s,e;

			for(int d=0;d<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;d++)
			{
				for(int j=0;j<4;j++)	scratch[j] = float( minorProlongationStencil.caseTable[jj].values[d] );
				res[d]=_mm_load_ps(scratch);
			}
			memset( localRAccum+highStart128 , 0 , sizeof(__m128)*(highEnd128-highStart128 ) );
			for(int d=0;d<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;d++)
				if(localRPtrs[d])
					for(int i=highStart128;i<highEnd128;i++) localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],localRPtrs[d][i]));

			// Offset 0
			{
				if( sRestrict==0 && periodicType==NO_PERIODIC )
				{
					dotSum=0;
					lB[0] = RestrictionUpdate0( prolongationStencil[0].matrixValues[0] , localRAccum , dotSum , 0 );
					s=2;
				}
				else
				{
					SetRestrictionDotSum(prolongationStencil[1].matrixValues[0],localRAccum,highStart128,dSum);
					_mm_store_ps(scratch,dSum);
					dotSum=scratch[3];
					s=sRestrict-_start/2;
				}

				if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-2;
				else									e=(eRestrict-_start/2);

				for(int i=s;i<e;i+=2)						lB[i]=RestrictionUpdate0(prolongationStencil[1].matrixValues[0],localRAccum,dotSum,i);
				int i=(eRestrict-_start/2)-2;
				if( eRestrict==major/2 && periodicType==NO_PERIODIC )	lB[i]=RestrictionUpdate0(prolongationStencil[2].matrixValues[0],localRAccum,dotSum,i);
			}
			// Offset 1
			{
				if( sRestrict==0 && periodicType==NO_PERIODIC )
				{
					SetRestrictionDotSum(prolongationStencil[0].matrixValues[1],localRAccum,0,dSum);
					s=1;
				}
				else
				{
					SetRestrictionDotSum(prolongationStencil[1].matrixValues[1],localRAccum,highStart128+1,dSum);
					s=(sRestrict-_start/2)+1;
				}
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[1]+scratch[2]+scratch[3];
				if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-4;
				else									e=(eRestrict-_start/2);
				for(int i=s;i<e;i+=2)	lB[i]=RestrictionUpdate1(prolongationStencil[1].matrixValues[1],localRAccum,dotSum,i);
				if( eRestrict==major/2 && periodicType==NO_PERIODIC )
				{
					int i=(eRestrict-_start/2)-3;
					lB[i]=RestrictionUpdate1(prolongationStencil[2].matrixValues[1],localRAccum,dotSum,i);
					i+=2;
					lB[i]=dotSum;
				}
			}
		}
	}
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetInteriorRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	int startY=FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(idx);
	int highStart = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(sRestrict);
	int highEnd = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(eRestrict)+FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;
	int highStart128 = (highStart-(RealPerWord-1))/RealPerWord;
	int highEnd128   = (highEnd  +(RealPerWord-1))/RealPerWord;
	highStart128 -= _start128;
	highEnd128 -= _start128;
	{
		ConstPointer( __m128 ) localRPtrs[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
		for(int y=0;y<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;y++)
			localRPtrs[y] = ( Pointer( __m128 ) )SocketedStreamingSolver< Channels , SyncType >::GetRRow(startY+y,c);

		ALIGN( float scratch[4] , 16 );
		__m128 res[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
		float dotSum;
		__m128 dSum;
		int s,e;
		for(int d=0;d<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;d++)
		{
			for(int j=0;j<4;j++)	scratch[j] = float( minorProlongationStencil.caseTable[Degree].values[d] );
			res[d]=_mm_load_ps(scratch);
		}
		{
			int d=0;
			for(int i=highStart128;i<highEnd128;i++)	localRAccum[i]=_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i]));
		}
		for(int d=1;d<=(Degree>>1);d++)
			for(int i=highStart128;i<highEnd128;i++)	localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i])));
		if(Degree&1)
		{
			int d=(Degree+1)>>1;
			for(int i=highStart128;i<highEnd128;i++)	localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i])));
		}
		// Offset 0
		{
			if( sRestrict==0 && periodicType==NO_PERIODIC )
			{
				dotSum=0;
				lB[0] = RestrictionUpdate0(prolongationStencil[0].matrixValues[0],localRAccum,dotSum,0);
				s=2;
			}
			else
			{
				SetRestrictionDotSum(prolongationStencil[1].matrixValues[0],localRAccum,highStart128,dSum);
				_mm_store_ps(scratch,dSum);
				dotSum = scratch[3];
				s = sRestrict-_start/2;
			}
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-2;
			else									e=(eRestrict-_start/2);

			for(int i=s;i<e;i+=2)					lB[i]=RestrictionUpdate0(prolongationStencil[1].matrixValues[0],localRAccum,dotSum,i);
			int i=(eRestrict-_start/2)-2;
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )	lB[i]=RestrictionUpdate0(prolongationStencil[2].matrixValues[0],localRAccum,dotSum,i);
		}
		// Offset 1
		{
			if( sRestrict==0 && periodicType==NO_PERIODIC )
			{
				SetRestrictionDotSum(prolongationStencil[0].matrixValues[1],localRAccum,0,dSum);
				s=1;
			}
			else
			{
				SetRestrictionDotSum(prolongationStencil[1].matrixValues[1],localRAccum,highStart128+1,dSum);
				s=(sRestrict-_start/2)+1;
			}
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-4;
			else									e=(eRestrict-_start/2);
			for(int i=s;i<e;i+=2)	lB[i]=RestrictionUpdate1(prolongationStencil[1].matrixValues[1],localRAccum,dotSum,i);
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )
			{
				int i=(eRestrict-_start/2)-3;
				lB[i]=RestrictionUpdate1(prolongationStencil[2].matrixValues[1],localRAccum,dotSum,i);
				i+=2;
				lB[i]=dotSum;
			}
		}
	}
}

template< int Channels , class StorageType , class SyncType >
bool SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::IterateRestriction( void )
{
	// index+OffsetR is the last index at which the residual is set after the solve.
	int restrictionIndex = FiniteElements1D< float , Type , Degree >::FullRestrictionStencil::RestrictionStencil::Start( index+OffsetR() );
	if( restrictionIndex>=minor/2 ) return false;
	int idx=index+OffsetB(iters);

	// If we have an initial guess, use it.
	if( idx+Degree>=0 && idx+Degree<minor )
	{
		if( inX )
		{
			Pointer( SyncType ) inPtr = ( Pointer( SyncType ) )(*inX)[idx+Degree];
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) x = SocketedStreamingSolver< Channels , SyncType >::GetXRow( idx+Degree , c );
				for( int i=0 ; i<_size ; i++ ) x[i] = float( inPtr[i*Channels+c] ) * laplacianScaleR;
			}
			inX->advance();
		}
	}

	// For the constraints, we can either get them directly, or we can get them by restricting from a child node
	if( idx>=0 && idx<minor )
	{
		if( inB )
		{
			Pointer( SyncType ) inPtr = ( Pointer( SyncType ) )(*inB)[idx];
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) bRow = SocketedStreamingSolver< Channels , SyncType >::GetBRow( idx , c );
				for( int i=0 ; i<_size ; i++ ) bRow[i] = float( inPtr[i*Channels+c] );
			}
			inB->advance();
		}
		else if( rChild )
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) bRow = SocketedStreamingSolver< Channels , SyncType >::GetBRow( idx , c );
				for( int l=0 ; l<laneNum ; l++ )
				{
					int b128 = _start128+((l  )*_size128)/laneNum;
					int e128 = _start128+((l+1)*_size128)/laneNum;
					int b = b128<<2;
					int e = e128<<2;
					rChild->SetRestriction( bRow , c , idx , b , e );
				}
			}
		}
		else fprintf( stderr , "[ERROR] SocketedMultiGridStreamingSolver::IterateRestriction: no input constraints @ [%d x %d]\n" , major , minor ) , exit( 0 );
	}

	// Run an iteration of the Gauss-Seidel solver
	if( idx>=0 ) SocketedStreamingSolver< Channels , SyncType >::Solve();

	// Write out the current solution
	if( SocketedStreamingSolver< Channels , SyncType >::template UpdateXOutput< StorageType >( _X ) ) _X->advance();
	if( SocketedStreamingSolver< Channels , SyncType >::template UpdateBOutput< StorageType >( _B ) ) _B->advance();

	// Write out the constraints
	if( outB && index>=0 && index<minor )
	{
		Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outB)[index];
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) inP = SocketedStreamingSolver< Channels , SyncType >::GetBRow( index , c );
			for( int i=0 ; i<_size ; i++ ) outPtr[i*Channels+c] = SyncType( inP[i] );
		}
		outB->advance();
	}

	// Write out the solution
	if( outX && index>=0 && index<minor )
	{
		Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outX)[index];
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) inP = SocketedStreamingSolver< Channels , SyncType >::GetXRow( index , c );
			for( int i=0 ; i<_size ; i++ ) outPtr[i*Channels+c] = SyncType( inP[i]*laplacianScale );
		}
		outX->advance();
	}

	// Output the restricted residual
	if( setResidual && outR && restrictionIndex>=0 && (index&1)==restrictionBit )
	{
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) outPtr = scratchR + c*_size/2;
			for( int l=0 ; l<laneNum ; l++ )
			{
				int b128 = (_start128/2) + ( (l  ) * (_size128/2) ) / laneNum;
				int e128 = (_start128/2) + ( (l+1) * (_size128/2) ) / laneNum;
				int b = (_start/2) + ( (l  ) * (_size/2) ) / laneNum;;
				int e = (_start/2) + ( (l+1) * (_size/2) ) / laneNum;
				SetRestriction( outPtr , c , restrictionIndex , b , e );
			}
		}
		Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outR)[restrictionIndex];
		for( int c=0 ; c<Channels ; c++ ) 
		{
			Pointer( float ) _outPtr = scratchR + c * _size/2;
			for( int i=0 ; i<_size/2 ; i++ ) outPtr[i*Channels+c] = SyncType( _outPtr[i] );
		}
		outR->advance();
	}

	// If we have done two iterations, have the parent do one
	if( parent ) if( index>=startRestriction && (index&1)==restrictionBit ) parent->IterateRestriction();
	return Increment();
}

template< int Channels , class StorageType , class SyncType >
bool SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::IterateProlongation( void )
{
	int prolongationIndex = FiniteElements1D< float , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( index + Degree - 1  );
	if( prolongationIndex>=minor*2 && index+OffsetR()>=minor ) return false;
	int idx = index + OffsetX( iters );
	if( idx>=0 )
	{
		if( idx<minor )
		{
			// Get the (external) solution
			if( inX )
			{
				Pointer( SyncType ) inPtr = ( Pointer( SyncType ) )(*inX)[idx];
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) x = SocketedStreamingSolver< Channels , SyncType >::GetXRow( idx , c );
					for( int i=0 ; i<_size ; i++ ) x[i] += float( inPtr[i*Channels+c] );
				}
				inX->advance();
			}
			else if( parent )
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) x = SocketedStreamingSolver< Channels , SyncType >::GetXRow( idx , c );
					parent->SetProlongation( x , c , idx , _start , _size , major , minor );
				}
			}

			// Rescale the external estimate of the solution into the internal
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) x = SocketedStreamingSolver< Channels , SyncType >::GetXRow( idx , c );
				for( int i=0 ; i<_size ; i++ ) x[i] *= laplacianScaleR;
			}
		}

		// Add in the (internal) solution and constraints from the restriction phase
		if( SocketedStreamingSolver< Channels , SyncType >::template UpdateBInput< StorageType >( _B ) ) _B->advance();
		if( SocketedStreamingSolver< Channels , SyncType >::template UpdateXInput< StorageType >( _X ) ) _X->advance();

		// Run an iteration of the Gauss-Seidel solver
		SocketedStreamingSolver< Channels , SyncType >::Solve();
	}

	// Output the solution
	if( index>=0 && index<minor && outX )
	{
		Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outX)[index];
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) x = SocketedStreamingSolver< Channels , SyncType >::GetXRow( index , c );
			for( int i=0 ; i<_size ; i++ ) outPtr[i*Channels+c] = SyncType( x[i]*laplacianScale );
		}
		outX->advance();
	}

	// Output the prolonged solution
	if( outP )
		for( int off=0 ; off<2 ; off++ )
			if( prolongationIndex+off>=0 && prolongationIndex+off<2*minor )
			{
				Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )( (*outP)[prolongationIndex+off] );
				memset( scratchP , 0 , _size * 2 * Channels * sizeof( float ) );
				for( int c=0 ; c<Channels ; c++ ) SetProlongation( scratchP + c*_size*2 , c , prolongationIndex+off , _start*2 , _size*2 , major*2 , minor*2 );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) _outPtr =  scratchP + c*_size*2;
					for( int i=0 ; i<_size*2 ; i++ ) outPtr[i*Channels+c] = SyncType( _outPtr[i] );
				}
				outP->advance();
			}

	// Advance the the child two steps
	if( pChild && index>=startProlongation ) if( pChild->IterateProlongation() ) pChild->IterateProlongation();
	return Increment();
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetProlongation( Pointer( float ) highX , int c , int highIdx , int highStart , int highSize , int highMajor , int highMinor )
{
	int startY = RestrictionStencil::RestrictionStencil::Start( highIdx );
	int jj;
	if( periodicType!=SPHERICAL_PERIODIC )
		if		( highIdx<Degree )				jj = highIdx;
		else if	( highIdx>highMinor-1-Degree )	jj = 2*Degree+1+(highIdx-(highMinor-1));
		else if ( (highIdx-Degree)&1 )			jj = Degree+1;
		else									jj = Degree;
	else
		if( (highIdx-Degree)&1 )	jj = Degree+1;
		else						jj = Degree;
	ConstPointer( float ) xPtrs[ RestrictionStencil::RestrictionStencil::Size ];
	for( int yy=0 ; yy<RestrictionStencil::RestrictionStencil::Size ; yy++ )
		if( (yy+startY>=0 && yy+startY<minor) || periodicType==SPHERICAL_PERIODIC )	xPtrs[yy] = SocketedStreamingSolver< Channels , SyncType >::GetXRow( yy+startY , c );
		else													xPtrs[yy] = NullPointer< float >( );
		for( int i=0 ; i<highSize ; i++ )
		{
			int ii;
			if( periodicType==NO_PERIODIC )
				if		( i+highStart<Degree)				ii = i+highStart;
				else if	( i+highStart>highMajor-1-Degree)	ii = 2*Degree+1+(i+highStart-(highMajor-1));
				else if ((i+highStart-Degree)&1)			ii = Degree+1;
				else										ii = Degree;
			else
				if ( (i+highStart-Degree)&1 )	ii = Degree+1;
				else							ii = Degree;

			double value=0;

			int startX = RestrictionStencil::RestrictionStencil::Start(i+highStart) - _start;
			for( int yy=0 ; yy<RestrictionStencil::RestrictionStencil::Size ; yy++ )
			{
				if( !xPtrs[yy] )	continue;
				double tValue=0;

				for( int xx=0 ; xx<RestrictionStencil::RestrictionStencil::Size ; xx++ )
					if( (xx+startX+_start>=0 && xx+startX+_start<major) || periodicType!=NO_PERIODIC )
						tValue += xPtrs[yy][startX+xx]*majorRestrictionStencil.caseTable[ii].values[xx];
				value += tValue*minorRestrictionStencil.caseTable[jj].values[yy];
			}
			highX[i] += float( value * laplacianScale );
		}
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SolveRestriction(void)
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// ...and finish up the trailing parent
	if( parent ) parent->SolveRestriction();
	if( !inCore ) // Why only for out-of-core?
	{
		volatile MultiStreamIOClient *b = (MultiStreamIOClient*) _B;
		volatile MultiStreamIOClient *x = (MultiStreamIOClient*) _X;
		while( ( b && b->server ) || ( x && x->server ) ) SleepThisThread( 0 );
	}
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SolveProlongation(void)
{
	// Run to completion...
	while( IterateProlongation() ){;}
	// ...and finish up the trailing child
	if( pChild ) pChild->SolveProlongation();
	if( !inCore ) // Why only for out-of-core?
	{
		volatile MultiStreamIOClient *b = (MultiStreamIOClient*) _B;
		volatile MultiStreamIOClient *x = (MultiStreamIOClient*) _X;
		while( ( b && b->server ) || ( x && x->server ) ) SleepThisThread( 0 );
	}
}

///////////////
// LabelData //
///////////////
template<class LabelType,int Channels>
inline bool LabelData<LabelType,Channels>::operator == (const LabelData& ld) const
{
	bool ret = true;
	for( int c=0 ; c<Channels ; c++ ) ret &= (ld.l[c]==l[c]);
	return ret;
}
template<> inline bool LabelData<unsigned char		,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<unsigned char		,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<unsigned char		,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template<> inline bool LabelData<uint16_t			,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<uint16_t			,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<uint16_t			,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template<> inline bool LabelData<unsigned int		,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<unsigned int		,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<unsigned int		,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template<> inline bool LabelData<unsigned long long	,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<unsigned long long	,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<unsigned long long	,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template< class LabelType , int Channels >
LabelData< LabelType , Channels > LabelData< LabelType , Channels >::UnknownLabel( void )
{
	LabelData ld;
#if 1 // Unknown pixels are marked white
	for( int c=0 ; c<Channels ; c++ ) ld.l[c] = LabelType( -1 );
#else // Unknown pixels are marked black
	for( int c=0 ; c<Channels ; c++ ) ld.l[c] = LabelType( 0 );
#endif
	return ld;
}
template< class LabelType , int Channels >
LabelData< LabelType , Channels >::operator size_t () const
{
	size_t hash = 0;
	for( int c=0 ; c<Channels ; c++ ) hash ^= l[c];
	return hash;
}
template< class LabelType , int Channels >
bool LabelData< LabelType , Channels >::operator ()( const LabelData& ld1 , const LabelData& ld2 ) const
{
	for( int c=0 ; c<Channels ; c++ ) if( ld1.l[c]<ld2.l[c] ) return true;
	return false;
}
template< class LabelType , int Channels >
bool LabelData< LabelType , Channels >::operator < ( const LabelData& ld ) const
{
	for( int c=0 ; c<Channels ; c++ ) if( l[c]<ld.l[c] ) return true;
	return false;
}
/////////////////////////////////
// SocketedStreamingDivergence //
/////////////////////////////////
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SocketedStreamingDivergence( void )
{
	unknownType = UNKNOWN_BLACK;
	for( int i=0 ; i<ISize ; i++ ) labels[i] = NullPointer< LabelData< LabelType , LabelChannels > >( );
	for( int i=0 ; i<Degree ; i++ ) localPAccum[i] = NullPointer< __m128 >( );
	for( int i=0 ; i<ISize*PixelChannels ; i++ ) pixels[i] = NullPointer< __m128 >( );
	for( int i=0 ; i<ISize*PixelChannels ; i++ ) lowPixels[i] = NullPointer< __m128 >( );
	for( int i=0 ; i<DXSize*PixelChannels ; i++ ) dx[i] = NullPointer< __m128 >( );
	for( int i=0 ; i<DYSize*PixelChannels ; i++ ) dy[i] = NullPointer< __m128 >( );

	parent = NULL;

	for( int c=0 ; c<PixelChannels ; c++ ) average[c] = 0;
	syncBuffer = NullPointer< ImageData< PixelChannels , LabelChannels , SyncType , LabelType > >( );
	localDMajorAccum = NullPointer< __m128 >( );
	localDMinorAccum = NullPointer< __m128 >( );

	lapTemplates = AlignedAllocPointer< TemplateSSE >( 3 * (2*Degree+1) , ALIGNMENT );
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::~SocketedStreamingDivergence( void )
{
	FreePointer( syncBuffer );
	for( int i=0 ; i<ISize ; i++ ) FreePointer( labels[i]  );
	for( int i=0 ; i<Degree ; i++ )
		if( localPAccum[i] )
		{
			localPAccum[i] -= _padSize128;
			AlignedFreePointer( localPAccum[i] );
		}
	for( int i=0 ; i< ISize*PixelChannels ; i++ ) AlignedFreePointer( pixels[i] );
	for( int i=0 ; i< ISize*PixelChannels ; i++ ) AlignedFreePointer( lowPixels[i] );
	for( int i=0 ; i<DXSize*PixelChannels ; i++ ) AlignedFreePointer( dx[i] );
	for( int i=0 ; i<DYSize*PixelChannels ; i++ ) AlignedFreePointer( dy[i] );

	if( localDMajorAccum )
	{
		localDMajorAccum -= _padSize128;
		AlignedFreePointer( localDMajorAccum );
	}
	if( localDMinorAccum )
	{
		localDMinorAccum -= _padSize128;
		AlignedFreePointer( localDMinorAccum );
	}
	AlignedFreePointer( lapTemplates );
}

template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::Initialize( StreamingGrid* lowPxls , StreamingGrid *pxls , StreamingGrid* lbls ,
																						double iWeight , bool lump , double gWeight , double gScale ,
																					    int start , int end , int major , int minor , int iters ,
																					    DataStream* leftStream,Socket* syncSockets,DataStream* rightStream,
																					    bool memoryMappedFile , int periodicType , MultiStreamIOServer* server ,
																						const std::vector< std::pair< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > > >* gradientAverageMap
																					  )
{
	_separateLaplacianComputation = ( lowPxls!=NULL ) || ( lbls!=NULL ) || ( gradientAverageMap!=NULL && gradientAverageMap->size()!=0 );
#if FIX_DERAMP
	LabelData< LabelType , LabelChannels > unknown = LabelData< LabelType , LabelChannels >::UnknownLabel( );
	if( gradientAverageMap ) for( int i=0 ; i<gradientAverageMap->size() ; i++ )
		if( !( (*gradientAverageMap)[i].first==unknown && unknownType==UNKNOWN_BLACK ) ) _gradientAverageMap[ (*gradientAverageMap)[i].first ] = (*gradientAverageMap)[i].second;
#else // !FIX_DERAMP
	if( gradientAverageMap ) for( int i=0 ; i<gradientAverageMap->size() ; i++ )
		_gradientAverageMap[ (*gradientAverageMap)[i].first ] = (*gradientAverageMap)[i].second;
#endif // FIX_DERAMP

	_iWeight = iWeight , _lump = lump , _gWeight = gWeight , _gScale = gScale;
	if( syncSockets && !SetSocketedStreamData( major , minor , start , end , leftStream , syncSockets[0] , rightStream , periodicType ) )	exit(0);
	else if(           !SetSocketedStreamData( major , minor , start , end , leftStream , _INVALID_SOCKET_ , rightStream , periodicType ) )	exit(0);

	this->_periodicType = periodicType;
	this->major = major , this->minor = minor;
	lowPixelStream = lowPxls , pixelStream = pxls , labelStream = lbls;

	if( parent )
		if( syncSockets )
			parent->Initialize( iWeight , lump , gWeight , start , end , major , minor , iters ,
			leftStream , syncSockets+1 , rightStream ,
			memoryMappedFile , periodicType , server
			);
		else
			parent->Initialize( iWeight , lump , gWeight , start , end , major , minor , iters ,
			leftStream , NULL , rightStream ,
			memoryMappedFile , periodicType , server
			);
	for( int i=0 ; i<ISize ; i++ ) labels[i] = AllocPointer< LabelData< LabelType , LabelChannels > >( _paddedSize );
	for( int i=0 ; i<ISize*PixelChannels ; i++ ) pixels[i] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
	for( int i=0 ; i<ISize*PixelChannels ; i++ ) lowPixels[i] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
	for( int i=0 ; i<Degree ; i++ )
	{
		localPAccum[i] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
		localPAccum[i] += _padSize128;
	}

	for( int i=0 ; i<DXSize*PixelChannels ; i++ ) dx[i] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
	for( int i=0 ; i<DYSize*PixelChannels ; i++ ) dy[i] = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );

	syncBuffer = AllocPointer< ImageData< PixelChannels , LabelChannels , SyncType , LabelType > >( _size * Degree );
	localDMajorAccum = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
	localDMinorAccum = AlignedAllocPointer< __m128 >( _paddedSize128 , ALIGNMENT );
	localDMajorAccum += _padSize128;
	localDMinorAccum += _padSize128;


	MatrixStencil lStencil , dStencil , stencil;
	DotProductStencil dotMajor , d2DotMajor , dotMinor , d2DotMinor;

	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil( major , dotMajor , 0 , 0 );
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil( major , d2DotMajor , 1 , 1 , false );
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil( minor , dotMinor , 0 , 0 );
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil( minor , d2DotMinor , 1 , 1 , false );

	if( lump )
	{
		for( int i=0 ; i<=2*Degree ; i++ )
		{
			double majorSum = 0 , minorSum = 0;
			for( int j=0 ; j<=2*Degree ; j++ )
			{
				majorSum += dotMajor.caseTable[i].values[j] , dotMajor.caseTable[i].values[j] = 0;
				minorSum += dotMinor.caseTable[i].values[j] , dotMinor.caseTable[i].values[j] = 0;
			}
			dotMajor.caseTable[i].values[Degree] = majorSum;
			dotMinor.caseTable[i].values[Degree] = minorSum;
		}
	}

	for( int i=0 ; i<=2*Degree ; i++ ) for( int j=0 ; j<=2*Degree ; j++ )
		for( int l=0 ; l<=2*Degree ; l++ ) for( int m=0 ; m<=2*Degree ; m++ )
		{
			lStencil.caseTable[i][j].values[l][m] = d2DotMajor.caseTable[i].values[l] * dotMinor.caseTable[j].values[m] + dotMajor.caseTable[i].values[l] * d2DotMinor.caseTable[j].values[m];
			dStencil.caseTable[i][j].values[l][m] = dotMajor.caseTable[i].values[l] * dotMinor.caseTable[j].values[m];
		}
	if( _separateLaplacianComputation )
		for( int i=0 ; i<=2*Degree ; i++ ) for( int j=0 ; j<=2*Degree ; j++ )
			for( int k=0 ; k<=2*Degree ; k++ ) for( int l=0 ; l<=2*Degree ; l++ )
				stencil.caseTable[i][j].values[k][l] = dStencil.caseTable[i][j].values[k][l]*iWeight;
	else
		for( int i=0 ; i<=2*Degree ; i++ ) for( int j=0 ; j<=2*Degree ; j++ )
			for( int k=0 ; k<=2*Degree ; k++ ) for( int l=0 ; l<=2*Degree ; l++ )
				stencil.caseTable[i][j].values[k][l] = dStencil.caseTable[i][j].values[k][l]*iWeight + lStencil.caseTable[i][j].values[k][l]*gScale*gWeight;

	ALIGN( float scratch[4] , 16 );
	for( int i=0 ; i<=2*Degree ; i++ )	// Iterate over the minor index in the mask
		for( int j=0 ; j<3 ; j++ )
			for( int k=0 ; k<=2*Degree ; k++ )
			{
				int jj;

				if(j==0)	jj=0;
				else		jj=Degree;
				scratch[0] = float( stencil.caseTable[i][jj].values[k][2] );
				scratch[1] = float( stencil.caseTable[i][jj].values[k][3] );
				scratch[2] = float( stencil.caseTable[i][jj].values[k][4] );
				if(j!=2)	scratch[3] = float( stencil.caseTable[i][Degree].values[k][1] );
				else		scratch[3] = float( 0 );
				lapTemplates[3*i+j].matrixValues[0][k] = _mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l] = float( stencil.caseTable[i][jj].values[k][l+1] );
				lapTemplates[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l] = float( stencil.caseTable[i][jj].values[k][l] );
				lapTemplates[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0] = float( stencil.caseTable[i][Degree].values[k][3] );
				else		scratch[0] = float( 0 );
				scratch[1] = float( stencil.caseTable[i][jj].values[k][0] );
				scratch[2] = float( stencil.caseTable[i][jj].values[k][1] );
				scratch[3] = float( stencil.caseTable[i][jj].values[k][2] );
				lapTemplates[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
			}

}

template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::InitRestriction( void )
{
	index = -Degree-1;
	if( pixelStream->rowSize()!=_size*PixelChannels*sizeof(PixelType) )
		fprintf( stderr , "[ERROR] Pixel width failure: %d != %d\n" , pixelStream->rowSize() , (int)(_size*PixelChannels*sizeof(PixelType)) ) , exit(0);
	if( lowPixelStream )
	{
		if( minor!=lowPixelStream->rows() )
			fprintf( stderr , "[ERROR] Low pixel height failure: %d != %d\n" , minor , lowPixelStream->rows() ) , exit(0);
		if( lowPixelStream->rowSize()!=_size*PixelChannels*sizeof(PixelType) )
			fprintf( stderr , "[ERROR] Low pixel width failure: %d != %d\n" , lowPixelStream->rowSize() , (int)(_size*PixelChannels*sizeof(PixelType)) ) , exit(0);
	}
	if( labelStream )
	{
		if( minor!=labelStream->rows() )
			fprintf( stderr , "[ERROR] Label height failure: %d != %d\n" , minor , labelStream->rows() ) , exit(0);
		if( labelStream->rowSize()!=_size*LabelChannels*sizeof(LabelType) )
			fprintf( stderr , "[ERROR] Label width failure: %d != %d\n" , labelStream->rowSize() , (int)(_size*LabelChannels*sizeof(LabelType)) ) , exit(0);
	}

	pixelStream->reset( true , 1 );
	if( lowPixelStream ) lowPixelStream->reset( true , 1 );
	if( labelStream ) labelStream->reset( true , 1 );
	if( _separateLaplacianComputation )
	{
		FiniteElements1D< float , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( major , dotMajorStencil , 0 , 0 );
		FiniteElements1D< float , Type , Degree >::DotProduct< Type , Degree >::DotProductStencil( minor , dotMinorStencil , 0 , 0 );
		FiniteElements1D< float , DERIVATIVE(Type) , Degree-1 >::DotProduct< Type , Degree >::DotProductStencil( FiniteElements1D< float , DERIVATIVE(Type) , Degree-1 >::Dimension( FiniteElements1D< float , Type , Degree >::DomainSize(major) ) , dDotMajorStencil , 0 , 1 );
		FiniteElements1D< float , DERIVATIVE(Type) , Degree-1 >::DotProduct< Type , Degree >::DotProductStencil( FiniteElements1D< float , DERIVATIVE(Type) , Degree-1 >::Dimension( FiniteElements1D< float , Type , Degree >::DomainSize(minor) ) , dDotMinorStencil , 0 , 1 );
		FiniteElements2D< float , Type , Degree> ::DivergenceStencil( major , minor , divergenceStencil );
		if( _lump )
		{
			for( int i=0 ; i<=2*Degree ; i++ )
			{
				double majorSum = 0 , minorSum = 0;
				for( int j=0 ; j<=2*Degree ; j++ )
				{
					majorSum += dotMajorStencil.caseTable[i].values[j] , dotMajorStencil.caseTable[i].values[j] = 0;
					minorSum += dotMinorStencil.caseTable[i].values[j] , dotMinorStencil.caseTable[i].values[j] = 0;
				}
				dotMajorStencil.caseTable[i].values[Degree] = (float)majorSum;
				dotMinorStencil.caseTable[i].values[Degree] = (float)minorSum;
			}
		}

		for( int i=0 ; i<2*Degree+1 ; i++ ) for( int k=0 ; k<2*Degree ; k++ )
		{
			dDotMajorStencil.caseTable[i].values[k] *= float( _gScale * _gWeight );
			dDotMinorStencil.caseTable[i].values[k] *= float( _gScale * _gWeight );
		}
		for( int i=0 ; i<2*Degree+1 ; i++ ) for( int j=0 ; j<2*Degree+2 ; j++ )
			for( int k=0 ; k<2*Degree ; k++ ) for( int l=0 ; l<2*Degree+1 ; l++ )
			{
				divergenceStencil.caseTable[i][j].values1[k][l] *= float( _gScale * _gWeight );
				divergenceStencil.caseTable[i][j].values2[l][k] *= float( _gScale * _gWeight );
			}
	}

	if( parent ) parent->InitRestriction();
}

template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::GetPixelRow(int row,int channel)
{
	return ( Pointer( float ) )( pixels[MyModIndex( row*PixelChannels + channel , ISize*PixelChannels )] ) + _padSize;
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::GetLowPixelRow( int row , int channel )
{
	return ( Pointer( float ) )( lowPixels[MyModIndex( row*PixelChannels + channel , ISize*PixelChannels )] ) + _padSize;
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( LabelData< LabelType , LabelChannels > ) SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::GetLabelRow( int row )
{
	return labels[ MyModIndex( row , ISize ) ] + _padSize;
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::GetDXRow(int row,int channel)
{
	return ( Pointer( float ) )( dx[MyModIndex(row*PixelChannels + channel , DXSize*PixelChannels)] + _padSize128 );
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::GetDYRow(int row,int channel)
{
	return ( Pointer( float ) )( dy[MyModIndex(row*PixelChannels + channel , DYSize*PixelChannels )] + _padSize128 );
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SyncImageHead( int idx , bool read )
{
	if( idx<0 || idx>=Degree )
	{
		fprintf( stderr , "Synchronizing divergence in non-head row: %d\n" , idx );
		return;
	}
	int r = -idx-1;
	if( syncXSocket!=_INVALID_SOCKET_ )
	{
		if( read )
		{
			ReceiveOnSocket( syncXSocket, syncBuffer , sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * _size , "SocketedStreamingDivergence::SyncImageHead" );
			// Copy the data back from the stream
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( r , c );
				for( int x=0 ; x<_size ; x++ ) pixelRow[x] = float( syncBuffer[x].pixel[c] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( r , c );
					for( int x=0 ; x<_size ; x++ ) lowPixelRow[x] = float( syncBuffer[x].lowPixel[c] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( r );
				for( int c=0 ; c<LabelChannels ; c++ ) for( int x=0 ; x<_size ; x++ ) labelRow[x].l[c] = syncBuffer[x].label[c];
			}
		}
		else
		{
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( idx , c );
				for( int x=0 ; x<_size ; x++ ) syncBuffer[x].pixel[c] = SyncType( pixelRow[x] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( idx , c );
					for( int x=0 ; x<_size ; x++ ) syncBuffer[x].lowPixel[c] = SyncType( lowPixelRow[x] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( idx );
				for( int c=0 ; c<LabelChannels ; c++ ) for( int x=0 ; x<_size ; x++ ) syncBuffer[x].label[c] = labelRow[x].l[c]; 
			}
			SendOnSocket( syncXSocket , ( ConstPointer( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) )syncBuffer, sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * _size , "SocketedStreamingDivergence::SyncImageHead" );
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SyncImageTail( int idx , bool read )
{
	if( idx<minor-Degree || idx>=minor )
	{
		fprintf( stderr , "Synchronizing divergence in non-tail row: %d\n" , idx );
		return;
	}
	int r = 2 * minor - 1 - idx;
	if( syncXSocket!=_INVALID_SOCKET_ )
	{
		if( read )
		{
			ReceiveOnSocket( syncXSocket, syncBuffer , sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * _size , "SocketedStreamingDivergence::SyncImageTail" );
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( r , c );
				for( int x=0 ; x<_size ; x++ ) pixelRow[x] = float( syncBuffer[x].pixel[c] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( r , c );
					for( int x=0 ; x<_size ; x++ ) lowPixelRow[x] = float( syncBuffer[x].lowPixel[c] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( r );
				for( int c=0 ; c<LabelChannels ; c++ ) for( int x=0 ; x<_size ; x++ ) labelRow[x].l[c] = syncBuffer[x].label[c];
			}
		}
		else
		{
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( idx , c );
				for( int x=0 ; x<_size ; x++ ) syncBuffer[x].pixel[c] = SyncType( pixelRow[x] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( idx , c );
					for( int x=0 ; x<_size ; x++ ) syncBuffer[x].lowPixel[c] = SyncType( lowPixelRow[x] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( idx );
				for( int c=0 ; c<LabelChannels ; c++ ) for( int x=0 ; x<_size ; x++ ) syncBuffer[x].label[c] = labelRow[x].l[c];
			}
			SendOnSocket( syncXSocket , ( ConstPointer( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) )syncBuffer , sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * _size , "SocketedStreamingDivergence::SyncImageTail" );
		}
	}
}

template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SyncImageLeft( int j , bool read )
{
	if( leftStream )
	{
		int xSize = 2*_padSize;
		Pointer( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) left = AllocPointer< ImageData< PixelChannels , LabelChannels , SyncType , LabelType > >( xSize );
		if( read )
		{
			if ( !leftStream->read( ( Pointer( byte ) )left , sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * xSize ) ) exit(0);
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j , c ) - _padSize;
				for( int i=0 ; i<_padSize ; i++ ) pixelRow[i] = float( left[i].pixel[c] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j , c ) - _padSize;
					for( int i=0 ; i<_padSize ; i++ ) lowPixelRow[i] = float( left[i].lowPixel[c] );
				}
			}
			if( labelStream ) 
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow(j) - _padSize;
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) labelRow[i].l[c] = left[i].label[c];
			}
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j+1 , c ) - _padSize;
				for( int i=0 ; i<_padSize ; i++ ) pixelRow[i] = float( left[_padSize+i].pixel[c] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j+1 , c ) - _padSize;
					for( int i=0 ; i<_padSize ; i++ ) lowPixelRow[i] = float( left[_padSize+i].lowPixel[c] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow(j+1) - _padSize;
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) labelRow[i].l[c] = left[_padSize+i].label[c];
			}
		}
		else
		{
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j , c );
				for( int i=0 ; i<_padSize ; i++ ) left[i].pixel[c] = SyncType( pixelRow[i] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j , c );
					for( int i=0 ; i<_padSize ; i++ ) left[i].lowPixel[c] = SyncType( lowPixelRow[i] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( j );
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) left[i].label[c] = labelRow[i].l[c];
			}
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j+1 , c );
				for( int i=0 ; i<_padSize ; i++ ) left[_padSize+i].pixel[c] = SyncType( pixelRow[i] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j+1 , c );
					for( int i=0 ; i<_padSize ; i++ ) left[_padSize+i].lowPixel[c] = SyncType( lowPixelRow[i] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( j+1 );
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) left[_padSize+i].label[c] = labelRow[i].l[c];
			}
			if ( !leftStream->write( ( Pointer( byte ) )left , sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * xSize ) )	exit(0);
		}
		FreePointer( left );
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SyncImageRight( int j , bool read )
{
	if( rightStream )
	{
		int xSize = 2*_padSize;
		Pointer( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) right = AllocPointer< ImageData< PixelChannels , LabelChannels , SyncType , LabelType > >( xSize );
		if( read )
		{
			if ( !rightStream->read( ( Pointer( byte ) )right , sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * xSize ) ) exit(0);
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j , c ) + _size;
				for( int i=0 ; i<_padSize ; i++ ) pixelRow[i] = float( right[i].pixel[c] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j , c ) + _size;
					for( int i=0 ; i<_padSize ; i++ ) lowPixelRow[i] = float( right[i].lowPixel[c] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( j ) + _size;
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) labelRow[i].l[c] = right[i].label[c];
			}
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j+1 , c ) + _size;
				for( int i=0 ; i<_padSize ; i++ ) pixelRow[i] = float( right[_padSize+i].pixel[c] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j+1 , c ) + _size;
					for( int i=0 ; i<_padSize ; i++ ) lowPixelRow[i] = float( right[_padSize+i].lowPixel[c] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( j+1 ) + _size;
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) labelRow[i].l[c] = right[_padSize+i].label[c];
			}
		}
		else
		{
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j , c ) + _size - _padSize;
				for( int i=0 ; i<_padSize ; i++ ) right[i].pixel[c] = SyncType( pixelRow[i] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j , c ) + _size - _padSize;
					for( int i=0 ; i<_padSize ; i++ ) right[i].lowPixel[c] = SyncType( lowPixelRow[i] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( j ) + _size - _padSize;
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) right[i].label[c] = labelRow[i].l[c];
			}
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( j+1 , c ) + _size - _padSize;
				for( int i=0 ; i<_padSize ; i++ ) right[_padSize+i].pixel[c] = SyncType( pixelRow[i] );
			}
			if( lowPixelStream )
			{
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( j+1 , c ) + _size - _padSize;
					for( int i=0 ; i<_padSize ; i++ ) right[_padSize+i].lowPixel[c] = SyncType( lowPixelRow[i] );
				}
			}
			if( labelStream )
			{
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( j+1 ) + _size - _padSize;
				for( int c=0 ; c<LabelChannels ; c++ ) for( int i=0 ; i<_padSize ; i++ ) right[_padSize+i].label[c] = labelRow[i].l[c];
			}
			if ( !rightStream->write( ( Pointer( byte ) )right , sizeof( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) * xSize ) ) exit(0);
		}
		FreePointer( right );
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::UnSetRestriction( void )
{
	if( parent ) parent->UnSetRestriction();
}

template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SetRestriction( void )
{
	if( parent ) parent->SetRestriction();
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SetRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	if( _periodicType==NO_PERIODIC )
	{
		if( sRestrict<0     ) sRestrict = 0;
		if( eRestrict>major ) eRestrict = major;
	}
	int myStart = sRestrict-parent->start() , myEnd = eRestrict-parent->start();
	memset( lB+myStart, 0 , sizeof( float ) * (myEnd-myStart) );

	if( (idx>Degree && idx<minor-Degree) || _periodicType==SPHERICAL_PERIODIC )
	{
		if( _separateLaplacianComputation )
		{
			if( _gWeight ) _AddInteriorDerivativeRestriction( lB , c , idx , sRestrict , eRestrict );
			if( _iWeight ) _AddInteriorStencilRestriction( lB , c , idx , sRestrict , eRestrict );
		}
		else _AddInteriorStencilRestriction( lB , c , idx , sRestrict , eRestrict );
	}
	else
	{
		if( _separateLaplacianComputation )
		{
			if( _gWeight ) _AddDerivativeRestriction( lB , c , idx , sRestrict , eRestrict );
			if( _iWeight ) _AddStencilRestriction( lB , c , idx , sRestrict , eRestrict );
		}
		else _AddStencilRestriction( lB , c , idx , sRestrict , eRestrict );
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_AddInteriorDerivativeRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	int  off = FiniteElements1D< float , Type , Degree >::DotProduct< Type , Degree>::Helper::StartOffset();
	int dOff = FiniteElements1D< float , DERIVATIVE(Type) , Degree-1 >::DotProduct< Type , Degree >::Helper::StartOffset();
	int  I = idx+ off;
	int dI = idx+dOff;
	Pointer( __m128 ) majorPtrs[2*Degree+1];
	Pointer( __m128 ) minorPtrs[2*Degree  ];

	int highStart = sRestrict-_start-Degree;
	int highEnd = eRestrict-_start+Degree;
	int highStart128 = (highStart-(RealPerWord-1))/RealPerWord;
	int highEnd128   = (highEnd  +(RealPerWord-1))/RealPerWord;
	int myStart = sRestrict-parent->start();
	int myEnd   = eRestrict-parent->start();

	{
		for( int xx=0 ; xx<=2*Degree ; xx++ ) majorPtrs[xx] = ( Pointer( __m128 ) )GetDXRow(  I+xx , c );
		for( int xx=0 ; xx< 2*Degree ; xx++ ) minorPtrs[xx] = ( Pointer( __m128 ) )GetDYRow( dI+xx , c );
		ALIGN( float scratch[4] , 16 );
		__m128 dot[Degree+1],dDot[Degree];
		for( int d=0 ; d<=Degree ; d++ )
		{
			for( int j=0 ; j<4 ; j++ ) scratch[j] = dotMinorStencil.caseTable[Degree].values[d];
			dot[d]=_mm_load_ps(scratch);
		}
		for( int d=0 ; d<Degree ; d++ )
		{
			for( int j=0 ; j<4 ; j++ ) scratch[j] = dDotMinorStencil.caseTable[Degree].values[d];
			dDot[d]=_mm_load_ps(scratch);
		}
		for(int i=highStart128;i<highEnd128;i++)
		{
			localDMajorAccum[i]=                                                  _mm_mul_ps(dot[Degree],majorPtrs[Degree][i]);
			localDMajorAccum[i]=_mm_add_ps(localDMajorAccum[i],_mm_mul_ps( dot[0],_mm_add_ps(majorPtrs[0][i],majorPtrs[2*Degree  ][i])));
			localDMajorAccum[i]=_mm_add_ps(localDMajorAccum[i],_mm_mul_ps( dot[1],_mm_add_ps(majorPtrs[1][i],majorPtrs[2*Degree-1][i])));
		}
		for(int i=highStart128;i<highEnd128;i++)
		{
			localDMinorAccum[i]=                               _mm_mul_ps(dDot[0],_mm_sub_ps(minorPtrs[0][i],minorPtrs[2*Degree-1][i]));
			localDMinorAccum[i]=_mm_add_ps(localDMinorAccum[i],_mm_mul_ps(dDot[1],_mm_sub_ps(minorPtrs[1][i],minorPtrs[2*Degree-2][i])));
		}
		Pointer( float ) localDX = ( ( Pointer( float ) )localDMinorAccum)+FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();
		Pointer( float ) localDY = ( ( Pointer( float ) )localDMajorAccum)+FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();

		if( _periodicType!=SPHERICAL_PERIODIC )
		{
			if( myStart+_start<Degree )
			{
				for( int j=myStart ; j<Degree-_start ; j++ )
				{
					double temp=0;
					const float* xValues =  dotMajorStencil.caseTable[j].values;
					const float* yValues = dDotMajorStencil.caseTable[j].values;

					// Partial w.r.t minor index
					int J = off+j;
					for(int dj=0;dj<=2*Degree;dj++)
					{
						if( J+dj+_start<0 ||  J+dj+_start>=major)	continue;
						temp+=localDX[j+dj] * xValues[dj];
					}
					// Partial w.r.t major index
					int DJ = dOff + j;
					for(int dj=0;dj<2*Degree;dj++)
					{
						if(DJ+dj+_start<0 || DJ+dj+_start>=major-1)	continue;
						temp+=localDY[j+dj] * yValues[dj];
					}
					lB[j] += float( temp );
				}
				myStart = Degree-_start;
			}
			if( myEnd+_start>major-Degree )
			{
				for( int j=major-Degree-_start ; j<myEnd ; j++ )
				{
					double temp=0;
					const float* xValues =  dotMajorStencil.caseTable[2*Degree+(j+_start-(major-1))].values;
					const float* yValues = dDotMajorStencil.caseTable[2*Degree+(j+_start-(major-1))].values;

					// Partial w.r.t minor index
					int J = off+j;
					for(int dj=0;dj<=2*Degree;dj++)
					{
						if( J+dj+_start<0 ||  J+dj+_start>=major)	continue;
						temp+=localDX[j+dj] * xValues[dj];
					}
					// Partial w.r.t major index
					int DJ = dOff + j;
					for(int dj=0;dj<2*Degree;dj++)
					{
						if(DJ+dj+_start<0 || DJ+dj+_start>=major-1)	continue;
						temp+=localDY[j+dj] * yValues[dj];
					}
					lB[j] += float( temp );
				}
				myEnd = major-Degree-_start;
			}
		}
		const float* xValues =  dotMajorStencil.caseTable[Degree].values;
		const float* yValues = dDotMajorStencil.caseTable[Degree].values;
		for( int j=myStart ; j<myEnd ; j++ )
		{
			double temp=0;
			// Partial w.r.t minor index
			temp = localDX[j+Degree]*xValues[Degree];
			for(int dj=0;dj<Degree;dj++) temp += (localDX[j+dj]+localDX[j+2*Degree-dj])*xValues[dj];

			// Partial w.r.t major index
			for(int dj=0;dj<Degree;dj++) temp += (localDY[j+dj]-localDY[j+2*Degree-1-dj])*yValues[dj];

			lB[j] += float( temp );
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_AddDerivativeRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	int ii;
	if( _periodicType==SPHERICAL_PERIODIC )	ii=Degree;
	else
		if     ( idx<Degree )        ii=idx;
		else if( idx>=minor-Degree ) ii=2*Degree+(idx-(minor-1));
		else                         ii=Degree;
	int dI=idx+FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
	int I =idx+FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();
	int myStart = sRestrict-parent->start();
	int myEnd = eRestrict-parent->start();

	for( int j=myStart ; j<myEnd ; j++ )
	{
		double temp=0;
		int jj;
		if( _periodicType!=NO_PERIODIC )	jj=Degree;
		else
			if(j+_start<Degree)				jj=j+_start;
			else if(j+_start>=major-Degree)	jj=2*Degree+(j+_start-(major-1));
			else							jj=Degree;
		int dJ=j+FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
		int J =j+FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

		// Partial w.r.t minor index
		for( int di=0 ; di<2*Degree ; di++ )
		{
			if( (dI+di<0 || dI+di>=minor-1) && _periodicType!=SPHERICAL_PERIODIC ) continue;
			Pointer( float ) localD = GetDYRow( dI+di , c ) + J;
			for( int dj=0 ; dj<=2*Degree ; dj++ )
			{
				if( (J+dj+_start<0 || J+dj+_start>=major) && _periodicType==NO_PERIODIC ) continue;
				temp += localD[dj]*divergenceStencil.caseTable[jj][ii].values2[dj][di];
			}
		}

		// Partial w.r.t major index
		for( int di=0 ; di<=2*Degree ; di++ )
		{
			if( (I+di<0 || I+di>=minor) && _periodicType!=SPHERICAL_PERIODIC )	continue;
			Pointer( float ) localD = GetDXRow( I+di , c ) + dJ;
			for(int dj=0;dj<2*Degree;dj++)
			{
				if( (dJ+dj+_start<0 || dJ+dj+_start>=major-1) && _periodicType==NO_PERIODIC )	continue;
				temp += localD[dj]*divergenceStencil.caseTable[jj][ii].values1[dj][di];
			}
		}
		lB[j] += float( temp );
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_AddInteriorStencilRestriction( Pointer( float ) lB,int c,int idx,int sRestrict,int eRestrict)
{
	int jj=Degree*3;
	{
		ConstPointer( float ) localPPtrs[2*Degree+1];
		if( lowPixelStream ) for( int y=0 ; y<2*Degree+1 ; y++ ) localPPtrs[y] = GetLowPixelRow( idx-Degree+y , c );
		else                 for( int y=0 ; y<2*Degree+1 ; y++ ) localPPtrs[y] =    GetPixelRow( idx-Degree+y , c );
		{
			ConstPointer( __m128 ) pPtrs[] =
			{
				(ConstPointer( __m128 ) )localPPtrs[0],
				(ConstPointer( __m128 ) )localPPtrs[1],
				(ConstPointer( __m128 ) )localPPtrs[3],
				(ConstPointer( __m128 ) )localPPtrs[4]
			};

			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for(int i=-WordPerDegree+((sRestrict-_start)>>2);i<((eRestrict-_start)>>2)+WordPerDegree;i++)
			{
				localPAccum[0][i]=_mm_add_ps(pPtrs[0][i],pPtrs[3][i]);
				localPAccum[1][i]=_mm_add_ps(pPtrs[1][i],pPtrs[2][i]);
			}
		}


		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) pPtrs[] = { localPAccum[0] , localPAccum[1] , ( ConstPointer( __m128 ) )localPPtrs[2] };
		int s,e;
		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
			{
				dotSum=0;
				lB[0] += GetInteriorLaplacianValue0( lapTemplates[jj].matrixValues[0] , pPtrs , dotSum , 0 );
				s=4;
			}
			else
			{
				SetInteriorDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
			else								e=eRestrict-_start;

			for(int i=s;i<e;i+=4)				lB[i] += GetInteriorLaplacianValue0(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-4;
			if( eRestrict==major && _periodicType==NO_PERIODIC )	lB[i] += GetInteriorLaplacianValue0(lapTemplates[jj+2].matrixValues[0],pPtrs,dotSum,i);
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
			{
				dotSum=0;
				lB[1] += GetInteriorLaplacianValue1(lapTemplates[jj].matrixValues[1],pPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetInteriorDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)			lB[i] += GetInteriorLaplacianValue1(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-3;
			if( eRestrict==major && _periodicType==NO_PERIODIC )lB[i] += GetInteriorLaplacianValue1(lapTemplates[jj+2].matrixValues[1],pPtrs,dotSum,i);
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[2],pPtrs,0,dSum);
			else							SetInteriorDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			if( eRestrict==major && _periodicType==NO_PERIODIC)	e=eRestrict-_start-8;
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i] += GetInteriorLaplacianValue2(mValues,pPtrs,dotSum,i);
			if( eRestrict==major && _periodicType==NO_PERIODIC )
			{
				int i=eRestrict-_start-6;
				lB[i] += GetInteriorLaplacianValue2(lapTemplates[jj+2].matrixValues[2],pPtrs,dotSum,i);
				i+=4;
				lB[i] += dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[3],pPtrs,0,dSum);
			else							SetInteriorDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-8;
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i] += GetInteriorLaplacianValue3(mValues,pPtrs,dotSum,i);
			if( eRestrict==major && _periodicType==NO_PERIODIC )
			{
				int i=eRestrict-_start-5;
				lB[i] += GetInteriorLaplacianValue3(lapTemplates[jj+2].matrixValues[3],pPtrs,dotSum,i);
				i+=4;
				lB[i] += dotSum;
			}
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_AddStencilRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	int jj;
	if		(idx<Degree)			jj=idx;
	else if	(idx>minor-1-Degree)	jj=2*Degree+(idx-(minor-1));
	else							jj=Degree;

	if( _periodicType==NO_PERIODIC )
	{
		if( sRestrict<0 )		sRestrict=0;
		if( eRestrict>major )	eRestrict=major;
	}
	if( jj==Degree || _periodicType==SPHERICAL_PERIODIC ) return _AddInteriorStencilRestriction( lB , c , idx , sRestrict , eRestrict );
	jj*=3;

	{
		ConstPointer( float ) localPPtrs[2*Degree+1];
		for(int y=0;y<2*Degree+1;y++)
			if(idx-Degree+y>=0 && idx-Degree+y<minor)
				if( lowPixelStream ) localPPtrs[y] = GetLowPixelRow( idx-Degree+y , c );
				else                 localPPtrs[y] =    GetPixelRow( idx-Degree+y , c );
			else
				localPPtrs[y] = NullPointer< float >( );
		float dotSum;
		__m128 dSum;
		ALIGN( float scratch[4] , 16 );
		ConstPointer( __m128 ) pPtrs[] =
		{
			( ConstPointer( __m128 ) )localPPtrs[0],
			( ConstPointer( __m128 ) )localPPtrs[1],
			( ConstPointer( __m128 ) )localPPtrs[2],
			( ConstPointer( __m128 ) )localPPtrs[3],
			( ConstPointer( __m128 ) )localPPtrs[4]
		};
		int s,e;
		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2],
				lapTemplates[jj+1].matrixValues[0][3],
				lapTemplates[jj+1].matrixValues[0][4]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
			{
				dotSum = 0;
				lB[0] += GetLaplacianValue0( lapTemplates[jj].matrixValues[0] , pPtrs , dotSum , 0 );
				s=4;
			}
			else
			{
				SetDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)				lB[i] += GetLaplacianValue0(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-4;
			if( eRestrict==major && _periodicType==NO_PERIODIC )	lB[i] += GetLaplacianValue0(lapTemplates[jj+2].matrixValues[0],pPtrs,dotSum,i);
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2],
				lapTemplates[jj+1].matrixValues[1][3],
				lapTemplates[jj+1].matrixValues[1][4]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
			{
				dotSum=0;
				lB[1] += GetLaplacianValue1(lapTemplates[jj].matrixValues[1],pPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
			if( eRestrict==major &&  _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)				lB[i] += GetLaplacianValue1(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-3;
			if( eRestrict==major && _periodicType==NO_PERIODIC )	lB[i] += GetLaplacianValue1(lapTemplates[jj+2].matrixValues[1],pPtrs,dotSum,i);
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2],
				lapTemplates[jj+1].matrixValues[2][3],
				lapTemplates[jj+1].matrixValues[2][4]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[2],pPtrs,0,dSum);
			else							SetDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-8;
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i] += GetLaplacianValue2(mValues,pPtrs,dotSum,i);
			if( eRestrict==major && _periodicType==NO_PERIODIC )
			{
				int i=eRestrict-_start-6;
				lB[i] += GetLaplacianValue2(lapTemplates[jj+2].matrixValues[2],pPtrs,dotSum,i);
				i+=4;
				lB[i] += dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2],
				lapTemplates[jj+1].matrixValues[3][3],
				lapTemplates[jj+1].matrixValues[3][4]
			};
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[3],pPtrs,0,dSum);
			else							SetDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-8;
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i] += GetLaplacianValue3(mValues,pPtrs,dotSum,i);
			if( eRestrict==major && _periodicType==NO_PERIODIC )
			{
				int i=eRestrict-_start-5;
				lB[i] += GetLaplacianValue3(lapTemplates[jj+2].matrixValues[3],pPtrs,dotSum,i);
				i+=4;
				lB[i] += dotSum;
			}
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_setPartials( int y )
{
	if( y&1 )
	{
		SyncImageLeft( y , false );
		if( _separateLaplacianComputation )
		{
			_setPartialsX( y , -_padSize , _size-_padSize );
			_setPartialsY( y , -_padSize , _size-_padSize );
		}
		SyncImageRight( y , true );
		if( _separateLaplacianComputation )
		{
			_setPartialsX( y , _size-_padSize-1 , _size+_padSize );
			_setPartialsY( y , _size-_padSize , _size+_padSize );
		}
	}
	else
	{
		SyncImageRight( y , false );
		if( _separateLaplacianComputation )
		{
			_setPartialsX( y , _padSize , _size+_padSize );
			_setPartialsY( y , _padSize , _size+_padSize );
		}
		SyncImageLeft( y , true );
		if( _separateLaplacianComputation )
		{
			_setPartialsX( y , -_padSize , _padSize+1 );
			_setPartialsY( y , -_padSize , _padSize );
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_setPartialsX( int y , int start , int end )
{
	Pointer( float ) pixelRow[PixelChannels];
	Pointer( float ) dx[PixelChannels];
	Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( y );
	for( int c=0 ; c<PixelChannels ; c++ )
	{
		pixelRow[c] = GetPixelRow( y , c );
		dx[c] = GetDXRow( y , c );
	}
	LabelData< LabelType , LabelChannels > unknown = LabelData< LabelType , LabelChannels >::UnknownLabel( );
	for( int x=start ; x<end-1 ; x++ )
	{
		bool useGradient = (labelStream==NULL) || (labelRow[x]==labelRow[x+1]);
		if( labelStream )
		{
			bool b0 = ( labelRow[x  ] == unknown );
			bool b1 = ( labelRow[x+1] == unknown );
			if     ( unknownType==UNKNOWN_BLACK    ) useGradient |=  b0 ||  b1;
			else if( unknownType==UNKNOWN_HARMONIC ) useGradient &= !b0 && !b1;
		}
		if( useGradient ) for( int c=0 ; c<PixelChannels ; c++ ) dx[c][x] = pixelRow[c][x+1]-pixelRow[c][x];
		else              for( int c=0 ; c<PixelChannels ; c++ ) dx[c][x] = 0;
		if( useGradient )
		{
			typename std::map< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > >::iterator iter = _gradientAverageMap.find( labelRow[x] );
			if( iter!=_gradientAverageMap.end() ) for( int c=0 ; c<PixelChannels ; c++ ) dx[c][x] -= float( iter->second.dx[c] );
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_setPartialsY( int y , int start , int end )
{
	Pointer( float ) dy[PixelChannels];
	Pointer( float ) oldPixelRow[PixelChannels];
	Pointer( float ) newPixelRow[PixelChannels];
	Pointer( LabelData< LabelType , LabelChannels > ) oldLabelRow = GetLabelRow( y-1 );
	Pointer( LabelData< LabelType , LabelChannels > ) newLabelRow = GetLabelRow( y   );
	for( int c=0 ; c<PixelChannels ; c++ )
	{
		oldPixelRow[c] = GetPixelRow( y-1 , c );
		newPixelRow[c] = GetPixelRow( y   , c );
		dy[c] = GetDYRow( y-1 , c );
	}

	LabelData< LabelType , LabelChannels > unknown = LabelData< LabelType , LabelChannels >::UnknownLabel( );
	for( int x=start ; x<end ; x++ )
	{
		bool useGradient = (labelStream==NULL) || (newLabelRow[x]==oldLabelRow[x]);
		if( labelStream )
		{
			bool b0 = ( oldLabelRow[x]==unknown );
			bool b1 = ( newLabelRow[x]==unknown );
			if     ( unknownType==UNKNOWN_BLACK    ) useGradient |=  b0 ||  b1;
			else if( unknownType==UNKNOWN_HARMONIC ) useGradient &= !b0 && !b1;
		}
		if( useGradient ) for( int c=0 ; c<PixelChannels ; c++ ) dy[c][x] = newPixelRow[c][x]-oldPixelRow[c][x];
		else              for( int c=0 ; c<PixelChannels ; c++ ) dy[c][x] = 0;
		if( useGradient )
		{
			typename std::map< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > >::iterator iter = _gradientAverageMap.find( newLabelRow[x] );
			if( iter!=_gradientAverageMap.end() ) for( int c=0 ; c<PixelChannels ; c++ ) dy[c][x] -= float( iter->second.dy[c] );
		}
	}
#if 0
	if( _periodicType==SPHERICAL_PERIODIC && y==0 )     if( (x>=major/4 && x<major/2) || (x>=3*major/4 && x<major ) ) useGradient = false;
	if( _periodicType==SPHERICAL_PERIODIC && y==minor ) if( (x>=major/2 && x<major) ) useGradient = false;
#endif
}

template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_setPartialsX( int y )
{
	float *dx[PixelChannels];
	float *pixelRow[PixelChannels];
	LabelData<LabelType,LabelChannels>* labelRow = GetLabelRow(y);
	for(int c=0;c<PixelChannels;c++)
	{
		pixelRow[c] = GetPixelRow(y,c);
		dx[c] = GetDXRow(y,c);
	}
	for(int x=-_padSize;x<_size+_padSize-1;x++)
	{
		bool useGradient = (labelStream==NULL) || (labelRow[x]==labelRow[x+1]);
		if( useGradient ) for( int c=0 ; c<PixelChannels ; c++ ) dx[c][x] = pixelRow[c][x+1]-pixelRow[c][x];
		else              for( int c=0 ; c<PixelChannels ; c++ ) dx[c][x] = 0;
		if( useGradient )
		{
			typename std::map< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > >::iterator iter = _gradientAverageMap.find( labelRow[x] );
			if( iter!=_gradientAverageMap.end() ) for( int c=0 ; c<PixelChannels ; c++ ) dx[c][x] -= iter->second.dx[c];
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::_setPartialsY( int y )
{
	float *dy[PixelChannels];
	float *oldPixelRow[PixelChannels],*newPixelRow[PixelChannels];
	LabelData< LabelType , LabelChannels > *oldLabelRow = GetLabelRow(y-1);
	LabelData< LabelType , LabelChannels > *newLabelRow = GetLabelRow(y  );
	for(int c=0;c<PixelChannels;c++)
	{
		oldPixelRow[c] = GetPixelRow(y-1,c);
		newPixelRow[c] = GetPixelRow(y  ,c);
		dy[c] = GetDYRow( y-1 , c );
	}
	for(int x=-_padSize;x<_size+_padSize;x++)
	{
		bool useGradient = (labelStream==NULL) || (newLabelRow[x]==oldLabelRow[x]);
		if(useGradient)		for(int c=0;c<PixelChannels;c++)	dy[c][x] = newPixelRow[c][x]-oldPixelRow[c][x];
		else				for(int c=0;c<PixelChannels;c++)	dy[c][x] = 0;
		if( useGradient )
		{
			typename std::map< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > >::iterator iter = _gradientAverageMap.find( newLabelRow[x] );
			if( iter!=_gradientAverageMap.end() ) for( int c=0 ; c<PixelChannels ; c++ ) dy[c][x] -= iter->second.dy[c];
		}
	}
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
bool SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::IterateRestriction( void )
{
	if( index>=minor ) return false;
	int idx = index+Degree+1;
	// Load the next row of pixel/label values
	if( idx<minor )
	{
		{
			if( labelStream )
			{
				Pointer( LabelType ) labels = ( Pointer( LabelType ) )(*labelStream)[idx];
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow(idx);
				for( int x=0 ; x<_size ; x++ ) for( int c=0 ; c<LabelChannels ; c++ ) labelRow[x].l[c] = labels[LabelChannels*x+c];
				labelStream->advance();
				// Set beyond the bounds, this should get over-written on synchronization (I hope)
				for( int x=-_padSize ; x<0 ; x++ ) labelRow[x] = labelRow[0];
				for( int x=_size ; x<_size+_padSize ; x++ ) labelRow[x] = labelRow[_size-1];
			}
			if( lowPixelStream )
			{
				Pointer( PixelType ) lowPixels = ( Pointer( PixelType ) )(*lowPixelStream)[idx];
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) lowPixelRow = GetLowPixelRow( idx , c );
					for( int x=0 ; x<_size ; x++ ) lowPixelRow[x] = float( lowPixels[ x*PixelChannels + c ] );
					// Set beyond the bounds, this should get over-written on synchronization (I hope)
					for( int x=-_padSize ; x<0 ; x++ ) lowPixelRow[x] = lowPixelRow[0];
					for( int x=_size ; x<_size+_padSize ; x++ ) lowPixelRow[x] = lowPixelRow[_size-1];
				}
				lowPixelStream->advance();
			}
			{
				Pointer( PixelType ) pixels = ( Pointer( PixelType ) )(*pixelStream)[idx];
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					Pointer( float ) pixelRow = GetPixelRow( idx , c );
					for( int x=0 ; x<_size ; x++ ) pixelRow[x] = float( pixels[ x*PixelChannels + c ] );
					// Set beyond the bounds, this should get over-written on synchronization (I hope)
					for( int x=-_padSize ; x<0 ; x++ ) pixelRow[x] = pixelRow[0];
					for( int x=_size ; x<_size+_padSize ; x++ ) pixelRow[x] = pixelRow[_size-1	];
					// Note: We should probably set the rest of the pixel rows to something consistent
					// so that we don't get unpredictable behavior when we're using padding on a spherical image
					// (not that it is reasonable to pad in the spherical case)
					if( labelStream && unknownType==UNKNOWN_BLACK )
					{
						LabelData< LabelType , LabelChannels > unknown = LabelData< LabelType , LabelChannels >::UnknownLabel( );
						Pointer( LabelData< LabelType , LabelChannels > ) labelRow = GetLabelRow( idx );
						for( int x=0 ; x<_size ; x++ )
						{
							if( labelRow[x]==unknown ) pixelRow[x] = 0.f;
							average[c] += pixelRow[x];
						}
					}
					else for( int x=0 ; x<_size ; x++ ) average[c] += pixelRow[x];
				}
				pixelStream->advance();
			}
		}
	}
	if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1 ) for( int d=0 ; d<Degree ; d++ ) SyncImageHead( d , false ) , SyncImageHead( d , true );
	if( _periodicType==SPHERICAL_PERIODIC && idx== minor-1 ) for( int d=minor-Degree ; d<minor ; d++ ) SyncImageTail( d , false ) , SyncImageTail( d , true );
	// If this is the first row, perform the preliminary left/right synchronization
	if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1 )
	{
		SyncImageLeft( -Degree , false ) , SyncImageRight( -Degree , false );
		SyncImageLeft( -Degree , true  ) , SyncImageRight( -Degree , true  );
	}
	if( _periodicType!=SPHERICAL_PERIODIC && idx==0 )
	{
		SyncImageLeft( 0 , false ) , SyncImageRight( 0 , false );
		SyncImageLeft( 0 , true  ) , SyncImageRight( 0 , true  );
	}
	idx = index+Degree;
	if( idx>=0 )
	{
		if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1 ) for( int y=-Degree ; y<=0 ; y++ ) _setPartials( y );
		if( _periodicType==SPHERICAL_PERIODIC && idx==minor-1 ) for( int y=minor; y<minor+Degree ; y++ ) _setPartials( y );
		_setPartials( idx );
		if( parent ) parent->IterateRestriction();
	}
	index++;
	return true;
}
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >::SolveRestriction( void )
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// ...and finish up the trailing parent
	if( parent ) parent->SolveRestriction();
}
//////////////////////////////////////
// SocketedStreamingGradientAverage //
//////////////////////////////////////
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void SetGradientAverageMap( StreamingGrid* pixelStream , StreamingGrid* labelStream , std::vector< std::pair< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > > >& averageMap , bool spherical , bool showProgress )
{
	typedef std::map< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > > AverageMap;
	std::map< LabelData< LabelType , LabelChannels > , int > unmap;
	AverageMap _averageMap;
	if( !pixelStream ) fprintf( stderr , "[ERROR] Pixel stream undefined\n" ) , exit(0);
	int height = pixelStream->rows() , width = pixelStream->rowSize() / ( PixelChannels * sizeof( PixelType ) );
	if( labelStream && height!=labelStream->rows() ) 
		fprintf( stderr , "[ERROR] Pixel and Label heights differ: %d != %d\n" , height , labelStream->rows() ) , exit( 0 );
	if( labelStream && width*LabelChannels*sizeof( LabelType )!=labelStream->rowSize() )
		fprintf( stderr , "[ERROR] Pixel and Label widths differ: %d != %zd\n" , width , labelStream->rowSize() / ( LabelChannels * sizeof( LabelType ) ) ) , exit( 0 );

	Pointer( float ) pixelRows[2];
	Pointer( LabelData< LabelType , LabelChannels > ) labelRows[2];
	pixelRows[0] = AllocPointer< float >( width*PixelChannels );
	pixelRows[1] = AllocPointer< float >( width*PixelChannels );

	if( labelStream )
	{
		labelRows[0] = AllocPointer< LabelData< LabelType , LabelChannels > >( width );
		labelRows[1] = AllocPointer< LabelData< LabelType , LabelChannels > >( width );
	}

	for( int y=0 ; y<=height ; y++ )
	{
		// Accumulate values per row and then combine. (There shouldn't be overflow with doubles, but just in case...)
		AverageMap _subMap;
		if( showProgress ) printf( "[%.1f%%]        \r" , float(y)/height*100 );
		// Read in the next row of pixels
		if( y<height )
		{
			if( labelStream )
			{
				Pointer( LabelType ) labels = ( Pointer( LabelType ) )(*labelStream)[y];
				Pointer( LabelData< LabelType , LabelChannels > ) labelRow = labelRows[y%2];
				for( int x=0 ; x<width ; x++ ) for( int c=0 ; c<LabelChannels ; c++ ) labelRow[x].l[c] = labels[LabelChannels*x+c];
				if( spherical && ( y==0 || y==height-1 ) ) for( int x=0 ; x<width ; x++ ) unmap[ labelRow[x] ] = 0;
				labelStream->advance();
			}
			{
				Pointer( PixelType ) pixels = ( Pointer( PixelType ) )(*pixelStream)[y];
				Pointer( float ) pixelRow = pixelRows[y%2];
				for( int x=0 ; x<width ; x++ ) for( int c=0 ; c<PixelChannels ; c++ ) pixelRow[x*PixelChannels+c] = float( pixels[x*PixelChannels+c] );
				pixelStream->advance();
			}
		}

		// Set the x-partials
		if( y<height )
		{
			Pointer( float ) pixelRow = pixelRows[y%2];
			Pointer( LabelData< LabelType , LabelChannels > ) labelRow = labelRows[y%2];
			for( int x=0 ; x<width-1 ; x++ )
			{
				if( (!labelStream) || ( labelRow[x]==labelRow[x+1] ) )
				{
					GradientAverage< PixelChannels > dx;
					dx.dxCount = 1;
					for( int c=0 ; c<PixelChannels ; c++ ) dx.dx[c] = pixelRow[PixelChannels*(x+1)+c] - pixelRow[PixelChannels*x+c];
					LabelData< LabelType , LabelChannels > label = labelStream!=NULL ? labelRow[x] : LabelData< LabelType , LabelChannels >::UnknownLabel();
					typename AverageMap::iterator iter = _subMap.find( label );
					if( iter!=_subMap.end() ) iter->second += dx;
					else _subMap[ label ] = dx;
				}
			}
		}
		// Set the y-partials
		if( y>0 )
		{
			Pointer( float ) oldPixelRow = pixelRows[(y+1)%2];
			Pointer( float ) newPixelRow = pixelRows[(y  )%2];
			Pointer( LabelData< LabelType , LabelChannels > ) oldLabelRow = labelRows[(y+1)%2];
			Pointer( LabelData< LabelType , LabelChannels > ) newLabelRow = labelRows[(y  )%2];
			for( int x=0 ; x<width ; x++ )
			{
				if( (!labelStream) || ( oldLabelRow[x]==newLabelRow[x] ) )
				{
					GradientAverage< PixelChannels > dy;
					dy.dyCount = 1;
					for( int c=0 ; c<PixelChannels ; c++ ) dy.dy[c] = newPixelRow[PixelChannels*x+c] - oldPixelRow[PixelChannels*x+c];
					LabelData< LabelType , LabelChannels > label = labelStream!=NULL ? oldLabelRow[x] : LabelData< LabelType , LabelChannels >::UnknownLabel();
					typename AverageMap::iterator iter = _subMap.find( label );
					if( iter!=_subMap.end() ) iter->second += dy;
					else _subMap[ label ] = dy;
				}
			}
		}
		// Now add the values computed from the row to the overall average
		for( typename AverageMap::iterator subIter=_subMap.begin() ; subIter!=_subMap.end() ; subIter++ )
		{
			typename AverageMap::iterator iter = _averageMap.find( subIter->first );
			if( iter==_averageMap.end() ) _averageMap[ subIter->first ] = subIter->second;
			else iter->second += subIter->second;
		}
	}
	if( showProgress ) printf( "\n" );

	FreePointer( pixelRows[0] );
	FreePointer( pixelRows[1] );
	if( labelStream )
	{
		FreePointer( labelRows[0] );
		FreePointer( labelRows[1] );
	}
	int idx = 0;
	averageMap.resize( _averageMap.size() );
	for( typename AverageMap::iterator iter=_averageMap.begin() ; iter!=_averageMap.end() ; iter++ )
	{
		averageMap[idx].first = iter->first;
		if( unmap.find( iter->first )==unmap.end() ) averageMap[idx].second = iter->second;
		else                                         averageMap[idx].second = GradientAverage< PixelChannels >();
		idx++;
	}
}
