#include "Util/XPlatform.h"
#include "Util/MultiStreamIO.h"
#include "Util/Socket.h"
#include <omp.h>

#define DEBUG_SOCKETS 0
#define LOW_MEMORY_HACK 0
#define USE_PARALLEL_CG 1

template< bool AddAverage >
void Multiply( const SparseMatrix< double >& A , ConstPointer( double ) in , Pointer( double ) out )
{
#if USE_PARALLEL_CG
#pragma omp parallel for
#endif // USE_PARALLEL_CG
	for( int i=0 ; i<A.groups ; i++ )
	{
		out[i] = 0.;
		for( int j=0 ; j<A.groupSizes[i] ; j++ ) out[i] += in[ A[i][j].N ] * A[i][j].Value;
	}
	if( AddAverage )
	{
		double average = 0.;
#if USE_PARALLEL_CG
#pragma omp parallel for reduction( + : average )
#endif // USE_PARALLEL_CG
		for( int i=0 ; i<A.groups ; i++ ) average += in[i];
		average /= A.groups;
#if USE_PARALLEL_CG
#pragma omp parallel for
#endif // USE_PARALLEL_CG
		for( int i=0 ; i<A.groups ; i++ ) out[i] += average;
	}
}
template< bool AddAverage >
int MySolveConjugateGradient( const SparseMatrix< double >& A , ConstPointer( double ) b , Pointer( double ) x , int iters )
{
	double eps=1e-16;
	int dim = A.groups;
	Pointer( double )  r = AllocPointer< double >( dim );
	Pointer( double )  d = AllocPointer< double >( dim );
	Pointer( double ) Ad = AllocPointer< double >( dim );

	double delta_new = 0 , delta_0;

	////////////////////////
	// d = r = b - M * x
	// \delta_new = ||r||^2
	Multiply< AddAverage >( A , x , d );
#if USE_PARALLEL_CG
#pragma omp parallel for
#endif // USE_PARALLEL_CG
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = b[i] - d[i];
#if USE_PARALLEL_CG
#pragma omp parallel for reduction( + : delta_new )
#endif // USE_PARALLEL_CG
	for( int i=0 ; i<dim ; i++ ) delta_new += r[i] * r[i];
	////////////////////////

	delta_0 = delta_new;
	if( delta_new<eps )
	{
		fprintf( stderr , "[WARNING] MySolveConjugateGradient: Initial residual too low: %g < %g\n" , delta_new , eps );
		FreePointer( r );
		FreePointer( d );
		FreePointer( Ad );
		return 0;
	}
	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		////////////////////////////////////
		// \alpha = ||r||^2 / (d^t * M * d)
		Multiply< AddAverage >( A , ( ConstPointer( double ) )d , Ad );
		double dDotMd = 0;
#if USE_PARALLEL_CG
#pragma omp parallel for reduction( + : dDotMd )
#endif // USE_PARALLEL_CG
		for( int i=0 ; i<dim ; i++ ) dDotMd += Ad[i] * d[i];
		double alpha = double( delta_new / dDotMd );
		double delta_old = delta_new;
		////////////////////////////////////

		delta_new = 0;

		if( (ii%50)==(50-1) )
		{
			////////////////////////////////
			// x = x + d * alpha
			// r = b - M * ( x + d * alpha )
			//   = r - M * d * alpha
			// \delta_new = ||r||^2
#if USE_PARALLEL_CG
#pragma omp parallel for
#endif // USE_PARALLEL_CG
			for( int i=0 ; i<dim ; i++ ) x[i] += d[i] * alpha;
			memset( r , 0 , sizeof(double) * dim );
			Multiply< AddAverage >( A , x , r );
#if USE_PARALLEL_CG
#pragma omp parallel for reduction( + : delta_new )
#endif // USE_PARALLEL_CG
			for( int i=0 ; i<dim ; i++ ) r[i] = b[i] - r[i] , delta_new += r[i] * r[i];
			////////////////////////////////
		}
		else
		{
			////////////////////////////////
			// x = x + d * alpha
			// r = r - M * d * alpha
			// \delta_new = ||r||^2
#if USE_PARALLEL_CG
#pragma omp parallel for reduction( + : delta_new )
#endif // USE_PARALLEL_CG
			for( int i=0 ; i<dim ; i++ ) r[i] -= Ad[i] * alpha , delta_new += r[i] * r[i] , x[i] += d[i] * alpha;
			////////////////////////////////
		}
		////////////////////////////////
		// beta = ||r||^2 / ||r_old||^2
		// d = r + d * beta
		double beta = delta_new / delta_old;
#if USE_PARALLEL_CG
#pragma omp parallel for
#endif // USE_PARALLEL_CG
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + d[i] * beta;
		////////////////////////////////
	}
	FreePointer( r );
	FreePointer( d );
	FreePointer( Ad );
	return ii;
}

void SetPaddedSize( int width , int height , int minRes , int& paddedWidth , int& paddedHeight , int& widthBlockSize, int& heightBlockSize )
{
	// Compute the size of the domain
	int domainW = FiniteElements1D< float , Type , Degree >::DomainSize( width  );
	int domainH = FiniteElements1D< float , Type , Degree >::DomainSize( height );
	if( domainW<minRes ) fprintf( stderr , "[ERROR] Width is smaller than minimum-resolution: %d<%d\n" , domainW , minRes ) , exit( 0 );
	if( domainH<minRes ) fprintf( stderr , "[ERROR] Height is smaller than minimum-resolution: %d<%d\n" , domainH , minRes ) , exit( 0 );
	// Want domainW <= _paddedWidth = a * 2^b for the smallest a \in [ minRes , 2*minRes )
	// b is as small as possible the property that domainW <= (2*minRes-1)^2^b
	// -- by assumption (2*minRes-1)/2 < domainW
	// => the smallest possible value of b is 0.
#if 1
	long long aW = (2*minRes-1) , aH = (2*minRes-1);
	int bW=0 , bH=0;
	// First find the smallest value of b
	while( domainW>(aW<<bW) ) bW++ ; while( domainH>(aH<<bH) ) bH++;
	widthBlockSize = 1<<bW , heightBlockSize = 1<<bH;
	// Now find the smallest value of a
	while( domainW<=( (aW-1)<<bW ) ) aW-- ; while( domainH<=( (aH-1)<<bH ) ) aH--;
	long long _paddedWidth = (aW<<bW) , _paddedHeight = (aH<<bH);
#else
	long long _paddedWidth = minRes*2-1 , _paddedHeight = minRes*2-1;
	widthBlockSize = heightBlockSize = 1;
	while( domainW>_paddedWidth  ) _paddedWidth  *= 2 , widthBlockSize  *= 2;
	while( domainH>_paddedHeight ) _paddedHeight *= 2 , heightBlockSize *= 2;
	// Now we have b such that (2*minRes-1)*2^{b-1} < domainW <= (2*minRes-1)*2^b

	while( _paddedWidth - widthBlockSize>=domainW ) _paddedWidth  -= widthBlockSize;
	while( _paddedHeight-heightBlockSize>=domainH ) _paddedHeight -= heightBlockSize;
#endif
	paddedWidth  = (int)( FiniteElements1D<float,Type,Degree>::Dimension( (int)( _paddedWidth  ) ) );
	paddedHeight = (int)( FiniteElements1D<float,Type,Degree>::Dimension( (int)( _paddedHeight ) ) );
}
void SetPaddedSize( int width , int height , int minRes , int& paddedWidth , int& paddedHeight )
{
	int widthBlockSize, heightBlockSize;
	SetPaddedSize( width , height , minRes , paddedWidth , paddedHeight , widthBlockSize , heightBlockSize );
}
bool IsDownSamplable( int width , int height , int iters , int& newWidth , int& newHeight )
{
	if( (width&3) || (height&3) || (width<20) )	return false;
	newWidth  = width>>1;
	newHeight = height>>1;
	return true;
}
bool IsDownSamplable(int width,int height,int& newWidth,int& newHeight)
{
	if( (width&3) || (height&3) || (width<20) )	return false;
	newWidth  = width>>1;
	newHeight = height>>1;
	return true;
}
bool IsDownSamplable(int dCount,ClientSocket* cSockets,int sockCount)
{
	for(int i=0;i<sockCount;i++)
		if ( (cSockets[i].clientData.start>>dCount)&3 || ((cSockets[i].clientData.end-cSockets[i].clientData.start)>>dCount)<20)	return false;
	return true;
}
template<class Real>
void SetRectangularLaplacianMatrix( DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
								    SparseMatrix< Real > & lap , int width , int height , double iWeight , double gWeight )
{
	lap.Resize( width*height );
	int ii , jj;
	int xStart , xEnd , yStart , yEnd;
	for( int y=0 ; y<height ; y++ )
	{
		if		(y<Degree)			jj = y							, yStart = -y		, yEnd = Degree;
		else if	(y>height-1-Degree)	jj = 2*Degree+(y-(height-1))	, yStart = -Degree	, yEnd = height-1-y;
		else						jj = Degree						, yStart = -Degree	, yEnd = Degree;

		for( int x=0 ; x<width ; x++ )
		{
			if		(x<Degree)			ii = x						, xStart = -x		, xEnd = Degree;
			else if	(x>width-1-Degree)	ii = 2*Degree+(x-(width-1))	, xStart = -Degree	, xEnd = width-1-x;
			else						ii = Degree					, xStart = -Degree	, xEnd = Degree;

			int idx = x+y*width;
			lap.SetGroupSize( idx , (yEnd-yStart+1)*(xEnd-xStart+1) );
			int _i = 0;
			for( int yy=yStart ; yy<=yEnd ; yy++ )
				for( int xx=xStart ; xx<=xEnd ; xx++ )
				{
					lap[idx][_i  ].N = (x+xx)+(y+yy)*width;
					lap[idx][_i++].Value = 
						( dotMajor.caseTable[ii].values[xx+Degree]*d2DotMinor.caseTable[jj].values[yy+Degree] + d2DotMajor.caseTable[ii].values[xx+Degree]*dotMinor.caseTable[jj].values[yy+Degree] ) * gWeight +
						( dotMajor.caseTable[ii].values[xx+Degree]*dotMinor.caseTable[jj].values[yy+Degree] ) * iWeight;
				}
		}
	}
}
template<class Real>
void SetCylindricalLaplacianMatrix( DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
								    SparseMatrix< Real > & lap , int width , int height , double iWeight , double gWeight )
{
	lap.Resize( width*height );

	int jj;
	int yStart , yEnd;
	for( int y=0 ; y<height ; y++ )
	{
		if		(y<Degree)			jj = y							, yStart = -y		, yEnd = Degree;
		else if	(y>height-1-Degree)	jj = 2*Degree+(y-(height-1))	, yStart = -Degree	, yEnd = height-1-y;
		else						jj = Degree						, yStart = -Degree	, yEnd = Degree;

		for( int x=0 ; x<width ; x++ )
		{
			int idx = x+y*width;
			lap.SetGroupSize( idx , (yEnd-yStart+1)*(2*Degree+1) );
			int _i = 0;
			for(int yy=yStart;yy<=yEnd;yy++)
				for( int i=-Degree ; i<=Degree ; i++ )
				{
					int xx=(x+i+width)%width;
					lap[idx][_i  ].N=xx+(y+yy)*width;
					lap[idx][_i++].Value=
						( dotMajor.caseTable[Degree].values[i+Degree]*d2DotMinor.caseTable[jj].values[yy+Degree] + d2DotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[jj].values[yy+Degree] ) * gWeight +
						( dotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[jj].values[yy+Degree] ) * iWeight;
				}
		}
	}
}
template<class Real>
void SetSphericalLaplacianMatrix(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
								 SparseMatrix<Real>& lap , int width , int height , double iWeight , double gWeight )
{
	lap.Resize(width*height);

	for(int x=0;x<width;x++) for(int y=0;y<height;y++)
	{
		int idx1 = x + y*width;
		lap.SetGroupSize( idx1 , (2*Degree+1)*(2*Degree+1) );
		for(int i=-Degree;i<=Degree;i++) for(int j=-Degree;j<=Degree;j++)
		{
			int idx2 = (i+Degree)+(j+Degree)*(2*Degree+1);
			int xx=(x+i+width)%width;
			int yy= y+j;
			if(yy<0)
			{
				yy = -yy-1;
				if(xx<width/2)	xx = width/2-1-xx;
				else			xx = width/2+width-1-xx;
			}
			else if(yy>=height)
			{
				yy = height-1-(yy-height);
				xx = width-1-xx;
			}
			lap[idx1][idx2].N = xx+yy*width;
			lap[idx1][idx2].Value = 
				( dotMajor.caseTable[Degree].values[i+Degree]*d2DotMinor.caseTable[Degree].values[j+Degree] + d2DotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree] ) * gWeight +
				( dotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree] ) * iWeight;
		}
	}
}
template< class Real >
void UpSampleSphericalSolution( const std::vector< Real >&low , int width , int height , std::vector< Real >& high )
{
	int width2 , height2;
	FiniteElements1D< double , Type , Degree >::FullProlongationStencil wStencil , hStencil;
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( width  , wStencil , width2  );
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( height , hStencil , height2 );
	high.resize( width2 * height2 );
	for( int i=0 ; i<width ; i++ )
	{
		int iStart = FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( i );
		for( int j=0 ; j<height ; j++ )
		{
			int jStart = FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( j );
			for( int x=0 ; x<FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size ; x++ )
			{
				for( int y=0 ; y<FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size ; y++ )
				{
					int xx = ( iStart + x + width2 ) % width2;
					int yy= jStart+y;
					if( yy<0 )
					{
						yy = -yy-1;
						if( xx<width2/2 )	xx = width2/2          - 1 - xx;
						else				xx = width2/2 + width2 - 1 - xx;
					}
					else if( yy>=height2 )
					{
						yy = height2-1-(yy-height2);
						xx = width2-1-xx;
					}
					high[ yy*width2 + xx ] += Real( low[j*width+i] * wStencil.caseTable[Degree].values[x] * hStencil.caseTable[Degree].values[y] );
				}
			}
		}
	}
}
template< class Real >
void UpSampleCylindricalSolution( const std::vector< Real >&low , int width , int height , std::vector< Real >& high )
{
	int width2 , height2;
	FiniteElements1D< double , Type , Degree >::FullProlongationStencil wStencil , hStencil;
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( width  , wStencil , width2 );
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( height , hStencil , height2 );
	high.resize( width2 * height2 );
	for( int i=0 ; i<width ; i++ )
	{
		int iStart = FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( i );
		for( int j=0 ; j<height ; j++ )
		{
			int jj;
			if( j<Degree )				jj = j;
			else if( j>=height-Degree)	jj = 2*Degree+(j-(height-1));
			else						jj = Degree;
			int jStart = FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( j );
			for( int x=0 ; x<FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size ; x++ )
			{
				int xx = ( iStart + x + width2 ) % width2;
				for( int y=0 ; y<FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size ; y++ )
					if( jStart+y<0 || jStart+y>=height2 ) continue;
					else 
						high[ (jStart+y)*width2 + xx ] += Real( low[j*width+i] * wStencil.caseTable[Degree].values[x] * hStencil.caseTable[jj].values[y] );
			}
		}
	}
}
template< class Real >
void UpSampleRectangularSolution( const std::vector< Real >&low , int width ,int height , std::vector< Real >& high )
{
	int width2 , height2;
	FiniteElements1D< double , Type , Degree >::FullProlongationStencil wStencil , hStencil;
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( width  , wStencil , width2 );
	FiniteElements1D< double , Type , Degree >::ProlongationStencil( height , hStencil , height2 );
	high.resize( width2 * height2 );
	for( int i=0 ; i<width ; i++ )
	{
		int iStart = FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( i );
		int ii;
		if( i<Degree )				ii = i;
		else if( i>=width-Degree)	ii = 2*Degree+(i-(width-1));
		else						ii = Degree;
		for( int j=0 ; j<height ; j++ )
		{
			int jj;
			if( j<Degree )				jj = j;
			else if( j>=height-Degree)	jj = 2*Degree+(j-(height-1));
			else						jj = Degree;
			int jStart = FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( j );
			for( int x=0 ; x<FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size ; x++ )
				if( iStart+x<0 || iStart+x>=width2 ) continue;
				else
					for( int y=0 ; y<FiniteElements1D< double , Type , Degree >::FullProlongationStencil::ProlongationStencil::Size ; y++ )
						if( jStart+y<0 || jStart+y>=height2 ) continue;
						else 
							high[ (jStart+y)*width2 + (iStart+x) ] += Real( low[j*width+i] * wStencil.caseTable[ii].values[x] * hStencil.caseTable[jj].values[y] );
		}
	}
}
///////////////////////////
// SphericalSynchronizer //
///////////////////////////
SphericalSynchronizer::SphericalSynchronizer(void)
{
	_cCount = 0;
	_clientSockets = NULL;
}
void SphericalSynchronizer::init( int width , const int* widths , int cCount , int rPasses , int pPasses , int vCycles )
{
	_width=width;
	_cCount=cCount;
	_rPasses=rPasses;
	_pPasses=pPasses;
	_vCycles=vCycles;

	if( _clientSockets ) delete[] _clientSockets;

	_clientSockets = new ClientSocket[cCount];
	if(!_clientSockets)
	{
		fprintf(stderr,"Failed to allocate SphericalSynchronizer::_clientSockets\n");
		exit(0);
	}
	for( int i=0 ; i<cCount ; i++ )
	{
		_clientSockets[i].listener = GetListenSocket(_clientSockets[i].port);
		if(_clientSockets[i].listener == _INVALID_ACCEPTOR_SOCKET_ ) fprintf(stderr,"Failed to set server in SphericalSynchronizer\n") , exit(0);
		if( i ) _clientSockets[i].clientData.start = _clientSockets[i-1].clientData.end;
		else    _clientSockets[i].clientData.start = 0;
		_clientSockets[i].clientData.end = _clientSockets[i].clientData.start+widths[i];
	}
}
SphericalSynchronizer::~SphericalSynchronizer(void)
{
	delete[] _clientSockets;
	_clientSockets=NULL;
}
int SphericalSynchronizer::clients(void)		const	{ return _cCount; }
int SphericalSynchronizer::port(int cIndex)	const	{ return _clientSockets[cIndex].port; }
template< class DataType >
void SphericalSynchronizer::Run( void )
{
	for(int i=0;i<_cCount;i++)
	{
		_clientSockets[i].socket = AcceptSocket( _clientSockets[i].listener );
		if ( _clientSockets[i].socket == _INVALID_SOCKET_ ) exit(0);
	}

	Pointer( DataType ) buffer1;
	Pointer( DataType ) buffer2;
	buffer1 = AllocPointer< DataType >( _width );
	buffer2 = AllocPointer< DataType >( _width );
	for ( int v = 0 ; v < _vCycles ; v++ )
	{
		//////////////////////////////////////////
		// First synchronize on the restriction //
		//////////////////////////////////////////

		// Synchronize the head of the stream
		for(int k=0;k<_rPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				ReceiveOnSocket( _clientSockets[i].socket , buffer1 , sizeof(DataType)*w , "Failed to get restriction head from client" );
				for( int x=0 ; x<w ; x++ ) buffer2[x+_clientSockets[i].clientData.start] = buffer1[x];
			}
			for(int x=0;x<_width/2;x++)
			{
				int xx = _width/2-1-x;
				buffer1[xx] = buffer2[x];
			}
			for(int x=_width/2;x<_width;x++)
			{
				int xx = _width/2+_width-1-x;
				buffer1[xx] = buffer2[x];
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				for( int x=0 ; x<w ; x++ ) buffer2[x] = buffer1[x+_clientSockets[i].clientData.start];
				SendOnSocket( _clientSockets[i].socket , ( ConstPointer( DataType ) )buffer2 , sizeof(DataType)*w , "Failed to send restriction head to client" ) ;
			}
		}
		// Synchronize the tail of the stream
		for(int k=0;k<_rPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				ReceiveOnSocket( _clientSockets[i].socket , buffer1 , sizeof(DataType)*w , "Failed to get restriction tail from client" );
				for( int x=0 ; x<w ; x++ ) buffer2[x+_clientSockets[i].clientData.start] = buffer1[x];
			}
			for( int x=0 ; x<_width ; x++ )
			{
				int xx = _width-1-x;
				buffer1[xx] = buffer2[x];
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				for( int x=0 ; x<w ; x++ ) buffer2[x] = buffer1[x+_clientSockets[i].clientData.start];
				SendOnSocket( _clientSockets[i].socket , ( ConstPointer( DataType ) )buffer2 , sizeof(DataType)*w , "Failed to send restriction tail to client" );
			}
		}

		/////////////////////////////////////////
		// Now synchronize on the prolongation //
		/////////////////////////////////////////

		// Synchronize the head of the stream
		for(int k=0;k<_pPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				ReceiveOnSocket( _clientSockets[i].socket , buffer1 , sizeof(DataType)*w , "Failed to get prolongation head from client" );
				for( int x=0 ; x<w ; x++ ) buffer2[x+_clientSockets[i].clientData.start] = buffer1[x];
			}
			for(int x=0;x<_width/2;x++)
			{
				int xx = _width/2-1-x;
				buffer1[xx] = buffer2[x];
			}
			for(int x=_width/2;x<_width;x++)
			{
				int xx = _width/2+_width-1-x;
				buffer1[xx] = buffer2[x];
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				for( int x=0 ; x<w ; x++ ) buffer2[x] = buffer1[x+_clientSockets[i].clientData.start];
				SendOnSocket( _clientSockets[i].socket , ( ConstPointer( DataType ) )buffer2 , sizeof(DataType)*w , "Failed to send prolongation head to client" );
			}
		}
		// Synchronize the tail of the stream
		for(int k=0;k<_pPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				ReceiveOnSocket( _clientSockets[i].socket , buffer1 , sizeof(DataType)*w , "Failed to get prolongation tail from client" );
				for( int x=0 ; x<w ; x++ ) buffer2[x+_clientSockets[i].clientData.start] = buffer1[x];
			}
			for(int x=0;x<_width;x++)
			{
				int xx = _width-1-x;
				buffer1[xx] = buffer2[x];
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].clientData.end-_clientSockets[i].clientData.start;
				for( int x=0 ; x<w ; x++ ) buffer2[x] = buffer1[x+_clientSockets[i].clientData.start];
				SendOnSocket( _clientSockets[i].socket , ( ConstPointer( DataType ) )buffer2 , sizeof(DataType)*w , "Failed to send prolongation tail to client" );
			}
		}
	}
	FreePointer( buffer1 );
	FreePointer( buffer2 );
}
template< class DataType >
int SphericalSynchronizer::RunThread( void* vparams )
{
	SphericalSynchronizer* synchronizer = (SphericalSynchronizer*)vparams;
	synchronizer->Run< DataType >();
	return 0;
}
#if MISHA_CODE_CLEAN_UP
template< int PixelChannels , int LabelChannels , class SyncType , class LabelType >
static void SocketedMultigridServer< PixelChannels , LabelChannels , SyncType , LabelType >::Run
(
	char* prefix , int port , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles , int minBandSize ,
	int tileWidth , int tileHeight , const char* tileExt , 
	bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType ,
	double iWeight , bool lump , double gWeight , double gScale ,
	bool removeAverageGradient , int unknownType , bool showProgress ,
	bool noCG , bool shortSync
)
{
	GlobalData globalData;
	char address[512];

	if( minBandSize > inCoreRes ) inCoreRes = minBandSize;

	globalData.iters = iters;
	globalData.vCycles = vCycles;
	globalData.tileWidth = tileWidth;
	globalData.tileHeight = tileHeight;
	if( tileExt ) strcpy( globalData.tileExt , tileExt );
	else globalData.tileExt[0] = 0;
	globalData.gammaCorrection = gammaCorrection;
	globalData.quality = quality;
	globalData.lanes = lanes;
	globalData.verbose = verbose;
	globalData.periodicType = periodicType;
	globalData.iWeight = iWeight;
	globalData.lump = lump;
	globalData.gWeight = gWeight;
	globalData.gScale = gScale;
	globalData.unknownType = unknownType;
	globalData.showProgress = showProgress;
	globalData.shortSync = shortSync;

	GetHostAddress( address , prefix );

	int subClientCount = 0;

	// Create a listening SOCKET for connecting to server
	AcceptorSocket listenSocket = _INVALID_ACCEPTOR_SOCKET_;
	listenSocket = GetListenSocket( port );
	if ( listenSocket == _INVALID_ACCEPTOR_SOCKET_ ) return false;
	printf( "Server Address: %s:%d\n", address , port ) , fflush( stdout );

	// Geth the client information
	{
		ClientSocket* clientSockets = new ClientSocket[clientCount];
		int coreCount = 0;

		// Establish a connection to the clients
		for( int i=0 ; i<clientCount ; i++ )
		{
			clientSockets[i].socket  = AcceptSocket( listenSocket );
			if( globalData.verbose ) printf( "Connected to process: %d /%d     \r" , i+1 , clientCount );
		}
		if( globalData.verbose ) printf( "\n" ) , fflush( stdout );

		// Get the client info
		int width = 0 , height;
		for( int i=0 ; i<clientCount ; i++ )
		{
			ReceiveOnSocket( clientSockets[i].socket , GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) , "Failed to get client data from client" );
			if(!i) height = clientSockets[i].clientData.height;
			else if( height!=clientSockets[i].clientData.height ) fprintf( stderr , "Band heights differ: %d != %d\n" , height , clientSockets[i].clientData.height ) , exit(0);
			width += clientSockets[i].clientData.width;
			subClientCount += clientSockets[i].clientData.subClients;
		}
		// Sort the clients by index order
		qsort( clientSockets , clientCount , sizeof( ClientSocket ) , ClientSocket::Sort );

		// Store the clipping width and then compute the padded dimensions
		globalData.cWidth  = width;
		globalData.cHeight = height;

#if 0
#pragma message( "[WARNING] For the life of me, why do I have to make this change" )
		globalData.width = width , globalData.height = height;
#else
		if     ( globalData.periodicType==         NO_PERIODIC ) globalData.width = width+1, globalData.height = height+1;
		else if( globalData.periodicType==CYLINDRICAL_PERIODIC ) globalData.width = width  , globalData.height = height+1;
		else                                                     globalData.width = width  , globalData.height = height;
#endif
		int w = globalData.width , h = globalData.height , wBS , hBS;
		SetPaddedSize( w , h , minMGRes , globalData.width , _globalData.height , wBS , hBS );
		if( globalData.periodicType==CYLINDRICAL_PERIODIC && globalData.width!=w ) fprintf( stderr , "Cylindrical width was not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , globalData.width , globalData.height , wBS , hBS ) , exit( 0 );
		if( globalData.periodicType==  SPHERICAL_PERIODIC && ( globalData.width!=w || globalData.height!=h ) ) fprintf( stderr , "Spherical dimensions were not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , globalData.width , globalData.height , wBS , hBS ) , exit( 0 );
		int oocLevels = 1;

	
#if 1 // [WARNING] Shouldn't this computed relative to the globalData.width/height?
		while( (long long)( globalData.width >> oocLevels ) * (long long)( globalData.height >> oocLevels ) >= inCoreRes*inCoreRes ) oocLevels++;
		if( !ProcessPartition< int >::IsValidBandSize( globalData.width>>oocLevels , globalData.height>>oocLevels , globalData.iters , minBandSize ) )
		{
			fprintf( stderr , "[ERROR] In-core resolution is too small: %d x %d\n" , globalData.width>>oocLevels , globalData.height>>oocLevels );
			return false;
		}
#else
		while( (long long)( width >> oocLevels ) * (long long)( height >> oocLevels ) >= inCoreRes*inCoreRes ) oocLevels++;
		if( !ProcessPartition<int>::IsValidBandSize( width>>oocLevels , height>>oocLevels , globalData.iters , minBandSize ) )
		{
			fprintf( stderr , "[ERROR] In-core resolution is too small: %d x %d\n" , width>>oocLevels , height>>oocLevels );
			return false;
		}
#endif
		int blockSize = 4<<oocLevels;
#if 1	// Since the coarsest level is a multiple of four, and in the next finer one pixels have a width of two, we only need to
		// double that to ensure multipler of four.
		blockSize >>= 1;
#endif

		for( int i=0 ; i<clientCount-1 ; i++ ) if( clientSockets[i].clientData.width % blockSize )
		{
			fprintf( stderr , "Non-paddable client[%d] width not multiple of block size: %d mod %d != 0 \n" , i , clientSockets[i].clientData.width , blockSize );
			exit(0);
		}

		if( globalData.verbose )
		{
			printf( "  Image Dimension: %d x %d -> %d x %d\n" , width , height , globalData.width , globalData.height ) , fflush( stdout );
			printf( "In-Core Dimension: %d x %d\n" , _globalData.width>>oocLevels , globalData.height>>oocLevels ) , fflush( stdout );
			printf( "       Block Size: %d\n" , blockSize ) , fflush( stdout );
		}

		// Set the local information for the bands
		int start = 0 , index = 0;
		for( int i=0 ; i<clientCount ; i++ )
		{
			int w = clientSockets[i].clientData.width;
			clientSockets[i].clientData.index = index;
			index += clientSockets[i].clientData.subClients;
			clientSockets[i].clientData.start = start;
			clientSockets[i].clientData.end = clientSockets[i].clientData.cEnd = start+w;
			if( i==clientCount-1 ) clientSockets[i].clientData.end = globalData.width;
			start = clientSockets[i].clientData.end;
	
			SendOnSocket( clientSockets[i].socket , ( ConstPointer( GlobalData ) )GetPointer( globalData ) , sizeof( globalData ), "Failed to send global data to client" );
			SendOnSocket( clientSockets[i].socket , ( ConstPointer( ClientData ) )GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) ,"Failed to send client data to client" );
			SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( blockSize ), sizeof( blockSize ) , "Failed to send block size to client" );
			int getGradientAverage = removeAverageGradient ? 1 : 0;
			SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( getGradientAverage ) , sizeof( getGradientAverage ) , "Failed to send gradient average status to client" );
			int useGrayImage = PixelChannels==3 ? 0 : 1;
			SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( useGrayImage ) , sizeof( useGrayImage ) , "Failed to send use gray status to client" );
		}
		delete[] clientSockets;
	}

	int clipDimensions[2];
	clipDimensions[0] = globalData.cWidth;
	clipDimensions[1] = globalData.cHeight;

	{
		ClientSocket* clientSockets = new ClientSocket[clientCount];

		// Establish a connection to each of the processes
		for( int i=0 ; i<subClientCount ; i++ )
		{
			clientSockets[i].socket = AcceptSocket( listenSocket );
			if( globalData.verbose ) printf( "Connected to thread: %d /%d     \r" , i+1 , subClientCount );
		}
		if( globalData.verbose ) printf( "\n" ) , fflush( stdout );

		// Get the client info
		int width = 0 , height;
		for( int i=0 ; i<subClientCount ; i++ )
		{
			ReceiveOnSocket( clientSockets[i].socket , GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) , "Failed to get client data from client" );
			if( !i ) height = clientSockets[i].clientData.height;
			else if( height!=clientSockets[i].clientData.height ) fprintf( stderr , "[ERROR] Band heights differ: %d != %d\n" , height , clientSockets[i].clientData.height ) , exit(0);
			width += clientSockets[i].clientData.width;
		}
		{
			std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > > averageMap;
			// Receive average gradient values
			std::map< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > _averageMap;

			for( int i=0 ; i<subClientCount ; i++ )
			{
				int count;
				ReceiveOnSocket( clientSockets[i].socket , GetPointer( count ) , sizeof( count ) , "Failed to get average gradient count from client" );
				if( count )
				{
					std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > > __averageMap( count );
					long long size = sizeof( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) * count;
					ReceiveOnSocket( clientSockets[i].socket , GetPointer( __averageMap ) , (int)(size) ,  "Failed to get average gradients from client" );
					for( int j=0 ; j<__averageMap.size() ; j++ )
					{
						typename std::map< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > >::iterator iter = _averageMap.find( __averageMap[j].first );
						if( iter==_averageMap.end() ) _averageMap[ __averageMap[j].first ] = __averageMap[j].second;
						else iter->second += __averageMap[j].second;
					}
				}
			}

			for( typename std::map< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > >::iterator iter=_averageMap.begin() ; iter!=_averageMap.end() ; iter++ )
				averageMap.push_back( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > >( iter->first , iter->second ) );
	
			for( int i=0 ; i<averageMap.size() ; i++ ) averageMap[i].second.average();

			// Send back average gradient values
			for( int i=0 ; i<clientCount ; i++ )
			{
				int count = (int)( averageMap.size() );
				SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( count ) , int (sizeof( count ) ) , "Failed to send average gradient count to client" );
				if( count )
				{
					long long size = sizeof( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) * count;
					SendOnSocket( clientSockets[i].socket , ( ConstPointer( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) )GetPointer( averageMap ) , (int)(size) ,  "Failed to send average gradients to client" );
				}
			}
		}

	// Sort the clients by index order
	qsort( clientSockets , clientCount , sizeof(ClientSocket) , ClientSocket::Sort );

	// Store the clipping width and then compute the padded dimensions
	// Only compute the padded dimensions if we don't already know them.
	if( clipDimensions )
	{
		_globalData.cWidth	= clipDimensions[0];
		_globalData.cHeight	= clipDimensions[1];
		_globalData.width	=  width;
		_globalData.height	= height;
	}
	else
	{
		_globalData.cWidth  = width;
		_globalData.cHeight = height;

#if 0
#pragma message( "[WARNING] For the life of me, why do I have to make this change" )
		_globalData.width = width  , _globalData.height = height;
#else
		if     ( periodicType==         NO_PERIODIC ) _globalData.width = width+1, _globalData.height = height+1;
		else if( periodicType==CYLINDRICAL_PERIODIC ) _globalData.width = width  , _globalData.height = height+1;
		else                                          _globalData.width = width  , _globalData.height = height  ;
#endif

		int w = _globalData.width , h = _globalData.height , wBS, hBS;
		SetPaddedSize( w , h , _minMGRes, _globalData.width , _globalData.height , wBS , hBS );
		if( periodicType==CYLINDRICAL_PERIODIC && _globalData.width!=w ) fprintf( stderr , "Cylindrical width was not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , _globalData.width , _globalData.height , wBS , hBS ) , exit( 0 );
		if( periodicType==  SPHERICAL_PERIODIC && ( _globalData.width!=w || _globalData.height!=h ) ) fprintf( stderr , "Spherical dimensions were not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , _globalData.width , _globalData.height , wBS , hBS ) , exit( 0 );
	}
	// Send the general solver information
	for( int i=0 ; i<clientCount ; i++ ) SendOnSocket( clientSockets[i].socket , ( ConstPointer( GlobalData ) )GetPointer( _globalData ) , sizeof( _globalData ) , "Failed to send global data to client" );

	// Set the local information for the bands
	int start=0;
	for( int i=0 ; i<clientCount ; i++ )
	{
		int w = clientSockets[i].clientData.width;
		clientSockets[i].clientData.start = start;
		clientSockets[i].clientData.end = clientSockets[i].clientData.cEnd = start+w;
		if( i==clientCount-1 ) clientSockets[i].clientData.end  = _globalData.width;
		start = clientSockets[i].clientData.end;
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( ClientData ) )GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) , "Failed to send client data to client" );
	}

	// Figure out how many different out-of-core levels are needed
	_oocLevels = 1;
	while( (long long)( width>>_oocLevels ) * (long long)( height>>_oocLevels ) >= _inCoreRes*_inCoreRes ) _oocLevels++;
	if( !ProcessPartition< int >::IsValidBandSize( width>>_oocLevels , height>>_oocLevels , iters , minBandSize ) )
	{
		fprintf( stderr , "In-core resolution is too small: %d x %d\n" , width>>_oocLevels , height>>_oocLevels );
		return false;
	}
	// WARNING!!! It's possible that _oocLevels = 0 if the image fits in-core

	ProcessPartition< int > rootPartition;
	rootPartition.resize( clientCount );
	for( int i=0 ; i<clientCount ; i++ )
	{
		rootPartition[i].width = clientSockets[i].clientData.end - clientSockets[i].clientData.start;
		if( !ProcessPartition< int >::IsValidBandSize( rootPartition[i].width , _globalData.height , _globalData.iters , minBandSize ) )
		{
			fprintf( stderr , "Client[%d] has invalid dimension: %d x %d\n" , i , rootPartition[i].width , _globalData.height );
			return false;
		}
	}
	// [WARNING] The "repeat" flag should be turned off if processors don't double-solve the same level
	if( !processPartition.Initialize( rootPartition , _globalData.height , _globalData.iters , _oocLevels-1 , false , minBandSize , _globalData.verbose ) )
	{
		fprintf( stderr , "Failed to create process partition\n" );
		return false;
	}
	for( int i=0 ; i<clientCount ; i++ )
	{
		int idx = i;
		int tCount = 1;
		for( int j=1 ; j<processPartition.size() && !(idx&1) && ( (idx>>1)<processPartition[j].size() ); j++ ) tCount++ , idx>>=1;
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( tCount ) , sizeof( tCount ) , "failed to send tCount to client" );
		Pointer( ProcessData ) pd = AllocPointer< ProcessData >( tCount );
		pd[0].index = 0 ;
		pd[0].maxIndex = processPartition.size()-1;
		pd[0].offset = i;
		pd[0].maxOffset = clientCount-1;
		pd[0].startDepth = processPartition[0].startDepth;
		pd[0].endDepth = processPartition[0].endDepth;
		processPartition[0].setBounds( i , pd[0].start , pd[0].stop );
		pd[0].width = _globalData.width;
		pd[0].height = _globalData.height;
		pd[0].children = (int)( processPartition[0][i].children.size() );

		idx = i;
		for( int j=1 ; j<tCount ; j++ )
		{
			idx >>= 1;
			pd[j].index  = j;
			pd[j].maxIndex = processPartition.size()-1;
			pd[j].offset = idx;
			pd[j].maxOffset = processPartition[j].size()-1;
			pd[j].startDepth = processPartition[j].startDepth;
			pd[j].endDepth = processPartition[j].endDepth;
			processPartition[j].setBounds( idx , pd[j].start , pd[j].stop );
			pd[j].width = _globalData.width >> processPartition[j].startDepth;
			pd[j].height = _globalData.height >> processPartition[j].startDepth;
			pd[j].children = (int)( processPartition[j][idx].children.size() );
		}
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( ProcessData ) )pd , sizeof(ProcessData) * tCount , "Failed to send process data to client" );
		FreePointer( pd );
	}
	delete[] clientSockets;

	// At this point we transition our communications so that we are talking to the threads, not processes
	for( int i=0 ; i<processPartition.size() ; i++ )
		for( int j=0 ; j<processPartition[i].size() ; j++ )
		{
			// Establish a connection to the process
			int index , offset;
			Socket sock = AcceptSocket( listenSocket );
			ReceiveOnSocket( sock , GetPointer( index  ) , sizeof( index  ) , "Failed to get client index" );
			ReceiveOnSocket( sock , GetPointer( offset ) , sizeof( offset ) , "Failed to get client offset" );
			processPartition[index][offset].data.sock = sock;
			if( index )
			{
				ReceiveOnSocket ( sock , GetPointer( processPartition[index][offset].data.childAddressX ) , sizeof( IPData ) , "Failed to get child address X" );
				ReceiveOnSocket ( sock , GetPointer( processPartition[index][offset].data.childAddressB ) , sizeof( IPData ) , "Failed to get child address B" );
			}
			if( periodicType!=NO_PERIODIC || offset )
			{
				ReceiveOnSocket ( sock , GetPointer( processPartition[index][offset].data.leftAddress ) , sizeof( IPData ) , "Failed to get left address from client" );
			}
		}

	// Create a listening socket to establish the connection between the in-core and out-of-core solver
	int listenerCount = 2;
	ThreadHandle listenerHandles[3];

	_inCoreConnectInfo.init( processPartition[processPartition.size()-1].size() );

	IPData inCoreAddressX , inCoreAddressB , inCoreAddressLeft;

	{
		IOServer::SystemLock lock;
		GetHostEndpointAddress( &inCoreAddressX.address );
		GetHostEndpointAddress( &inCoreAddressB.address );
		GetHostEndpointAddress( &inCoreAddressLeft.address );
	}
	_inCoreConnectInfo.listenChildX = GetListenSocket( inCoreAddressX.port );
	_inCoreConnectInfo.listenChildB = GetListenSocket( inCoreAddressB.port );
	if ( _inCoreConnectInfo.listenChildX == _INVALID_ACCEPTOR_SOCKET_ ) return false;
	if ( _inCoreConnectInfo.listenChildB == _INVALID_ACCEPTOR_SOCKET_ ) return false;
	listenerHandles[0] = SpawnListenerSocket( _inCoreConnectInfo.listenChildX , _inCoreConnectInfo.childX , processPartition[processPartition.size()-1].size() );
	listenerHandles[1] = SpawnListenerSocket( _inCoreConnectInfo.listenChildB , _inCoreConnectInfo.childB , processPartition[processPartition.size()-1].size() );
	if( !TestThreadHandle( listenerHandles[0] ) ) return false;
	if( !TestThreadHandle( listenerHandles[1] ) ) return false;

	Socket leftSocket = _INVALID_SOCKET_;
	if( periodicType!=NO_PERIODIC )
	{
		_inCoreConnectInfo.listenLeft = GetListenSocket( inCoreAddressLeft.port );
		if ( _inCoreConnectInfo.listenLeft == _INVALID_ACCEPTOR_SOCKET_ ) return false;
		listenerCount++;
		listenerHandles[2] = SpawnListenerSocket( _inCoreConnectInfo.listenLeft , &leftSocket , 1 );
		if( !TestThreadHandle( listenerHandles[2] ) ) return false;
	}

	// Now send everybody information about who to connect to for their data
	for( int i=0 ; i<processPartition.size() ; i++ )
		for( int j=0 ; j<processPartition[i].size() ; j++ )
		{
			if( i )
				for( size_t k=0 ; k<processPartition[i][j].children.size() ; k++ )
				{
					Socket& sock = processPartition[i-1][processPartition[i][j].children[k]].data.sock;
					SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( processPartition[i][j].data.childAddressX ) , sizeof( IPData ) , "Failed to send child address X to client" );
					SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( processPartition[i][j].data.childAddressB ) , sizeof( IPData ) , "Failed to send child address B to client" );
				}
		}
	{
		int i = processPartition.size()-1;
		for( int j=0; j<processPartition[i].size() ; j++ )
		{
			Socket& sock = processPartition[i][j].data.sock;
			SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( inCoreAddressX ) , sizeof( IPData ) , "Failed to send address X to client" );
			SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( inCoreAddressB ) , sizeof( IPData ) , "Failed to send address B to client" );
		}
	}
	for( int i=0 ; i<processPartition.size() ; i++ )
		for( int j=0 ; j<processPartition[i].size() ; j++ )
			if( periodicType!=NO_PERIODIC || j )
			{
				Socket& sock = processPartition[i][(j-1+processPartition[i].size())%processPartition[i].size()].data.sock;
				SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( processPartition[i][j].data.leftAddress ) , sizeof( IPData ), "Failed to send left address to client" );
			}
	if( periodicType!=NO_PERIODIC )
		_inCoreConnectInfo.right = new SocketStream( inCoreAddressLeft.address , inCoreAddressLeft.port , 5 , false );
	WaitOnThreads( listenerHandles , listenerCount , 10000 , "SocketedMultigridServer::SetUp" );
	if( periodicType!=NO_PERIODIC ) _inCoreConnectInfo.left = new SocketStream( leftSocket );

	// Figure out how many different in-core resolutions are needed
	_icLevels = 1;
	int w = _globalData.width >> _oocLevels , h = _globalData.height >> _oocLevels;
	while( ProcessPartition<int>::IsValidBandSize( w>>_icLevels , h>>_icLevels , _globalData.iters , 0 ) && (w>>_icLevels)>=_minMGRes && (h>>_icLevels)>=_minMGRes )
		_icLevels++;

	// Now for the spherical synchronizers
	if( _globalData.periodicType==SPHERICAL_PERIODIC )
	{
		_icSynchronizers = new SphericalSynchronizer[ _icLevels ];
		if( !_icSynchronizers )
		{
			fprintf(stderr,"Failed to allocate in-core synchronizers\n");
			return false;
		}
		_oocSynchronizers.resize( processPartition.size() );
		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			if( !i ) _oocSynchronizers[i] = new SphericalSynchronizer[ (processPartition[i].endDepth-processPartition[i].startDepth+1)+1 ];
			else	 _oocSynchronizers[i] = new SphericalSynchronizer[ (processPartition[i].endDepth-processPartition[i].startDepth+1)   ];
		}

		int *oocWidths , icWidths;
		oocWidths = new int[clientCount];

		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			if(!i)
			{
				for( int j=0 ; j<processPartition[i].size() ; j++ ) oocWidths[j] = processPartition[i][j].width;
				imageSynchronizer.init( _globalData.width , oocWidths , clientCount , Degree , 0 , 1 );
			}
			for( int d=processPartition[i].startDepth ; d<=processPartition[i].endDepth ; d++ )
			{
				int depth = d-processPartition[i].startDepth;
				for( int j=0 ; j<processPartition[i].size() ; j++ ) oocWidths[j] = processPartition[i][j].width >> depth;
				if( _globalData.verbose )
					_oocSynchronizers[i][depth].init( _globalData.width>>d , oocWidths , processPartition[i].size() , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+Degree , _globalData.vCycles );
				else
					_oocSynchronizers[i][depth].init( _globalData.width>>d , oocWidths , processPartition[i].size() , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+     0 , _globalData.vCycles );
			}
		}
		delete[] oocWidths;
		for( int d=0 ; d<_icLevels ; d++ )
		{
			int depth = d+_oocLevels;
			icWidths = _globalData.width >> depth;
			if( _globalData.verbose )
				_icSynchronizers[d].init( icWidths , &icWidths , 1 , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+Degree , _globalData.vCycles );
			else
				_icSynchronizers[d].init( icWidths , &icWidths , 1 , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+     0 , _globalData.vCycles );
		}
		_synchronizerHandles = new ThreadHandle[ _oocLevels + 1 + _icLevels ];

		int idx = 0;
		{
			// [WARNING] The spherical synchronizer should only be enabled depending on whether or not there are
			// labels, not the interpolation weight
			if( _globalData.iWeight )
				_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< SyncType > , &imageSynchronizer );
			else
				_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< ImageData< PixelChannels , LabelChannels , SyncType , uint16_t > > , &imageSynchronizer );
			if( !TestThreadHandle( _synchronizerHandles[idx++] ) ) fprintf( stderr , "Failed to create sync thread\n" ) , exit(0);
		}
		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			int depths = processPartition[i].endDepth - processPartition[i].startDepth + 1;
			for( int d=0 ; d<depths ; d++ )
			{
				_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< SyncType > , &_oocSynchronizers[i][d] );
				if( !TestThreadHandle( _synchronizerHandles[idx++] ) ) fprintf( stderr , "Failed to create sync thread\n" ) , exit(0);
			}
		}
		for(int i=0;i<_icLevels;i++)
		{
			_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< float > , &_icSynchronizers[i] );
			if( !TestThreadHandle( _synchronizerHandles[idx++] ) ) fprintf( stderr , "Failed to create sync thread\n" ) , exit(0);
		}

		int* icPorts = new int[_icLevels];
		for( int d=0 ; d<_icLevels ; d++ )	icPorts[d] = _icSynchronizers[d].port(0);

		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			int depths = processPartition[i].endDepth - processPartition[i].startDepth + 1;
			Pointer( int ) oocPorts = AllocPointer< int >( depths );
			for( int j=0 ; j< processPartition[i].size() ; j++ )
			{
				if( !i )
				{
					int imagePort = imageSynchronizer.port(j);
					SendOnSocket( processPartition[i][j].data.sock , ( ConstPointer( int ) )GetPointer( imagePort ) , sizeof( int ) , "Failed to send image port to client" );
				}
				for( int d=0 ; d<depths ; d++ ) oocPorts[d] = _oocSynchronizers[i][d].port(j);
				SendOnSocket( processPartition[i][j].data.sock , ( ConstPointer( int ) )oocPorts , sizeof(int)*depths , "Failed to send out of core ports to client" );
			}
			FreePointer( oocPorts );
		}
		_icSyncSockets = new Socket[_icLevels];
		for( int i=0 ; i<_icLevels ; i++ )
		{
			_icSyncSockets[i] = GetConnectSocket( address , icPorts[i] );
			if( _icSyncSockets[i] == _INVALID_SOCKET_ )	return false;
		}
		delete[] icPorts;
	}
	CloseAcceptorSocket( listenSocket );
	return true;
	}
	//////////////////////
	SocketedMultigridServer< PixelChannels , LabelChannels , SyncType > server;
	if( !server.template SetUp< LabelType >
		(
		_address , listenSocket , subClientCount , _globalData.iters , _inCoreRes , _minMGRes , _globalData.vCycles , _minBandSize ,
		_globalData.tileWidth , _globalData.tileHeight , _globalData.tileExt ,
		_globalData.gammaCorrection , _globalData.quality , _globalData.lanes , _globalData.verbose , _globalData.periodicType , _globalData.iWeight , _globalData.lump , _globalData.gWeight , _globalData.gScale ,
		_removeAverageGradient , _globalData.unknownType , _globalData.showProgress ,
		_noCG , _globalData.shortSync , clipDimensions )
		) return false;
	server.Run();
	return true;

	/************************/
	if( globalData.periodicType==SPHERICAL_PERIODIC )
	{
		WaitOnThreads( synchronizerHandles , oocLevels+1+icLevels , 10000 , "SocketedMultigridServer::~SocketedMultigridServer" );
		for( int i=0 ; i<oocLevels+1+icLevels ; i++ ) CloseThreadHandle( synchronizerHandles[i] );
		delete[] synchronizerHandles;
		synchronizerHandles = NULL;
		if( icSyncSockets )
		{
			for( int i=0;i<icLevels ; i++ ) CloseSocket( icSyncSockets[i] );
			delete[] icSyncSockets;
		}
		icSyncSockets = NULL;
	}
	if( icSynchronizers )	delete[] icSynchronizers;
	for( sizet i=0 ; i<oocSynchronizers.size() ; i++ ) if( oocSynchronizers[i] ) delete[] oocSynchronizers[i];
	oocSynchronizers.clear();
	icSynchronizers = NULL;
}

#else // !MISHA_CODE_CLEAN_UP
/////////////////////////////
// SocketedMultigridServer //
/////////////////////////////
template< int PixelChannels , int LabelChannels , class SyncType >
SocketedMultigridServer< PixelChannels , LabelChannels , SyncType >::SocketedMultigridServer(void)
{
	_synchronizerHandles = NULL;
	_icSynchronizers = NULL;
	_icSyncSockets = NULL;
}
template< int PixelChannels , int LabelChannels , class SyncType >
SocketedMultigridServer< PixelChannels , LabelChannels , SyncType >::~SocketedMultigridServer(void)
{
	if( _globalData.periodicType==SPHERICAL_PERIODIC )
	{
		WaitOnThreads( _synchronizerHandles , _oocLevels+1+_icLevels , 10000 , "SocketedMultigridServer::~SocketedMultigridServer" );
		for( int i=0 ; i<_oocLevels+1+_icLevels ; i++ ) CloseThreadHandle( _synchronizerHandles[i] );
		delete[] _synchronizerHandles;
		_synchronizerHandles = NULL;
		if( _icSyncSockets )
		{
			for( int i=0;i<_icLevels ; i++ ) CloseSocket( _icSyncSockets[i] );
			delete[] _icSyncSockets;
		}
		_icSyncSockets = NULL;
	}
	if( _icSynchronizers )	delete[] _icSynchronizers;
	for( size_t i=0 ; i<_oocSynchronizers.size() ; i++ ) if( _oocSynchronizers[i] ) delete[] _oocSynchronizers[i];
	_oocSynchronizers.clear();
	_icSynchronizers = NULL;
}
template< int PixelChannels , int LabelChannels , class SyncType >
template< class LabelType >
bool SocketedMultigridServer< PixelChannels , LabelChannels , SyncType >::SetUp( char* address , int port , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles, int minBandSize ,
											   int tileWidth , int tileHeight , const char* tileExt , 
											   bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType , double iWeight , bool lump , double gWeight , double gScale ,
											   bool removeAverageGradient , int unknownType , bool showProgress , bool noCG , bool shortSync , int* clipDimensions )
{
	// Create a listening Socket for connecting to server
	Socket listenSocket = _INVALID_SOCKET_;
	listenSocket = GetListenSocket( port );
	if( listenSocket == _INVALID_SOCKET_ ) return false;
	if( address )
	{
		printf( "Server Address: %s:%d\n", address , port ) , fflush( stdout );
		return SetUp( address , listenSocket , clientCount , iters , inCoreRes , minMGRes , vCycles , minBandSize , tileWidth , tileHeight , tileExt , gammaCorrection , quality , lanes , verbose , periodicType , iWeight , lump , gWeight , gScale , removeAverageGradient , unknownType , showProgress , noCG , shortSync , clipDimensions );
	}
	else
	{
		char hostAddress[512];
		GetHostAddress( hostAddress );
		printf( "Host Address: %s:%d\n", hostAddress , port ) , fflush( stdout );
		return SetUp( hostAddress , listenSocket , clientCount , iters , inCoreRes , minMGRes , vCycles , minBandSize , tileWidth , tileHeight , tileExt , gammaCorrection , quality , lanes , verbose , periodicType , iWeight , lump ,  gWeight , gScale , removeAverageGradient , unknownType , showProgress , noCG , shortSync , clipDimensions );
	}
}
template< int PixelChannels , int LabelChannels , class SyncType >
template< class LabelType >
bool SocketedMultigridServer< PixelChannels , LabelChannels , SyncType >::SetUp( char* address , AcceptorSocket listenSocket , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles , int minBandSize ,
														   int tileWidth , int tileHeight , const char* tileExt , 
											   bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType , double iWeight , bool lump , double gWeight , double gScale ,
											   bool removeAverageGradient , int unknownType , bool showProgress , bool noCG , bool shortSync , int* clipDimensions )
{
	ClientSocket* clientSockets = new ClientSocket[clientCount];
	_noCG		= noCG;
	_inCoreRes	= inCoreRes;
	_minMGRes	= minMGRes;

	_globalData.periodicType= periodicType;
	_globalData.verbose		= verbose;
	_globalData.iters		= iters;
	_globalData.vCycles		= vCycles;
	_globalData.iWeight		= iWeight;
	_globalData.lump		= lump;
	_globalData.gWeight		= gWeight;
	_globalData.gScale		= gScale;
	_globalData.tileWidth	= tileWidth;
	_globalData.tileHeight	= tileHeight;
	_globalData.quality		= quality;
	_globalData.lanes		= lanes;
	_globalData.unknownType	= unknownType;
	_globalData.showProgress= showProgress;
	_globalData.shortSync	= shortSync;
	if( tileExt ) strcpy( _globalData.tileExt , tileExt );
	else          _globalData.tileExt[0] = 0;

	// Establish a connection to each of the processes
	for( int i=0 ; i<clientCount ; i++ )
	{
		clientSockets[i].socket = AcceptSocket( listenSocket );
		if( _globalData.verbose ) printf( "Connected to thread: %d /%d     \r" , i+1 , clientCount );
	}
	if( _globalData.verbose ) printf( "\n" ) , fflush( stdout );

	// Get the client info
	int width = 0 , height;
	for( int i=0 ; i<clientCount ; i++ )
	{
		ReceiveOnSocket( clientSockets[i].socket , GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) , "Failed to get client data from client" );
		if( !i ) height = clientSockets[i].clientData.height;
		else if( height!=clientSockets[i].clientData.height ) fprintf( stderr , "[ERROR] Band heights differ: %d != %d\n" , height , clientSockets[i].clientData.height ) , exit(0);
		width += clientSockets[i].clientData.width;
	}
	{
		std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > > averageMap;
		// Receive average gradient values
		std::map< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > _averageMap;

		for( int i=0 ; i<clientCount ; i++ )
		{
			int count;
			ReceiveOnSocket( clientSockets[i].socket , GetPointer( count ) , sizeof( count ) , "Failed to get average gradient count from client" );
			if( count )
			{
				std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > > __averageMap( count );
				long long size = sizeof( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) * count;
				ReceiveOnSocket( clientSockets[i].socket , GetPointer( __averageMap ) , (int)(size) ,  "Failed to get average gradients from client" );
				for( int j=0 ; j<__averageMap.size() ; j++ )
				{
					typename std::map< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > >::iterator iter = _averageMap.find( __averageMap[j].first );
					if( iter==_averageMap.end() ) _averageMap[ __averageMap[j].first ] = __averageMap[j].second;
					else iter->second += __averageMap[j].second;
				}
			}
		}
		for( typename std::map< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > >::iterator iter=_averageMap.begin() ; iter!=_averageMap.end() ; iter++ )
			averageMap.push_back( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > >( iter->first , iter->second ) );

		for( int i=0 ; i<averageMap.size() ; i++ ) averageMap[i].second.average();

		// Send back average gradient values
		for( int i=0 ; i<clientCount ; i++ )
		{
			int count = (int)( averageMap.size() );
			SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( count ) , int (sizeof( count ) ) , "Failed to send average gradient count to client" );
			if( count )
			{
				long long size = sizeof( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) * count;
				SendOnSocket( clientSockets[i].socket , ( ConstPointer( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) )GetPointer( averageMap ) , (int)(size) ,  "Failed to send average gradients to client" );
			}
		}
	}

	// Sort the clients by index order
	qsort( clientSockets , clientCount , sizeof(ClientSocket) , ClientSocket::Sort );

	// Store the clipping width and then compute the padded dimensions
	// Only compute the padded dimensions if we don't already know them.
	if( clipDimensions )
	{
		_globalData.cWidth	= clipDimensions[0];
		_globalData.cHeight	= clipDimensions[1];
		_globalData.width	=  width;
		_globalData.height	= height;
	}
	else
	{
		_globalData.cWidth  = width;
		_globalData.cHeight = height;

#if 0
#pragma message( "[WARNING] For the life of me, why do I have to make this change" )
		_globalData.width = width  , _globalData.height = height;
#else
		if     ( periodicType==         NO_PERIODIC ) _globalData.width = width+1, _globalData.height = height+1;
		else if( periodicType==CYLINDRICAL_PERIODIC ) _globalData.width = width  , _globalData.height = height+1;
		else                                          _globalData.width = width  , _globalData.height = height  ;
#endif

		int w = _globalData.width , h = _globalData.height , wBS, hBS;
		SetPaddedSize( w , h , _minMGRes, _globalData.width , _globalData.height , wBS , hBS );
		if( periodicType==CYLINDRICAL_PERIODIC && _globalData.width!=w ) fprintf( stderr , "Cylindrical width was not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , _globalData.width , _globalData.height , wBS , hBS ) , exit( 0 );
		if( periodicType==  SPHERICAL_PERIODIC && ( _globalData.width!=w || _globalData.height!=h ) ) fprintf( stderr , "Spherical dimensions were not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , _globalData.width , _globalData.height , wBS , hBS ) , exit( 0 );
	}
	// Send the general solver information
	for( int i=0 ; i<clientCount ; i++ ) SendOnSocket( clientSockets[i].socket , ( ConstPointer( GlobalData ) )GetPointer( _globalData ) , sizeof( _globalData ) , "Failed to send global data to client" );

	// Set the local information for the bands
	int start=0;
	for( int i=0 ; i<clientCount ; i++ )
	{
		int w = clientSockets[i].clientData.width;
		clientSockets[i].clientData.start = start;
		clientSockets[i].clientData.end = clientSockets[i].clientData.cEnd = start+w;
		if( i==clientCount-1 ) clientSockets[i].clientData.end  = _globalData.width;
		start = clientSockets[i].clientData.end;
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( ClientData ) )GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) , "Failed to send client data to client" );
	}

	// Figure out how many different out-of-core levels are needed
	_oocLevels = 1;
	while( (long long)( width>>_oocLevels ) * (long long)( height>>_oocLevels ) >= _inCoreRes*_inCoreRes ) _oocLevels++;
	if( !ProcessPartition< int >::IsValidBandSize( width>>_oocLevels , height>>_oocLevels , iters , minBandSize ) )
	{
		fprintf( stderr , "In-core resolution is too small: %d x %d\n" , width>>_oocLevels , height>>_oocLevels );
		return false;
	}
	// WARNING!!! It's possible that _oocLevels = 0 if the image fits in-core

	ProcessPartition< int > rootPartition;
	rootPartition.resize( clientCount );
	for( int i=0 ; i<clientCount ; i++ )
	{
		rootPartition[i].width = clientSockets[i].clientData.end - clientSockets[i].clientData.start;
		if( !ProcessPartition< int >::IsValidBandSize( rootPartition[i].width , _globalData.height , _globalData.iters , minBandSize ) )
		{
			fprintf( stderr , "Client[%d] has invalid dimension: %d x %d\n" , i , rootPartition[i].width , _globalData.height );
			return false;
		}
	}
	// [WARNING] The "repeat" flag should be turned off if processors don't double-solve the same level
	if( !processPartition.Initialize( rootPartition , _globalData.height , _globalData.iters , _oocLevels-1 , false , minBandSize , _globalData.verbose ) )
	{
		fprintf( stderr , "Failed to create process partition\n" );
		return false;
	}
	for( int i=0 ; i<clientCount ; i++ )
	{
		int idx = i;
		int tCount = 1;
		for( int j=1 ; j<processPartition.size() && !(idx&1) && ( (idx>>1)<processPartition[j].size() ); j++ ) tCount++ , idx>>=1;
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( tCount ) , sizeof( tCount ) , "failed to send tCount to client" );
		Pointer( ProcessData ) pd = AllocPointer< ProcessData >( tCount );
		pd[0].index = 0 ;
		pd[0].maxIndex = processPartition.size()-1;
		pd[0].offset = i;
		pd[0].maxOffset = clientCount-1;
		pd[0].startDepth = processPartition[0].startDepth;
		pd[0].endDepth = processPartition[0].endDepth;
		processPartition[0].setBounds( i , pd[0].start , pd[0].stop );
		pd[0].width = _globalData.width;
		pd[0].height = _globalData.height;
		pd[0].children = (int)( processPartition[0][i].children.size() );

		idx = i;
		for( int j=1 ; j<tCount ; j++ )
		{
			idx >>= 1;
			pd[j].index  = j;
			pd[j].maxIndex = processPartition.size()-1;
			pd[j].offset = idx;
			pd[j].maxOffset = processPartition[j].size()-1;
			pd[j].startDepth = processPartition[j].startDepth;
			pd[j].endDepth = processPartition[j].endDepth;
			processPartition[j].setBounds( idx , pd[j].start , pd[j].stop );
			pd[j].width = _globalData.width >> processPartition[j].startDepth;
			pd[j].height = _globalData.height >> processPartition[j].startDepth;
			pd[j].children = (int)( processPartition[j][idx].children.size() );
		}
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( ProcessData ) )pd , sizeof(ProcessData) * tCount , "Failed to send process data to client" );
		FreePointer( pd );
	}
	delete[] clientSockets;

	// At this point we transition our communications so that we are talking to the threads, not processes
	for( int i=0 ; i<processPartition.size() ; i++ )
		for( int j=0 ; j<processPartition[i].size() ; j++ )
		{
			// Establish a connection to the process
			int index , offset;
			Socket sock = AcceptSocket( listenSocket );
			ReceiveOnSocket( sock , GetPointer( index  ) , sizeof( index  ) , "Failed to get client index" );
			ReceiveOnSocket( sock , GetPointer( offset ) , sizeof( offset ) , "Failed to get client offset" );
			processPartition[index][offset].data.sock = sock;
			if( index )
			{
				ReceiveOnSocket ( sock , GetPointer( processPartition[index][offset].data.childAddressX ) , sizeof( IPData ) , "Failed to get child address X" );
				ReceiveOnSocket ( sock , GetPointer( processPartition[index][offset].data.childAddressB ) , sizeof( IPData ) , "Failed to get child address B" );
			}
			if( periodicType!=NO_PERIODIC || offset )
			{
				ReceiveOnSocket ( sock , GetPointer( processPartition[index][offset].data.leftAddress ) , sizeof( IPData ) , "Failed to get left address from client" );
			}
		}

	// Create a listening socket to establish the connection between the in-core and out-of-core solver
	int listenerCount = 2;
	ThreadHandle listenerHandles[3];

	_inCoreConnectInfo.init( processPartition[processPartition.size()-1].size() );

	IPData inCoreAddressX , inCoreAddressB , inCoreAddressLeft;

	{
		IOServer::SystemLock lock;
		GetHostEndpointAddress( &inCoreAddressX.address );
		GetHostEndpointAddress( &inCoreAddressB.address );
		GetHostEndpointAddress( &inCoreAddressLeft.address );
	}
	_inCoreConnectInfo.listenChildX = GetListenSocket( inCoreAddressX.port );
	_inCoreConnectInfo.listenChildB = GetListenSocket( inCoreAddressB.port );
	if ( _inCoreConnectInfo.listenChildX == _INVALID_ACCEPTOR_SOCKET_ ) return false;
	if ( _inCoreConnectInfo.listenChildB == _INVALID_ACCEPTOR_SOCKET_ ) return false;
	listenerHandles[0] = SpawnListenerSocket( _inCoreConnectInfo.listenChildX , _inCoreConnectInfo.childX , processPartition[processPartition.size()-1].size() );
	listenerHandles[1] = SpawnListenerSocket( _inCoreConnectInfo.listenChildB , _inCoreConnectInfo.childB , processPartition[processPartition.size()-1].size() );
	if( !TestThreadHandle( listenerHandles[0] ) ) return false;
	if( !TestThreadHandle( listenerHandles[1] ) ) return false;

	Socket leftSocket = _INVALID_SOCKET_;
	if( periodicType!=NO_PERIODIC )
	{
		_inCoreConnectInfo.listenLeft = GetListenSocket( inCoreAddressLeft.port );
		if ( _inCoreConnectInfo.listenLeft == _INVALID_ACCEPTOR_SOCKET_ ) return false;
		listenerCount++;
		listenerHandles[2] = SpawnListenerSocket( _inCoreConnectInfo.listenLeft , &leftSocket , 1 );
		if( !TestThreadHandle( listenerHandles[2] ) ) return false;
	}

	// Now send everybody information about who to connect to for their data
	for( int i=0 ; i<processPartition.size() ; i++ )
		for( int j=0 ; j<processPartition[i].size() ; j++ )
		{
			if( i )
				for( size_t k=0 ; k<processPartition[i][j].children.size() ; k++ )
				{
					Socket& sock = processPartition[i-1][processPartition[i][j].children[k]].data.sock;
					SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( processPartition[i][j].data.childAddressX ) , sizeof( IPData ) , "Failed to send child address X to client" );
					SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( processPartition[i][j].data.childAddressB ) , sizeof( IPData ) , "Failed to send child address B to client" );
				}
		}
	{
		int i = processPartition.size()-1;
		for( int j=0; j<processPartition[i].size() ; j++ )
		{
			Socket& sock = processPartition[i][j].data.sock;
			SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( inCoreAddressX ) , sizeof( IPData ) , "Failed to send address X to client" );
			SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( inCoreAddressB ) , sizeof( IPData ) , "Failed to send address B to client" );
		}
	}
	for( int i=0 ; i<processPartition.size() ; i++ )
		for( int j=0 ; j<processPartition[i].size() ; j++ )
			if( periodicType!=NO_PERIODIC || j )
			{
				Socket& sock = processPartition[i][(j-1+processPartition[i].size())%processPartition[i].size()].data.sock;
				SendOnSocket ( sock , ( ConstPointer( IPData ) )GetPointer( processPartition[i][j].data.leftAddress ) , sizeof( IPData ), "Failed to send left address to client" );
			}
	if( periodicType!=NO_PERIODIC )
		_inCoreConnectInfo.right = new SocketStream( inCoreAddressLeft.address , inCoreAddressLeft.port , 5 , false );
	WaitOnThreads( listenerHandles , listenerCount , 10000 , "SocketedMultigridServer::SetUp" );
	if( periodicType!=NO_PERIODIC ) _inCoreConnectInfo.left = new SocketStream( leftSocket );

	// Figure out how many different in-core resolutions are needed
	_icLevels = 1;
	int w = _globalData.width >> _oocLevels , h = _globalData.height >> _oocLevels;
	while( ProcessPartition<int>::IsValidBandSize( w>>_icLevels , h>>_icLevels , _globalData.iters , 0 ) && (w>>_icLevels)>=_minMGRes && (h>>_icLevels)>=_minMGRes )
		_icLevels++;

	// Now for the spherical synchronizers
	if( _globalData.periodicType==SPHERICAL_PERIODIC )
	{
		_icSynchronizers = new SphericalSynchronizer[ _icLevels ];
		if( !_icSynchronizers )
		{
			fprintf(stderr,"Failed to allocate in-core synchronizers\n");
			return false;
		}
		_oocSynchronizers.resize( processPartition.size() );
		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			if( !i ) _oocSynchronizers[i] = new SphericalSynchronizer[ (processPartition[i].endDepth-processPartition[i].startDepth+1)+1 ];
			else	 _oocSynchronizers[i] = new SphericalSynchronizer[ (processPartition[i].endDepth-processPartition[i].startDepth+1)   ];
		}

		int *oocWidths , icWidths;
		oocWidths = new int[clientCount];

		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			if(!i)
			{
				for( int j=0 ; j<processPartition[i].size() ; j++ ) oocWidths[j] = processPartition[i][j].width;
				imageSynchronizer.init( _globalData.width , oocWidths , clientCount , Degree , 0 , 1 );
			}
			for( int d=processPartition[i].startDepth ; d<=processPartition[i].endDepth ; d++ )
			{
				int depth = d-processPartition[i].startDepth;
				for( int j=0 ; j<processPartition[i].size() ; j++ ) oocWidths[j] = processPartition[i][j].width >> depth;
				if( _globalData.verbose )
					_oocSynchronizers[i][depth].init( _globalData.width>>d , oocWidths , processPartition[i].size() , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+Degree , _globalData.vCycles );
				else
					_oocSynchronizers[i][depth].init( _globalData.width>>d , oocWidths , processPartition[i].size() , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+     0 , _globalData.vCycles );
			}
		}
		delete[] oocWidths;
		for( int d=0 ; d<_icLevels ; d++ )
		{
			int depth = d+_oocLevels;
			icWidths = _globalData.width >> depth;
			if( _globalData.verbose )
				_icSynchronizers[d].init( icWidths , &icWidths , 1 , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+Degree , _globalData.vCycles );
			else
				_icSynchronizers[d].init( icWidths , &icWidths , 1 , (_globalData.iters*Degree+Degree)+Degree , (_globalData.iters*Degree+Degree)+     0 , _globalData.vCycles );
		}
		_synchronizerHandles = new ThreadHandle[ _oocLevels + 1 + _icLevels ];

		int idx = 0;
		{
			// [WARNING] The spherical synchronizer should only be enabled depending on whether or not there are
			// labels, not the interpolation weight
			if( _globalData.iWeight )
				_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< SyncType > , &imageSynchronizer );
			else
				_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< ImageData< PixelChannels , LabelChannels , SyncType , uint16_t > > , &imageSynchronizer );
			if( !TestThreadHandle( _synchronizerHandles[idx++] ) ) fprintf( stderr , "Failed to create sync thread\n" ) , exit(0);
		}
		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			int depths = processPartition[i].endDepth - processPartition[i].startDepth + 1;
			for( int d=0 ; d<depths ; d++ )
			{
				_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< SyncType > , &_oocSynchronizers[i][d] );
				if( !TestThreadHandle( _synchronizerHandles[idx++] ) ) fprintf( stderr , "Failed to create sync thread\n" ) , exit(0);
			}
		}
		for(int i=0;i<_icLevels;i++)
		{
			_synchronizerHandles[idx] = RunThread( SphericalSynchronizer::RunThread< float > , &_icSynchronizers[i] );
			if( !TestThreadHandle( _synchronizerHandles[idx++] ) ) fprintf( stderr , "Failed to create sync thread\n" ) , exit(0);
		}

		int* icPorts = new int[_icLevels];
		for( int d=0 ; d<_icLevels ; d++ )	icPorts[d] = _icSynchronizers[d].port(0);

		for( int i=0 ; i<processPartition.size() ; i++ )
		{
			int depths = processPartition[i].endDepth - processPartition[i].startDepth + 1;
			Pointer( int ) oocPorts = AllocPointer< int >( depths );
			for( int j=0 ; j< processPartition[i].size() ; j++ )
			{
				if( !i )
				{
					int imagePort = imageSynchronizer.port(j);
					SendOnSocket( processPartition[i][j].data.sock , ( ConstPointer( int ) )GetPointer( imagePort ) , sizeof( int ) , "Failed to send image port to client" );
				}
				for( int d=0 ; d<depths ; d++ ) oocPorts[d] = _oocSynchronizers[i][d].port(j);
				SendOnSocket( processPartition[i][j].data.sock , ( ConstPointer( int ) )oocPorts , sizeof(int)*depths , "Failed to send out of core ports to client" );
			}
			FreePointer( oocPorts );
		}
		_icSyncSockets = new Socket[_icLevels];
		for( int i=0 ; i<_icLevels ; i++ )
		{
			_icSyncSockets[i] = GetConnectSocket( address , icPorts[i] );
			if( _icSyncSockets[i] == _INVALID_SOCKET_ )	return false;
		}
		delete[] icPorts;
	}
	CloseAcceptorSocket( listenSocket );
	return true;
}
template< int PixelChannels , int LabelChannels , class SyncType >
void SocketedMultigridServer< PixelChannels , LabelChannels , SyncType >::Run( void )
{
#if MISHA_DENORMAL_CONTROL
	_MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
#endif // MISHA_DENORMAL_CONTROL
#if DEBUG_SOCKETS
	int testCount = 0;
	auto testCommunication = [&]( const char* message )
	{
		char str[512];
		printf( "Testing: (%s)\n" , message );
		printf( "\tReceiving... \n" );
		for( int i=0 ; i<processPartition[0].size() ; i++ ) ReceiveOnSocket( processPartition[0][i].data.sock , str , (size_t)512 , "Failed to get test string %d" , testCount );
		printf( "\tReceived: %s\n" , str );
		sprintf( str , "test %d" , testCount );
		printf( "\tSending: %s\n" , str );
		for( int i=0 ; i<processPartition[0].size() ; i++ ) SendOnSocket( processPartition[0][i].data.sock , (const char*)str , (size_t)512 , "Failed to send test string %d" , testCount );
		testCount++;
	};
#endif // DEBUG_SOCKETS


	StreamingGrid *inB , *outX;
	Pointer( SolverInfo< PixelChannels > )  solverInfo = AllocPointer< SolverInfo< PixelChannels > >( _oocLevels );
	Pointer( SolverInfo< PixelChannels > ) tSolverInfo = AllocPointer< SolverInfo< PixelChannels > >( _oocLevels );
	double zeroAverage[PixelChannels];
	for( int c=0 ; c<PixelChannels ; c++ )	zeroAverage[c]=0;
	std::vector< float > globalB , globalX;
	int gWidth  = _globalData.width  >> _oocLevels;
	int gHeight = _globalData.height >> _oocLevels;

	globalB.resize( gWidth * gHeight * PixelChannels );
	globalX.resize( gWidth * gHeight * PixelChannels );


	memset( _average , 0 , sizeof( _average ) );
	inB  = GetMultiSocketBackedGrid< SyncType , PixelChannels >( _inCoreConnectInfo.childB , _inCoreConnectInfo.children , gHeight    );
	outX = GetMultiSocketBackedGrid< SyncType , PixelChannels >( _inCoreConnectInfo.childX , _inCoreConnectInfo.children , gHeight<<1 );

	for( int ii=0 ; ii<_globalData.vCycles ; ii++ )
	{
#if DEBUG_SOCKETS
		testCommunication( "Starting V-Cycle" );
#endif // DEBUG_SOCKETS
		double t;

		t=Time();
		memset( solverInfo , 0 , sizeof( SolverInfo< PixelChannels > ) * _oocLevels );
		inB->reset( true , 1 );

		for( int j=0 ; j<gHeight ; j++ )
		{
			Pointer( SyncType ) row = ( Pointer( SyncType ) )(*inB)[j];
			for( int i=0 ; i<PixelChannels*gWidth ; i++ ) globalB[ gWidth*j*PixelChannels + i ] = float( row[i] );
			inB->advance();
		}
#if DEBUG_SOCKETS
		testCommunication( "Receiving solver info" );
#endif // DEBUG_SOCKETS
		for( int i=0 ; i<processPartition.size() ; i++ )
			for( int j=0 ; j<processPartition[i].size() ; j++ )
			{
				ReceiveOnSocket( processPartition[i][j].data.sock , tSolverInfo + processPartition[i].startDepth , sizeof(SolverInfo<PixelChannels>) * (processPartition[i].endDepth-processPartition[i].startDepth+1) , "Failed to get process data from process" );
				for( int d=processPartition[i].startDepth ; d<=processPartition[i].endDepth ; d++ )
				{
					solverInfo[d].bSquareNorm += tSolverInfo[d].bSquareNorm;
					solverInfo[d].rSquareNorm += tSolverInfo[d].rSquareNorm;
					solverInfo[d].xSquareNorm += tSolverInfo[d].xSquareNorm;
					for( int c=0 ; c<PixelChannels ; c++ ) solverInfo[d].solutionSum[c] += tSolverInfo[d].solutionSum[c];
				}
			}
#if DEBUG_SOCKETS
		testCommunication( "Receiving average" );
#endif // DEBUG_SOCKETS
		if( ii==0 )
		{
			for( int j=0 ; j<processPartition[0].size() ; j++ )
			{
				double avg[PixelChannels];
				ReceiveOnSocket( processPartition[0][j].data.sock , GetPointer( avg ) , sizeof( avg ) , "Failed to get average from process" );
				for( int c=0 ; c<PixelChannels ; c++ )	_average[c] += avg[c];
			}
			for( int c=0 ; c<PixelChannels ; c++ )	_average[c] /= _globalData.cWidth , _average[c] /= _globalData.cHeight;
		}
		if( _globalData.verbose )
		{
			printf( "Client Restriction:    %f\n" , Time()-t ) , fflush( stdout );
			for( int i=0 ; i<processPartition.size() ; i++ )
				for( int d=processPartition[i].startDepth ; d<=processPartition[i].endDepth ; d++ )
				{
					int dd = processPartition[i].endDepth - ( d - processPartition[i].startDepth );
					printf( "\tError[%d x %d] %g -> %g\t(%.8g)\n" , _globalData.width >> d , _globalData.height >> d , sqrt(solverInfo[dd].bSquareNorm) , sqrt(solverInfo[dd].rSquareNorm) , sqrt(solverInfo[dd].xSquareNorm) ) , fflush( stdout );
				}
		}
#if DEBUG_SOCKETS
		testCommunication( "Finished restriction" );
#endif // DEBUG_SOCKETS

		// Get the base solution
		DotProductStencil dotMajor , d2DotMajor , dotMinor ,d2DotMinor;
		SetDownSampledStencil( _globalData.width  , _oocLevels , dotMajor , d2DotMajor , _globalData.lump );
		SetDownSampledStencil( _globalData.height , _oocLevels , dotMinor , d2DotMinor , _globalData.lump );
		if( ii==_globalData.vCycles-1 ) SolveInCore( dotMajor , d2DotMajor , dotMinor , d2DotMinor , globalB , globalX , _average);
		else							SolveInCore( dotMajor , d2DotMajor , dotMinor , d2DotMinor , globalB , globalX , zeroAverage);

		t=Time();
		int gWidth2 = gWidth<<1 , gHeight2 = gHeight<<1;
		outX->reset( false , 1 );
		// Get the up-sampled solution
		float *upX = new float[ PixelChannels * gWidth2 * gHeight2 ];

		for( int c=0 ; c<PixelChannels ; c++ )
		{
			std::vector< float > low , high;
			low.resize( gWidth * gHeight );
			for( int j=0 ; j<gHeight ; j++ ) for( int i=0 ; i<gWidth ; i++ ) low[ j*gWidth + i ] = globalX[ j*gWidth*PixelChannels + i*PixelChannels + c ];
			if     ( _globalData.periodicType==  SPHERICAL_PERIODIC ) UpSampleSphericalSolution  ( low , gWidth , gHeight , high );
			else if( _globalData.periodicType==CYLINDRICAL_PERIODIC ) UpSampleCylindricalSolution( low , gWidth , gHeight , high );
			else                                                      UpSampleRectangularSolution( low , gWidth , gHeight , high );
			for( int j=0 ; j<gHeight2 ; j++ ) for( int i=0 ; i<gWidth2 ; i++ ) upX[ j*gWidth2*PixelChannels + i*PixelChannels + c ] = high[ j*gWidth2 + i ];
		}

#if DEBUG_SOCKETS
		testCommunication( "Starting prolongation" );
#endif // DEBUG_SOCKETS

		// Send the up-sampled solution to the socket
		for( int j=0 ; j<gHeight2 ; j++ )
		{
			Pointer( SyncType ) row = ( Pointer( SyncType ) )(*outX)[j];
			for( int i=0 ; i<PixelChannels*gWidth2 ; i++ ) row[i] = SyncType( upX[ gWidth2*j*PixelChannels+i ] );
			outX->advance();
		}
		delete upX;

#if DEBUG_SOCKETS
		testCommunication( "Receiving solver info" );
#endif // DEBUG_SOCKETS
		memset( solverInfo , 0 , sizeof(SolverInfo< PixelChannels >) * _oocLevels );
		for( int i=0 ; i<processPartition.size() ; i++ )
			for( int j=0 ; j<processPartition[i].size() ; j++ )
			{
				ReceiveOnSocket ( processPartition[i][j].data.sock , tSolverInfo + processPartition[i].startDepth , sizeof(SolverInfo<PixelChannels>) * (processPartition[i].endDepth-processPartition[i].startDepth+1) , "Failed to get process data from process" );
				for( int d=processPartition[i].startDepth ; d<=processPartition[i].endDepth ; d++ )
				{
					solverInfo[d].bSquareNorm += tSolverInfo[d].bSquareNorm;
					solverInfo[d].rSquareNorm += tSolverInfo[d].rSquareNorm;
					solverInfo[d].xSquareNorm += tSolverInfo[d].xSquareNorm;
					for( int c=0 ; c<PixelChannels ; c++ ) solverInfo[d].solutionSum[c] += tSolverInfo[d].solutionSum[c];
				}
			}
		if( _globalData.verbose )
		{
			printf( "Client Prolongation:    %f\n" , Time()-t ) , fflush( stdout );
			for( int i=processPartition.size()-1 ; i>=0 ; i-- )
				for( int d=processPartition[i].endDepth ; d>=processPartition[i].startDepth ; d-- )
				{
					int dd = processPartition[i].startDepth - ( d - processPartition[i].endDepth );
					printf("\tError[%d x %d] %g -> %g\t(%.8g)\n" , _globalData.width >> d , _globalData.height >> d , sqrt(solverInfo[dd].bSquareNorm) , sqrt(solverInfo[dd].rSquareNorm) , sqrt(solverInfo[dd].xSquareNorm) ) , fflush( stdout );
				}
		}
#if DEBUG_SOCKETS
		testCommunication( "Done prolongation" );
#endif // DEBUG_SOCKETS
	}
	FreePointer( solverInfo  );
	FreePointer( tSolverInfo );
	delete inB;
	delete outX;
}
template< int PixelChannels , int LabelChannels , class SyncType >
void SocketedMultigridServer< PixelChannels , LabelChannels , SyncType >::SolveInCore( DotProductStencil& dotMajor , DotProductStencil& d2DotMajor ,
																 DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
																 std::vector< float >& in , std::vector< float >& out , double average[PixelChannels] )
{
	double t;
	int w = _globalData.width >> _oocLevels , h = _globalData.height >> _oocLevels;
	SocketedMultiGridStreamingSolver< PixelChannels , float , float >* solvers;
	DotProductStencil lowDotMajor , lowD2DotMajor , lowDotMinor , lowD2DotMinor;

	SetDownSampledStencil( dotMajor , d2DotMajor , _globalData.width  >> _oocLevels , _icLevels-1 , lowDotMajor , lowD2DotMajor );
	SetDownSampledStencil( dotMinor , d2DotMinor , _globalData.height >> _oocLevels , _icLevels-1 , lowDotMinor , lowD2DotMinor );
	std::vector< float > lowB , lowX;
	StreamingGrid *B , *X;
	solvers = new SocketedMultiGridStreamingSolver< PixelChannels , float , float >[_icLevels];

	for( int i=1 ; i<_icLevels   ; i++ ) solvers[i].parent = &solvers[i-1];
	for( int i=0 ; i<_icLevels-1 ; i++ ) solvers[i].rChild = solvers[i].pChild = &solvers[i+1];
	for( int i=0 ; i<_icLevels   ; i++ )
	{
		solvers[i].bSquareNorm = solvers[i].rSquareNorm = solvers[i].xSquareNorm = 0;
		for( int c=0 ; c<PixelChannels ; c++ ) solvers[i].solutionSum[c] = 0;
		solvers[i].setResidual = true;
	}
	// BADNESS!!! For spherical domains, the server socket has not been set here.
#if LOW_MEMORY_HACK // WARNING!!!: Big hack over here: To avoid memory overflow we are using an ooc solver here
	MultiStreamIOServer server;
	solvers[_icLevels-1].Initialize( dotMajor , d2DotMajor , dotMinor , d2DotMinor , _globalData.iWeight , _globalData.gWeight , 0 , w , w , h , _globalData.iters , _inCoreConnectInfo.left , _icSyncSockets , _inCoreConnectInfo.right , true , _globalData.spherical , &server );
#else // !LOW_MEMORY_HACK
	solvers[_icLevels-1].Initialize( dotMajor , d2DotMajor , dotMinor , d2DotMinor , _globalData.iWeight , _globalData.gWeight , 0 , w , w , h , _globalData.iters , _inCoreConnectInfo.left , _icSyncSockets , _inCoreConnectInfo.right , false , _globalData.periodicType , NULL );
#endif // LOW_MEMORY_HACK

	lowX.resize( solvers[0].major*solvers[0].minor*PixelChannels );
	lowB.resize( solvers[0].major*solvers[0].minor*PixelChannels );
	B = new MemoryBackedGrid( ( Pointer( byte ) )GetPointer( lowB ) , solvers[0].major * PixelChannels * sizeof(float) , solvers[0].minor );
	X = new MemoryBackedGrid( ( Pointer( byte ) )GetPointer( lowX ) , solvers[0].major * PixelChannels * sizeof(float) , solvers[0].minor );
	solvers[_icLevels-1].inB = new MemoryBackedGrid( ( Pointer( byte ) )GetPointer(  in ) , w*PixelChannels*sizeof(float) , h );
	solvers[_icLevels-1].inX = new MemoryBackedGrid( ( Pointer( byte ) )GetPointer( out ) , w*PixelChannels*sizeof(float) , h );

	solvers[0].outX = X;
	solvers[0].outB = B;

	// Data for the interleaved streaming multigrid
	t=Time();
	// Initialize
	solvers[_icLevels-1].InitRestriction();
	solvers[_icLevels-1].SetRestriction();
#if LOW_MEMORY_HACK
	server.Reset();
	server.StartIO();
#endif // LOW_MEMORY_HACK
	solvers[_icLevels-1].SolveRestriction();
#if LOW_MEMORY_HACK
	server.Reset();
#endif // LOW_MEMORY_HACK

	// Clean up
	delete solvers[_icLevels-1].inB ; solvers[_icLevels-1].inB = NULL;
	delete solvers[_icLevels-1].inX ; solvers[_icLevels-1].inX = NULL;

	if( _globalData.verbose )
	{
		printf( "Server Restriction:        %f\n" , Time()-t ) , fflush( stdout );
		for( int i=_icLevels-1 ; i>=0 ; i-- )
			printf( "\tError[%d x %d] %g -> %g\n" , solvers[i].major , solvers[i].minor , sqrt(solvers[i].bSquareNorm) , sqrt(solvers[i].rSquareNorm)) , fflush( stdout );
	}

	solvers[_icLevels-1].UnSetRestriction();
	// Get the base solution
	SparseMatrix< double > lap;
	t=Time();
	if     ( _globalData.periodicType==SPHERICAL_PERIODIC )
		SetSphericalLaplacianMatrix  ( lowDotMajor , lowD2DotMajor , lowDotMinor , lowD2DotMinor , lap,solvers[0].major , solvers[0].minor , _globalData.iWeight , _globalData.gWeight );
	else if( _globalData.periodicType==CYLINDRICAL_PERIODIC )
		SetCylindricalLaplacianMatrix( lowDotMajor , lowD2DotMajor , lowDotMinor , lowD2DotMinor , lap,solvers[0].major , solvers[0].minor , _globalData.iWeight , _globalData.gWeight );
	else
		SetRectangularLaplacianMatrix( lowDotMajor , lowD2DotMajor , lowDotMinor , lowD2DotMinor , lap,solvers[0].major , solvers[0].minor , _globalData.iWeight , _globalData.gWeight );
	double squareRNorm = 0 , squareBNorm = 0;
	{
		std::vector< double > myLowX , myLowB , temp;
		myLowB.resize( solvers[0].major * solvers[0].minor );
		myLowX.resize( solvers[0].major * solvers[0].minor );
		temp.resize  ( solvers[0].major * solvers[0].minor );
		for( int c=0 ; c<PixelChannels ; c++ )
		{
			for( int i=0 ; i<solvers[0].major ; i++ ) for( int j=0 ; j<solvers[0].minor ; j++ )
			{
				myLowX[i+j*solvers[0].major] = lowX[ c + ( i + j*solvers[0].major ) * PixelChannels ];
				myLowB[i+j*solvers[0].major] = lowB[ c + ( i + j*solvers[0].major ) * PixelChannels ];
			}
			{
				Multiply< false >( lap , ( ConstPointer( double ) )GetPointer( myLowX ) , GetPointer( temp ) );
#if USE_PARALLEL_CG
#pragma omp parallel for
#endif // USE_PARALLEL_CG
				for( int i=0 ; i<myLowB.size() ; i++ ) myLowB[i] -= temp[i] , myLowX[i] = 0.;
			}
			if( !_noCG )
			{
				if( _globalData.iWeight==0 ) MySolveConjugateGradient< true  >( lap , ( ConstPointer( double ) )GetPointer( myLowB ) , GetPointer( myLowX ) , 20*(int)( sqrt(lap.groups+1.0) ) );
				else						 MySolveConjugateGradient< false >( lap , ( ConstPointer( double ) )GetPointer( myLowB ) , GetPointer( myLowX ) , 20*(int)( sqrt(lap.groups+1.0) ) );
				if( _globalData.verbose )
				{
					Multiply< false >( lap , ( ConstPointer( double ) )GetPointer( myLowX ) , GetPointer( temp ) );
					for( int i=0 ; i<myLowB.size() ; i++ ) squareBNorm += myLowB[i] * myLowB[i] , squareRNorm += (myLowB[i]-temp[i]) * (myLowB[i]-temp[i]);
				}
			}

			for( int i=0 ; i<solvers[0].major ; i++) for( int j=0 ; j<solvers[0].minor ; j++ )
			{
				lowX[ c + ( i + j*solvers[0].major) * PixelChannels ] = float( myLowX[i+j*solvers[0].major] );
				myLowX[i+j*solvers[0].major] = 0.0;
			}
		}
	}
	if( _globalData.verbose )
	{
		printf( "Server Conjugate Gradient: %f\n" , Time()-t ) , fflush( stdout );
		printf( "\tError[%d x %d] %g -> %g\n" , solvers[0].major , solvers[0].minor , sqrt(squareBNorm) , sqrt(squareRNorm) ) , fflush( stdout );
	}

	for( int i=0 ; i<_icLevels ; i++ )
	{
		solvers[i].bSquareNorm = solvers[i].rSquareNorm = solvers[i].xSquareNorm=0;
		for( int c=0 ; c<PixelChannels ; c++ ) solvers[i].solutionSum[c] = 0;
		solvers[i].setResidual = _globalData.verbose;
	}
	// Solve the prolongation
	t=Time();
	solvers[0].inX = X;
	// Set the child dependencies
	solvers[_icLevels-1].outX = new MemoryBackedGrid( ( Pointer( byte ) )GetPointer( out ) , w*PixelChannels*sizeof(float) , h );
	solvers[_icLevels-1].InitProlongation();
	solvers[0].SetProlongation();
	// Solve
#if LOW_MEMORY_HACK
	server.StartIO();
#endif // LOW_MEMORY_HACK
	solvers[0].SolveProlongation();
#if LOW_MEMORY_HACK
	server.Reset();
#endif // LOW_MEMORY_HACK
	if( _globalData.verbose )
	{
		printf( "Server Prolongation:       %f\n" , Time()-t ) , fflush( stdout );
		for( int i=0 ; i<_icLevels ; i++ ) printf( "\tError[%d x %d] %g -> %g\n" , solvers[i].major , solvers[i].minor , sqrt(solvers[i].bSquareNorm) , sqrt(solvers[i].rSquareNorm) ) , fflush( stdout );
	}
	solvers[0].UnSetProlongation();
	delete solvers[_icLevels-1].outX;
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;

	int ww=( (long long)(_globalData.cWidth) *w+_globalData.width -1)/_globalData.width;
	int hh=( (long long)(_globalData.cHeight)*h+_globalData.height-1)/_globalData.height;
	for( int c=0 ; c<PixelChannels ; c++ )
	{
		double newAverage = 0;
		int pCount = 0;

		for( int j=0 ; j<hh ; j++ ) for( int i=0 ; i<ww ; i++ ) newAverage += out[(j*w+i)*PixelChannels+c] , pCount++;

		newAverage /= pCount;
		newAverage = average[c]-newAverage;

		if( _globalData.iWeight==0 ) for( int j=0 ; j<h ; j++ ) for( int i=0 ; i<w ; i++ ) out[(j*w+i)*PixelChannels+c] += float( newAverage );
	}
}

//////////////////////////////////
// SocketedSuperMultigridServer //
//////////////////////////////////
bool SocketedSuperMultigridServer::SetUp( char* prefix , int port , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles , int minBandSize ,
												    int tileWidth , int tileHeight , const char* tileExt , 
												    bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType ,
													double iWeight , bool lump , double gWeight , double gScale ,
													bool removeAverageGradient , int unknownType , bool showProgress ,
												    bool noCG , bool shortSync )
{
	_noCG = noCG;
	_port = port;
	_clientCount = clientCount;
	_inCoreRes = inCoreRes;
	_minMGRes = minMGRes;
	_minBandSize = minBandSize;
	if( _minBandSize > _inCoreRes ) _inCoreRes = _minBandSize;
	_removeAverageGradient = removeAverageGradient;
	_globalData.iters = iters;
	_globalData.vCycles = vCycles;
	_globalData.tileWidth = tileWidth;
	_globalData.tileHeight = tileHeight;
	if( tileExt ) strcpy( _globalData.tileExt , tileExt );
	else _globalData.tileExt[0] = 0;
	_globalData.gammaCorrection = gammaCorrection;
	_globalData.quality = quality;
	_globalData.lanes = lanes;
	_globalData.verbose = verbose;
	_globalData.periodicType = periodicType;
	_globalData.iWeight = iWeight;
	_globalData.lump = lump;
	_globalData.gWeight = gWeight;
	_globalData.gScale = gScale;
	_globalData.unknownType = unknownType;
	_globalData.showProgress = showProgress;
	_globalData.shortSync = shortSync;

	GetHostAddress( _address , prefix );

	return true;
}
template< int PixelChannels , int LabelChannels , class SyncType , class LabelType >
bool SocketedSuperMultigridServer::Run( void )
{
	int subClientCount = 0;

	// Create a listening SOCKET for connecting to server
	AcceptorSocket listenSocket = _INVALID_ACCEPTOR_SOCKET_;
	listenSocket = GetListenSocket( _port );
	if ( listenSocket == _INVALID_ACCEPTOR_SOCKET_ ) return false;
	printf( "Server Address: %s:%d\n", _address , _port ) , fflush( stdout );

	ClientSocket* clientSockets = new ClientSocket[_clientCount];
	int coreCount = 0;

	// Establish a connection to the clients
	for( int i=0 ; i<_clientCount ; i++ )
	{
		clientSockets[i].socket  = AcceptSocket( listenSocket );
		if( _globalData.verbose ) printf( "Connected to process: %d /%d     \r" , i+1 , _clientCount );
	}
	if( _globalData.verbose ) printf( "\n" ) , fflush( stdout );
	// Get the client info
	int width = 0 , height;
	for( int i=0 ; i<_clientCount ; i++ )
	{
		ReceiveOnSocket( clientSockets[i].socket , GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) , "Failed to get client data from client" );
		if(!i) height = clientSockets[i].clientData.height;
		else if( height!=clientSockets[i].clientData.height ) fprintf( stderr , "Band heights differ: %d != %d\n" , height , clientSockets[i].clientData.height ) , exit(0);
		width += clientSockets[i].clientData.width;
		subClientCount += clientSockets[i].clientData.subClients;
	}
	// Sort the clients by index order
	qsort( clientSockets , _clientCount , sizeof(ClientSocket) , ClientSocket::Sort );

	// Store the clipping width and then compute the padded dimensions
	_globalData.cWidth  = width;
	_globalData.cHeight = height;

#if 0
#pragma message( "[WARNING] For the life of me, why do I have to make this change" )
	_globalData.width = width , _globalData.height = height;
#else
	if     ( _globalData.periodicType==         NO_PERIODIC ) _globalData.width = width+1, _globalData.height = height+1;
	else if( _globalData.periodicType==CYLINDRICAL_PERIODIC ) _globalData.width = width  , _globalData.height = height+1;
	else                                                      _globalData.width = width  , _globalData.height = height;
#endif
	int w = _globalData.width , h = _globalData.height , wBS , hBS;
	SetPaddedSize( w , h , _minMGRes , _globalData.width , _globalData.height , wBS , hBS );
	if( _globalData.periodicType==CYLINDRICAL_PERIODIC && _globalData.width!=w ) fprintf( stderr , "Cylindrical width was not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , _globalData.width , _globalData.height , wBS , hBS ) , exit( 0 );
	if( _globalData.periodicType==  SPHERICAL_PERIODIC && ( _globalData.width!=w || _globalData.height!=h ) ) fprintf( stderr , "Spherical dimensions were not preserved: (%d , %d) -> (%d , %d)\t(%d , %d)\n" , w , h , _globalData.width , _globalData.height , wBS , hBS ) , exit( 0 );
	int oocLevels = 1;

	
#if 1 // [WARNING] Shouldn't this computed relative to the globalData.width/height?
	while( (long long)( _globalData.width >> oocLevels ) * (long long)( _globalData.height >> oocLevels ) >= _inCoreRes*_inCoreRes ) oocLevels++;
	if( !ProcessPartition< int >::IsValidBandSize( _globalData.width>>oocLevels , _globalData.height>>oocLevels , _globalData.iters , _minBandSize ) )
	{
		fprintf( stderr , "In-core resolution is too small: %d x %d\n" , _globalData.width>>oocLevels , _globalData.height>>oocLevels );
		return false;
	}
#else
	while( (long long)( width >> oocLevels ) * (long long)( height >> oocLevels ) >= _inCoreRes*_inCoreRes ) oocLevels++;
	if( !ProcessPartition<int>::IsValidBandSize( width>>oocLevels , height>>oocLevels , _globalData.iters , _minBandSize ) )
	{
		fprintf( stderr , "In-core resolution is too small: %d x %d\n" , width>>oocLevels , height>>oocLevels );
		return false;
	}
#endif
	int blockSize = 4<<oocLevels;
#if 1	// Since the coarsest level is a multiple of four, and in the next finer one pixels have a width of two, we only need to
		// double that to ensure multipler of four.
	blockSize >>= 1;
#endif

	for( int i=0 ; i<_clientCount-1 ; i++ )
		if( clientSockets[i].clientData.width % blockSize )
		{
			fprintf( stderr , "Non-paddable client[%d] width not multiple of block size: %d mod %d != 0 \n" , i , clientSockets[i].clientData.width , blockSize );
			exit(0);
		}
	if( _globalData.verbose )
	{
		printf( "  Image Dimension: %d x %d -> %d x %d\n" , width , height , _globalData.width , _globalData.height ) , fflush( stdout );
		printf( "In-Core Dimension: %d x %d\n" , _globalData.width>>oocLevels , _globalData.height>>oocLevels ) , fflush( stdout );
		printf( "       Block Size: %d\n" , blockSize ) , fflush( stdout );
	}

	// Set the local information for the bands
	int start = 0 , index = 0;
	for( int i=0 ; i<_clientCount ; i++ )
	{
		int w = clientSockets[i].clientData.width;
		clientSockets[i].clientData.index = index;
		index += clientSockets[i].clientData.subClients;
		clientSockets[i].clientData.start = start;
		clientSockets[i].clientData.end = clientSockets[i].clientData.cEnd = start+w;
		if( i==_clientCount-1 ) clientSockets[i].clientData.end  = _globalData.width;
		start = clientSockets[i].clientData.end;
	
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( GlobalData ) )GetPointer( _globalData ) , sizeof( _globalData ), "Failed to send global data to client" );
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( ClientData ) )GetPointer( clientSockets[i].clientData ) , sizeof( ClientData ) ,"Failed to send client data to client" );
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( blockSize ), sizeof( blockSize ) , "Failed to send block size to client" );
		int getGradientAverage = _removeAverageGradient ? 1 : 0;
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( getGradientAverage ) , sizeof( getGradientAverage ) , "Failed to send gradient average status to client" );
		int useGrayImage = PixelChannels==3 ? 0 : 1;
		SendOnSocket( clientSockets[i].socket , ( ConstPointer( int ) )GetPointer( useGrayImage ) , sizeof( useGrayImage ) , "Failed to send use gray status to client" );
	}
	delete[] clientSockets;
	int clipDimensions[2];
	clipDimensions[0] = _globalData.cWidth;
	clipDimensions[1] = _globalData.cHeight;
	SocketedMultigridServer< PixelChannels , LabelChannels , SyncType > server;
	if( !server.template SetUp< LabelType >
		(
		_address , listenSocket , subClientCount , _globalData.iters , _inCoreRes , _minMGRes , _globalData.vCycles , _minBandSize ,
		_globalData.tileWidth , _globalData.tileHeight , _globalData.tileExt ,
		_globalData.gammaCorrection , _globalData.quality , _globalData.lanes , _globalData.verbose , _globalData.periodicType , _globalData.iWeight , _globalData.lump , _globalData.gWeight , _globalData.gScale ,
		_removeAverageGradient , _globalData.unknownType , _globalData.showProgress ,
		_noCG , _globalData.shortSync , clipDimensions )
		) return false;
	server.Run();
	return true;
}
#endif // MISHA_CODE_CLEAN_UP