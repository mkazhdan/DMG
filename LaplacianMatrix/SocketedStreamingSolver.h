#ifndef SOCKETED_STREAMING_SOLVER_INCLUDED
#define SOCKETED_STREAMING_SOLVER_INCLUDED

#include <map>

#define STREAMING_GRID_BUFFER_MULTIPLIER 2
#define ALIGNMENT 16

#define SAME_SIZE_BUFFERS 1
#define MISHA_DENORMAL_CONTROL 1
#define FIX_RESTRICTION


enum
{
	UNKNOWN_NONE ,
	UNKNOWN_BLACK ,
	UNKNOWN_HARMONIC
};

#define MyModIndex( idx , mod ) ( ModIndex( idx , mod ) )


#include <emmintrin.h>
#include "LinearAlgebra/SparseMatrix.h"
#include "LaplacianMatrix/LaplacianMatrix1D.h"
#include "LaplacianMatrix/LaplacianMatrix2D.h"
#include "Util/Socket.h"
#include "Util/GridStream.h"
#include "Util/MultiStreamIO.h"
#include "Util/Array.h"

const int Degree=2;
const int Type=ZERO_DERIVATIVE;
const int RealPerWord=sizeof(__m128)/sizeof(float);
const int WordPerDegree = (Degree+Degree+RealPerWord-1)/RealPerWord;


template< int Channels >
class AverageColor
{
	double values[Channels];
public:
	double& operator[] ( int idx ){ return values[idx]; }
};

class SocketedStreamData
{
protected:
	DataStream *leftStream , *rightStream;
	Socket syncXSocket , syncRSocket;
	int _start128,_end128,_size128,_paddedSize128,_padSize128;	// Sizes in terms of __m128's
	int _start,_end,_size,_paddedSize,_padSize;					// Sizes in terms of float's
	bool _set( int width , int height , int start , int end , DataStream* leftStream , Socket syncSocket , DataStream* rightStream , int periodicType );
public:
	bool SetSocketedStreamData( int width , int height , int start , int end ,int iters , DataStream* leftStream , Socket syncSocket , DataStream* rightStream , int periodicType );
	bool SetSocketedStreamData( int width , int height , int start , int end ,DataStream* leftStream , Socket syncSocket , DataStream* rightStream , int periodicType );
	int start	(void)	const;
	int end		(void)	const;
	int size	(void)	const;
};

class TemplateSSE
{
public:
	ALIGN( __m128 matrixValues[4][5] , 16);
	float diagonalR[4];
};
typedef FiniteElements2D< double , Type , Degree >::FullMatrixStencil MatrixStencil;
typedef FiniteElements1D< double , Type , Degree >::DotProduct< Type , Degree >::FullDotProductStencil DotProductStencil;

void SetDownSampledStencil( int dim , int iters , DotProductStencil& outDot , DotProductStencil& outD2Dot , bool lump );
void SetDownSampledStencil( const DotProductStencil& inDot , const DotProductStencil& inD2Dot, int dim , int iters , DotProductStencil& outDot , DotProductStencil& outD2Dot );

template< int Channels , class SyncType >
class SocketedStreamingSolver : public SocketedStreamData
{
protected:
	using SocketedStreamData::_start;
	using SocketedStreamData::_size;
private:
	Pointer( SyncType ) syncBuffer;
	Pointer( Pointer( __m128 ) ) RStream;
	Pointer( Pointer( __m128 ) ) BStream;
	Pointer( Pointer( __m128 ) ) XStream;
	Pointer( __m128 ) localXAccum[Degree];
	Pointer( __m128 ) RBuffer[Degree*Channels];
	Pointer( __m128 ) XBuffer[Degree*Channels];

	void freeStreams(void);

	void SyncSolverLeft ( int idx , bool read , bool overlapped=true );
	void SyncSolverRight( int idx , bool read , bool overlapped=true );
	void SyncSolverHead( int idx , bool read );
	void SyncSolverTail( int idx , bool read );
	void SyncResidualHead( int idx , bool read , bool syncPeriodic );
	void SyncResidualTail( int idx , bool read , bool syncPeriodic );
	void Solve        ( int idx , int c , int sSolve , int eSolve );
	void SolveInterior( int idx , int c , int sSolve , int eSolve );
	void SolveReverse        ( int idx , int c , int sSolve , int eSolve );
	void SolveInteriorReverse( int idx , int c , int sSolve , int eSolve );
	void SetResidual        ( int idx , int c , int sSolve , int eSolve );
	void SetInteriorResidual( int idx , int c , int sSolve , int eSolve );
protected:
	bool _deleteServer;
	MultiStreamIOServer *_server;
	int periodicType;
	bool clearX , clearB;
	int xSize,bSize,rSize,iters,index;
	Pointer( TemplateSSE ) lapTemplates;
	Pointer( TemplateSSE ) zeroLapTemplates;


	static int OffsetX( int iters );
	static int OffsetB( int iters );
	static int OffsetR( void );
public:
	int laneNum;
	int progressCount;
	int pastIndex;
	double startProgressTime , pastProgressTime;
	bool showProgress;
	int major,minor;
	float laplacianScale,laplacianScaleR;
	bool setResidual;
	double bSquareNorm , rSquareNorm , xSquareNorm , solutionSum[Channels];

	SocketedStreamingSolver(void);
	~SocketedStreamingSolver(void);

	void Init( int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , Socket syncSocket , DataStream *rightStream ,
		int periodicType , MultiStreamIOServer* server
		);
	void Init( const MatrixStencil& lStencil , int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , Socket syncSocket , DataStream *rightStream ,
		int periodicType , MultiStreamIOServer* server
		);

	void Set(int rSize=1);
	void Set(int start,int bStart,int xStart,int bEnd,int xEnd,int rSize=1);
	void UnSet(void);
	Pointer( float ) GetRRow( int row , int channel );
	Pointer( float ) GetXRow( int row , int channel );
	Pointer( float ) GetBRow( int row , int channel );
	bool Increment(void);
	template< class StorageType > bool UpdateXInput ( StreamingGrid* X );
	template< class StorageType > bool UpdateBInput ( StreamingGrid* B );
	template< class StorageType > bool UpdateXOutput( StreamingGrid* X );
	template< class StorageType > bool UpdateBOutput( StreamingGrid* B );

	void Solve(void);
	void SetResidual(void);
};

template< int Channels >
class SocketedMultiGridRestrictionNode
{
public:
	virtual void InitRestriction			( void ){;}													// Generic Initializer
	virtual void SetRestriction				( void ){;}													// Sets up interleaved dependencies
	virtual void SetRestriction				( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict ){;}	// Set parent's row
	virtual void SolveRestriction			( void ){;}
	virtual bool IterateRestriction			( void ){ return false; }
};

template< int Channels , class StorageType , class SyncType >
class SocketedMultiGridStreamingSolver : public SocketedStreamingSolver< Channels , SyncType > , public SocketedMultiGridRestrictionNode< Channels >
{
	typedef typename FiniteElements1D< double , Type , Degree >::FullProlongationStencil ProlongationStencil;
	typedef typename FiniteElements1D< double , Type , Degree >::FullRestrictionStencil RestrictionStencil;
	class ProlongationStencilSSE
	{
	public:
		ALIGN( __m128 matrixValues[2][4] , 16 );
	};
	class ProlongationStencilSSE2
	{
	public:
		ALIGN( __m128 matrixValues[2] , 16 );
	};
	Pointer( __m128 ) localRAccum;
	Pointer( ProlongationStencilSSE2 ) prolongationStencil;
	ProlongationStencil majorProlongationStencil,minorProlongationStencil;
	RestrictionStencil majorRestrictionStencil,minorRestrictionStencil;
	int prolongationOffset;
	int startProlongation,startRestriction;

	int restrictionBit;
	StreamingGrid *_B , *_X;	// Defined to be storage type by default
	bool inCore;
protected:
	using SocketedStreamData::_start;
	using SocketedStreamData::_size;
	using SocketedStreamData::_start128;
	using SocketedStreamData::_size128;
	using SocketedStreamData::_paddedSize128;
	using SocketedStreamData::_padSize128;
	using SocketedStreamingSolver< Channels , SyncType >::_server;
	using SocketedStreamingSolver< Channels , SyncType >::periodicType;
	using SocketedStreamingSolver< Channels , SyncType >::index;
	using SocketedStreamingSolver< Channels , SyncType >::iters;
	using SocketedStreamingSolver< Channels , SyncType >::xSize;
	using SocketedStreamingSolver< Channels , SyncType >::bSize;
	using SocketedStreamingSolver< Channels , SyncType >::rSize;
	using SocketedStreamingSolver< Channels , SyncType >::OffsetX;
	using SocketedStreamingSolver< Channels , SyncType >::OffsetB;
	using SocketedStreamingSolver< Channels , SyncType >::OffsetR;
public:
	using SocketedStreamingSolver< Channels , SyncType >::minor;
	using SocketedStreamingSolver< Channels , SyncType >::major;
	using SocketedStreamingSolver< Channels , SyncType >::laplacianScale;
	using SocketedStreamingSolver< Channels , SyncType >::laplacianScaleR;
	using SocketedStreamingSolver< Channels , SyncType >::setResidual;
	using SocketedStreamingSolver< Channels , SyncType >::laneNum;
	using SocketedStreamingSolver< Channels , SyncType >::Increment;

	DotProductStencil dotMajor , d2DotMajor , dotMinor , d2DotMinor;

	// There is some messiness in terms of how in/out are defined and used (e.g. do they increment or do they set?)
	StreamingGrid *inX , *outX , *inB , *outB , *outR , *outP;
	Pointer( float ) scratchR;
	Pointer( float ) scratchP;

	SocketedMultiGridRestrictionNode< Channels > *rChild;
	SocketedMultiGridStreamingSolver *parent , *pChild;

	SocketedMultiGridStreamingSolver( void );
	~SocketedMultiGridStreamingSolver( void );

	void Initialize(
		double iWeight , bool lump , double gWeight , int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , Socket* syncSockets , DataStream* rightStream ,
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
		);
	void Initialize(
		DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor , double iWeight , double gWeight ,
		int start , int end , int major , int minor , int iters,
		DataStream* leftStream , Socket* syncSockets , DataStream* rightStream ,
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
		);
	void InitProlongation	( void );
	void InitRestriction	( void );

	void SetProlongation	( void );
	void SetRestriction		( void );
	void UnSetProlongation	( void );
	void UnSetRestriction	( void );
	bool IterateProlongation( void );
	bool IterateRestriction	( void );
	void SolveProlongation	( void );
	void SolveRestriction	( void );
	void SetRestriction			( Pointer( float ) lowB , int c , int idx , int sRestrict , int eRestrict );
#ifdef FIX_RESTRICTION
#else // !FIX_RESTRICTION
	void SetInteriorRestriction ( Pointer( float ) lowB , int c , int idx , int sRestrict , int eRestrict );
#endif // FIX_RESTRICTION
	void SetProlongation		( Pointer( float ) highX , int c , int highIdx , int highStart , int highSize , int highMajor , int highMinor );
};

template< int PixelChannels , int LabelChannels , class PixelType , class LabelType >
class ImageData
{
public:
	PixelType pixel[PixelChannels];
	PixelType lowPixel[PixelChannels];
	LabelType label[LabelChannels];
};

template< class LabelType , int Channels >
class LabelData
{
public:
	LabelType l[Channels];
	inline bool operator == (const LabelData& ld) const;
	static LabelData UnknownLabel( void );
	operator size_t () const;
	bool operator() ( const LabelData& ld1 , const LabelData& ld2 ) const;
	bool operator < ( const LabelData& ld ) const;
};
template< int Channels >
struct GradientAverage
{
	long long dxCount , dyCount;
	double dx[Channels] , dy[Channels];
	GradientAverage( void ) { dxCount = dyCount = 0 , memset( dx , 0 , sizeof(double)*Channels ) , memset( dy , 0 , sizeof(double)*Channels ); }
	GradientAverage& operator *= ( double s )
	{
		for( int c=0 ; c<Channels ; c++ ) dx[c] *= s , dy[c] *= s;
		return *this;
	}
	GradientAverage& operator += ( const GradientAverage& a )
	{
		for( int c=0 ; c<Channels ; c++ ) dx[c] += a.dx[c] , dy[c] += a.dy[c];
		dxCount += a.dxCount , dyCount += a.dyCount;
		return *this;
	}
	void average( void )
	{
		for( int c=0 ; c<Channels ; c++ )
		{
			if( dxCount ) dx[c] /= dxCount;
			if( dyCount ) dy[c] /= dyCount;
		}
		dxCount = dyCount = 1;
	}
};
template< int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void SetGradientAverageMap( StreamingGrid* pixelStream , StreamingGrid* labelStream , std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& averageMap , bool spherical , bool showProgress );


template< int PixelChannels , int LabelChannels , class PixelType , class LabelType , class StorageType , class SyncType >
class SocketedStreamingDivergence : public SocketedStreamData , public SocketedMultiGridRestrictionNode< PixelChannels >
{
	static const int ISize  = 3*Degree+2;
	static const int DXSize = 3*Degree+1;
	static const int DYSize = 3*Degree;

	double _iWeight , _gWeight , _gScale;
	bool _lump;
	int _periodicType;
	int index;
	Pointer( LabelData< LabelType , LabelChannels > ) labels[ISize];
	Pointer( __m128 ) pixels[ISize*PixelChannels];
	Pointer( __m128 ) lowPixels[ISize*PixelChannels];
	bool _separateLaplacianComputation;
	std::map< class LabelData< LabelType , LabelChannels > , struct GradientAverage< PixelChannels > >  _gradientAverageMap;

	Pointer( __m128 ) dx[DXSize*PixelChannels];
	Pointer( __m128 ) dy[DYSize*PixelChannels];
	Pointer( float ) GetPixelRow( int row , int channel );
	Pointer( float ) GetLowPixelRow( int row , int channel );
	Pointer( LabelData< LabelType , LabelChannels > ) GetLabelRow( int row );
	Pointer( float ) GetDXRow( int row , int channel );
	Pointer( float ) GetDYRow( int row , int channel );
	Pointer( ImageData< PixelChannels , LabelChannels , SyncType , LabelType > ) syncBuffer;

	void SyncImageLeft ( int idx , bool read );
	void SyncImageRight( int idx , bool read );
	void SyncImageHead ( int idx , bool read );
	void SyncImageTail ( int idx,  bool read );
	StreamingGrid *lowPixelStream , *pixelStream , *labelStream;

	Pointer( __m128 ) localDMajorAccum;
	Pointer( __m128 ) localDMinorAccum;

	typename FiniteElements1D< float , Type,Degree >::DotProduct< Type , Degree >::FullDotProductStencil dotMajorStencil , dotMinorStencil;
	typename FiniteElements1D< float , DERIVATIVE(Type) , Degree-1 >::DotProduct< Type , Degree >::FullDotProductStencil dDotMajorStencil , dDotMinorStencil;
	typename FiniteElements2D< float , Type,Degree >::FullDivergenceStencil divergenceStencil;
	Pointer( __m128 ) localPAccum[Degree];
	Pointer( TemplateSSE ) lapTemplates;

	void _setPartials ( int y );
	void _setPartialsX( int y , int start , int end );
	void _setPartialsY( int y , int start , int end );
	void _setPartialsX( int y );
	void _setPartialsY( int y );
public:
	int unknownType;
	int major , minor;
	AverageColor< PixelChannels > average;

	SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType > *parent;
	SocketedStreamingDivergence( void );
	~SocketedStreamingDivergence( void );

	void Initialize( StreamingGrid* lowPixels , StreamingGrid *pixels , StreamingGrid* labels , double iWeight , bool lump , double gWeight , double gScale , int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , Socket* syncSockets , DataStream* rightStream ,
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server ,
		const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >* gradientAverageMap
		);

	void InitRestriction		( void );
	void SetRestriction			( void );
	void UnSetRestriction		( void );
	void SetRestriction			( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void _AddDerivativeRestriction			( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void _AddInteriorDerivativeRestriction	( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void _AddStencilRestriction				( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void _AddInteriorStencilRestriction		( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void SolveRestriction		( void );
	bool IterateRestriction		( void );
};

#include "StreamingSolver128.inl"
#include "SocketedStreamingSolver.inl"
#endif // SOCKETED_STREAMING_SOLVER_INCLUDED
