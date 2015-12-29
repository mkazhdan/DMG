#ifndef MULTIGRID_PROCESS_INCLUDED
#define MULTIGRID_PROCESS_INCLUDED

#include "LaplacianMatrix/SocketedStreamingSolver.h"
#include "LaplacianMatrix/SocketedMultigrid/SocketData.h"

#define STREAMING_GRID_BUFFER_MULTIPLIER 2


template< class PData >
class ProcessPartition
{
public:
	class MyData
	{
	public:
		int width;
		std::vector< int > children;
		PData data;
	};

	int startDepth , endDepth;
	std::vector< MyData > processData;

	void setBounds( int idx , int& start , int& stop );
	int size( void ) const;
	void resize( int );
	MyData& operator[] ( const int& idx );
	const MyData& operator[] ( const int& idx ) const;
	static bool IsValidBandSize( int width , int height , int iters , int minSize );
};

template< class PData >
class ProcessPartitionHierarchy
{
public:
	std::vector< ProcessPartition< PData > > levels;
	int size( void ) const;
	void resize( int );
	ProcessPartition< PData >& operator[] ( const int& idx );
	const ProcessPartition< PData >& operator[] ( const int& idx ) const;
	template< class TData >
	bool Initialize( const ProcessPartition< TData >& initialPartition , int height , int iters , int depths , bool repeat , int minSize , bool showPartition );
	int leaves( const int& depth , const int& offset ) const;
};

class ProcessingBlockData
{
public:
	ProcessData pData;
	Socket serverSocket;
	Pointer( Socket ) syncSockets;
	DataStream *leftStream , *rightStream;
	StreamingGrid *outHighX , *outHighP , *inHighB , *outLowR , *inLowX;
	ProcessingBlockData( void )
	{
		outHighX = outHighP = inHighB = outLowR = inLowX = NULL;
		syncSockets = NullPointer< Socket >();
		serverSocket = _INVALID_SOCKET_;
		leftStream = rightStream = NULL;
	}
	~ProcessingBlockData( void )
	{
#if 0
		if( outHighX ) delete outHighX , outHighX = NULL;
		if( outHighP ) delete outHighP , outHighP = NULL;
		if( inHighB )  delete inHighB  , inHighB  = NULL;
		if( outLowR )  delete outLowR  , outLowR  = NULL;
		if( inLowX )   delete inLowX   , inLowX   = NULL;
		if( leftStream )  delete leftStream  , leftStream  = NULL;
		if( rightStream ) delete rightStream , rightStream = NULL;
#endif
	}
};

template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
class MultigridThread
{
	Pointer( SolverInfo< PixelChannels > ) _solverInfo;

	SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >* _sDivergence;
	SocketedMultiGridRestrictionNode< PixelChannels >* _sRestriction;
	Pointer( Pointer( SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType > ) ) _solvers;
	int _vCycles;
	bool _verbose , _outOfCore;
	int _periodicType;
	Pointer( ProcessingBlockData ) _blockData;
	int _blockCount;	// Specifies how the hierarchy of depths is decomposed (depending on which other threads merge in at coarser levels)
	void _init(	Pointer( ProcessingBlockData ) blockData , int blockCount , const class GlobalData& globalData , bool showProgress , bool inCore , bool lowPixels , bool pixels , bool labels , bool outGuess );
	MultiStreamIOServer* multiStreamIOServer;
public:
	StreamingGrid *lowPixels , *pixels , *labels , *outGuess;
	MultigridThread( void );
	~MultigridThread( void );

	void Initialize(
		const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& gradientAverage ,
		MultiStreamIOServer* ioServer ,
		Pointer( ProcessingBlockData ) blockData , int blockCount ,
		const GlobalData& globalData , bool showProgress , bool inCore );

	static int RunThread( void* );
};

#include "MultigridThread.inl"
#endif // MULTIGRID_PROCESS_INCLUDED