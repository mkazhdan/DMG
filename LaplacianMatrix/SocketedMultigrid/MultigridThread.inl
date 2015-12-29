#include "Util/MultiStreamIO.h"
#define DEBUG_SOCKETS 0

//////////////////////
// ProcessPartition //
//////////////////////
template< class PData >
bool ProcessPartition< PData >::IsValidBandSize( int width , int height , int iters , int minSize )
{
	if( height< 2*Degree ) return false;
	if( width < 16*(iters+2)*WordPerDegree || width<minSize ) return false;
#if 1 // [QUESTION] Don't we require that the row width be a multiple of four?
	if( width&3 ) return false;
#endif
	return true;
}
template< class PData >
void ProcessPartition< PData >::setBounds( int idx , int& start , int& stop )
{
	start = 0;
	for( int i=0 ; i<idx ; i++ ) start += processData[i].width;
	stop = start + processData[idx].width;
}
template< class PData > int ProcessPartition< PData >::size( void ) const { return int( processData.size() ); }
template< class PData > void ProcessPartition< PData >::resize( int sz ) { processData.resize( sz ); }
template< class PData > typename ProcessPartition< PData >::MyData& ProcessPartition< PData >::operator[] ( const int& idx ) { return processData[idx]; }
template< class PData > const typename ProcessPartition< PData >::MyData& ProcessPartition< PData >::operator[] ( const int& idx ) const { return processData[idx]; }

///////////////////////////////
// ProcessPartitionHierarchy //
///////////////////////////////

template< class PData > int ProcessPartitionHierarchy< PData >::size( void ) const { return int( levels.size() ); }
template< class PData > void ProcessPartitionHierarchy< PData >::resize( int sz ){ levels.resize( sz ); }
template< class PData > ProcessPartition< PData >& ProcessPartitionHierarchy< PData >::operator[] ( const int& idx ) { return levels[idx]; }
template< class PData > const ProcessPartition< PData >& ProcessPartitionHierarchy< PData >::operator[] ( const int& idx ) const { return levels[idx]; }

template< class PData >
template< class TData >
bool ProcessPartitionHierarchy< PData >::Initialize( const ProcessPartition< TData >& initialPartition , int height , int iters , int depths , bool repeat , int minSize , bool showPartition )
{
	ProcessPartitionHierarchy temp;
	if( repeat ) temp.resize( (depths+1)*2 );
	else		 temp.resize( (depths+1) );
	temp[0].resize( initialPartition.size() );
	if( showPartition )
	{
		IOServer::StdoutLock lock;
		printf( "Process Partition:\n" ) , fflush( stdout );
		for( int j=0 ; j<initialPartition.size() ; j++ )
		{
			printf( "[%6d" , initialPartition[j].width );
			int spaces = 6*1 + 2*(1-1) - 6;
			for( int k=0 ; k<spaces ; k++ ) printf( " " );
			printf( "]" );
		}
		printf( "\n" ) , fflush( stdout );
	}


	for( int j=0 ; j<initialPartition.size() ; j++ )
	{
		if( !ProcessPartition< PData >::IsValidBandSize( initialPartition[j].width , height , iters , minSize ) ) return false;
		temp[0][j].width = initialPartition[j].width;
	}
	int myDepth = depths;
	int depth = 0;
	for( int i=1 ; i<=myDepth ; i++ )
	{
		bool downSample = false;
		for( int j=0 ; j<temp[i-1].size() ; j++ ) if( !ProcessPartition< PData >::IsValidBandSize( temp[i-1][j].width>>1 , height>>(depth+1) , iters , minSize ) ) downSample = true;
		if( downSample )
		{
			if( repeat ) myDepth++;
			else		 depth++;

			if( temp[i-1].size()<2 ) fprintf( stderr , "[ERROR] Attempting to merge with only one thread\n" ) , exit(0);
			temp[i].resize( temp[i-1].size()>>1 );

			for( int j=0 ; j<temp[i  ].size() ; j++ ) temp[i][j].width =  0;
			for( int j=0 ; j<temp[i-1].size() ; j++ )
			{
				int jj = j>>1;
				if( jj>=temp[i].size() ) jj = temp[i].size()-1;
				temp[i][jj].width += temp[i-1][j].width;
				temp[i][jj].children.push_back( j );
			}
			for( int j=0 ; j<temp[i].size() ; j++ )
			{
				if( repeat )
				{
					if( !ProcessPartition< PData >::IsValidBandSize( temp[i][j].width>>1 , height>>(depth+1) , iters , minSize ) )
					{
						IOServer::printfID( "ProcessPartition::IsValidBandSize( %d , %d , %d , %d ) failed\n" , temp[i][j].width>>1 , height>>(depth+1) , iters , minSize ) , fflush( stdout );
						return false;
					}
				}
				else
				{
					temp[i][j].width >>= 1;
					if( !ProcessPartition< PData >::IsValidBandSize( temp[i][j].width , height>>(depth+1) , iters , minSize ) )
					{
						IOServer::printfID( "ProcessPartition::IsValidBandSize( %d , %d , %d , %d ) failed\n" , temp[i][j].width , height>>(depth+1) , iters , minSize ) , fflush( stdout );
						return false;
					}
				}
			}
		}
		else
		{
			depth++;
			temp[i].resize( temp[i-1].size() );
			for( int j=0 ; j<temp[i-1].size() ; j++ )
			{
				temp[i][j].width = temp[i-1][j].width>>1;
				temp[i][j].children.push_back( j );
			}
		}
		if( showPartition )
		{
			IOServer::StdoutLock lock;
			for( int j=0 ; j<temp[i].size() ; j++ )
			{
				printf( "[%6d" , temp[i][j].width );
				int spaces = 6*temp.leaves(i,j) + 2*(temp.leaves(i,j)-1) - 6;
				for( int k=0 ; k<spaces ; k++ ) printf( " " );
				printf( "]" );
			}
			printf("\n") , fflush( stdout );
		}
	}
	int count = 1;

	for( int i=1 ; i<=myDepth ; i++ ) if( temp[i].size() != temp[i-1].size() ) count++;
	levels.resize( count );
	levels[0].startDepth = 0;
	levels[0].resize( temp[0].size() );
	for( int j=0 ; j<levels[0].size() ; j++ )
	{
		levels[0][j].width = temp[0][j].width;
		levels[0][j].children.resize( temp[0][j].children.size() );
		for( size_t k=0 ; k<(*this)[0][j].children.size() ; k++ ) (*this)[0][j].children[k] = temp[0][j].children[k];
	}

	count = 1;
	for( int i=1 ; i<=myDepth ; i++ )
		if( temp[i].size() != temp[i-1].size() )
		{
			if( repeat )
			{
				(*this)[count-1].endDepth   = i-count;
				(*this)[count  ].startDepth = i-count;
			}
			else
			{
				(*this)[count-1].endDepth   = i-1;
				(*this)[count  ].startDepth = i;
			}
			(*this)[count].resize( temp[i].size() );
			for( int j=0 ; j<(*this)[count].size() ; j++ )
			{
				(*this)[count][j].width = temp[i][j].width;
				(*this)[count][j].children.resize( temp[i][j].children.size() );
				for( size_t k=0 ; k<(*this)[count][j].children.size() ; k++ ) (*this)[count][j].children[k] = temp[i][j].children[k];
			}
			count++;
		}
	(*this)[count-1].endDepth = depths;
	return true;
}
template< class PData >
int ProcessPartitionHierarchy< PData >::leaves( const int& depth , const int& offset ) const
{
	int sum = 0;
	if( !depth ) return 1;
	else for( size_t i=0 ; i<(*this)[depth][offset].children.size() ; i++ ) sum += leaves( depth-1 , (*this)[depth][offset].children[i] );
	return sum;
}

/////////////////////
// MultigridThread //
/////////////////////
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
MultigridThread< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::MultigridThread( void )
{
	lowPixels = pixels = labels = outGuess = NULL;
	_solverInfo = NullPointer< SolverInfo< PixelChannels > >( );
	_solvers = NullPointer< Pointer( SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType > ) >();
	_sRestriction = NULL;
	_sDivergence = NULL;
}
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
MultigridThread< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::~MultigridThread( void  )
{
	FreePointer( _solverInfo );
	if( _solvers )
	{
		for( int i=0 ; i<_blockCount ; i++ ) DeletePointer( _solvers[i] );
		DeletePointer( _solvers );
	}
	if( _sDivergence ) delete _sDivergence , _sDivergence = NULL;
	_sRestriction = NULL;
}
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
void MultigridThread< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::_init
(
	Pointer( ProcessingBlockData ) blockData , int blockCount , const GlobalData& globalData , bool showProgress , bool inCore , bool lowPixels , bool pixels , bool labels , bool outGuess
)
{
	_outOfCore  = !inCore;
	_blockData  = blockData;
	_blockCount = blockCount;
	_vCycles   = globalData.vCycles;
	_verbose   = globalData.verbose;
	_periodicType = globalData.periodicType;
	int depths = 0;
	for( int i=0 ; i<_blockCount ; i++ ) depths += _blockData[i].pData.depths();
	_solverInfo = AllocPointer< SolverInfo< PixelChannels > >( depths );
	_solvers = NewPointer< Pointer( SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType > ) >( _blockCount );
	// A block is a (nesting) multi-resolution sequence
	// The sub-blocks within a block are ordered from coarsest to finest
	for( int i=0 ; i<_blockCount ; i++ ) 
	{
		_solvers[i] = NewPointer< SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType > >( _blockData[i].pData.depths() );
		for( int j=1 ; j<_blockData[i].pData.depths()   ; j++ )	_solvers[i][j].parent = &_solvers[i][j-1];
		for( int j=0 ; j<_blockData[i].pData.depths()-1 ; j++ )	_solvers[i][j].rChild =  _solvers[i][j].pChild = &_solvers[i][j+1];
		for( int j=0 ; j<_blockData[i].pData.depths()   ; j++ )	_solvers[i][j].laneNum = globalData.lanes;
		_solvers[i][_blockData[i].pData.depths()-1].rChild = _solvers[i][_blockData[i].pData.depths()-1].pChild = NULL;
	}
	_solvers[0][_blockData[0].pData.depths()-1].showProgress = showProgress;

	if( pixels )
	{
		_sDivergence = new SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >();
		_sDivergence->unknownType = globalData.unknownType;
		_sDivergence->parent = &_solvers[0][_blockData[0].pData.depths()-1];
		_sRestriction = _sDivergence;
		_solvers[0][_blockData[0].pData.depths()-1].rChild = _sRestriction;
	}
}

template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
void MultigridThread< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::Initialize
(
	const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& gradientAverage ,
	MultiStreamIOServer* multiStreamIOServer ,
	Pointer( ProcessingBlockData ) blockData , int blockCount ,
	const GlobalData& globalData , bool showProgress , bool inCore
)
{
	double iWeight = globalData.iWeight;
	bool lump = globalData.lump;
	double gWeight = globalData.gWeight;
	double gScale = globalData.gScale;
	int iters = globalData.iters;
	_init( blockData , blockCount , globalData , showProgress , inCore , lowPixels!=NULL , pixels!=NULL , labels!=NULL , outGuess!=NULL );

	this->multiStreamIOServer = multiStreamIOServer;
	if( lowPixels ) lowPixels->SetServer( multiStreamIOServer ); // QUESTION
	if( pixels ) pixels->SetServer( multiStreamIOServer ); // QUESTION
	if( labels ) labels->SetServer( multiStreamIOServer ); // QUESTION
	if( outGuess ) outGuess->SetServer( multiStreamIOServer ); // QUESTION
	for( int i=0 ; i<_blockCount ; i++ )
	{
		DotProductStencil dotMajor , dotMinor , d2DotMajor , d2DotMinor;
		SetDownSampledStencil( globalData.width  , blockData[i].pData.startDepth , dotMajor , d2DotMajor , lump );
		SetDownSampledStencil( globalData.height , blockData[i].pData.startDepth , dotMinor , d2DotMinor , lump );
		if( pixels && !i )
			((SocketedStreamingDivergence< PixelChannels , LabelChannels , PixelType , LabelType , StorageType , SyncType >*)_sRestriction )->Initialize
				(
				lowPixels , pixels , labels , iWeight , lump , gWeight , gScale ,
				blockData[0].pData.start , blockData[0].pData.stop , blockData[0].pData.width , blockData[0].pData.height , iters ,
				blockData[0].leftStream , blockData[0].syncSockets , blockData[0].rightStream ,
				_outOfCore , _periodicType ,
				multiStreamIOServer , &gradientAverage
				);
		else
			_solvers[i][_blockData[i].pData.depths()-1].Initialize
				( dotMajor , d2DotMajor , dotMinor , d2DotMinor , iWeight , gWeight ,
				blockData[i].pData.start , blockData[i].pData.stop , blockData[i].pData.width , blockData[i].pData.height , iters ,
				blockData[i].leftStream , blockData[i].syncSockets , blockData[i].rightStream ,
				_outOfCore , _periodicType
				, multiStreamIOServer
				);
	}

	// Link the blocks together
	for( int i=0 ; i<_blockCount-1 ; i++ )
	{
		_solvers[i][0].parent = &_solvers[i+1][_blockData[i+1].pData.depths()-1];
		_solvers[i+1][_blockData[i+1].pData.depths()-1].rChild = _solvers[i+1][_blockData[i+1].pData.depths()-1].pChild = &_solvers[i][0];
	}
}
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int MultigridThread< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::RunThread( void* vparams )
{
#if MISHA_DENORMAL_CONTROL
	_MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
#endif // MISHA_DENORMAL_CONTROL
	MultigridThread* thread = ( MultigridThread* )vparams;
#if DEBUG_SOCKETS
	int testCount = 0;
	auto testCommunication = [&]( const char* message )
	{
		char str[512];
		printf( "Testing: (%s)\n" , message );
		sprintf( str , "test %d" , testCount );
		printf( "\tSending: %s\n" , str );
		SendOnSocket( thread->_blockData[0].serverSocket , (const char*)str , (size_t)512 , "Failed to send test string %d" , testCount );
		printf( "\tReceiving...\n" );
		ReceiveOnSocket( thread->_blockData[0].serverSocket , str , (size_t)512 , "Failed to receive test string %d" , testCount );
		printf( "\tReceived: %s\n" , str );
		testCount++;
	};
#endif // DEBUG_SOCKETS

	double t;

	MultiStreamIOClient* _X = NULL;
	MultiStreamIOClient* _B = NULL;

	// [Q] Don't we need to store the solution even if we don't have use _sRestriction?
	// [Q] Why does the usage of _X depend on whether on the block index is 0?
	// [A] Because block index 0 corresponds to the component of depths containing the highest resolution.
	//     In contrast, block thread->blockCount-1 is the component with the coarsest resolution.

	// If we are performing more than one V-cycle, allocate buffers for the constraints and solution
	if( thread->_vCycles>1 && thread->_sRestriction )
	{
		_X = new MultiStreamIOClient( (thread->_blockData->pData.stop-thread->_blockData->pData.start) * sizeof( SyncType )*PixelChannels , thread->_blockData->pData.height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
		_B = new MultiStreamIOClient( (thread->_blockData->pData.stop-thread->_blockData->pData.start) * sizeof( SyncType )*PixelChannels , thread->_blockData->pData.height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
	}

	SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* highSolver = &thread->_solvers[0][thread->_blockData[0].pData.depths()-1];
	SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* lowSolver  = &thread->_solvers[thread->_blockCount-1][0];

	for( int ii=0 ; ii<thread->_vCycles ; ii++ )
	{
#if DEBUG_SOCKETS
		testCommunication( "Starting V-Cycle" );
#endif // DEBUG_SOCKETS

		// On the first v-cycle, get the data from _sRestriction
		if( !ii ) highSolver->rChild = thread->_sRestriction;
		else      highSolver->rChild = NULL;

		/////////////////
		// RESTRICTION //
		/////////////////
		for( int b=0 ; b<thread->_blockCount ; b++ )
		{
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* solvers = thread->_solvers[b];
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >& lSolver = solvers[0];
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >& hSolver = solvers[thread->_blockData[b].pData.depths()-1];

			// Initialize the norms and averages to 0
			for( int i=0 ; i<thread->_blockData[b].pData.depths() ; i++ )
			{
				solvers[i].bSquareNorm = solvers[i].rSquareNorm = solvers[i].xSquareNorm = 0;
				for ( int c = 0 ; c < PixelChannels ; c++ )	solvers[i].solutionSum[c] = 0;
				solvers[i].setResidual = true;
			}

			// Clear everything
			lSolver.inX = lSolver.outX = hSolver.inX = hSolver.outX = NULL;
			lSolver.inB = lSolver.outB = hSolver.inB = hSolver.outB = NULL;
			lSolver.outR = hSolver.outR = NULL;

			// If there are multiple blocks the constraints for the finest levels of the coarser blocks are streamed in
			// and the restricted residual from the coarsest level is streamed out
			hSolver.inB = thread->_blockData[b].inHighB;
			lSolver.outR = thread->_blockData[b].outLowR;
		}

		// At the finest resolution the constraints/solution are either:
		// -- obtained from _sRestriction/initial-guess on the first pass
		// -- buffered into temporary storage for subsequent passes
		if( thread->_sRestriction )
		{
			if( ii ) highSolver->inX = _X , highSolver->inB = _B;
			else     highSolver->inX = thread->outGuess , highSolver->outB = _B;
		}

		// Data for the interleaved streaming multigrid
		t=Time();
		// Initialize
		if( ii || !thread->_sRestriction ) highSolver->InitRestriction() , highSolver->SetRestriction();
		else thread->_sRestriction->InitRestriction() , thread->_sRestriction->SetRestriction();
		// Solve
		// [QUESTION] Why do I have to comment this out?
//		if( ii ) in->SetServer( &SocketedStreamingSolver< Channels >::server );
		if( ii || !thread->_sRestriction ) highSolver->SolveRestriction();
		else							   thread->_sRestriction->SolveRestriction();

		t = Time()-t;
		int idx = 0;
		for( int b=0 ; b<thread->_blockCount ; b++ )
		{
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* solvers = thread->_solvers[b];
			for( int i=0 ; i<thread->_blockData[b].pData.depths() ; i++ )
			{
				thread->_solverInfo[idx].bSquareNorm = solvers[i].bSquareNorm;
				thread->_solverInfo[idx].rSquareNorm = solvers[i].rSquareNorm;
				thread->_solverInfo[idx].xSquareNorm = solvers[i].xSquareNorm;
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					thread->_solverInfo[idx].solutionSum[c]  = solvers[i].solutionSum[c];
					thread->_solverInfo[idx].solutionSum[c] /= solvers[i].major;
					thread->_solverInfo[idx].solutionSum[c] /= solvers[i].minor;
				}
				idx++;
			}
		}
		for( int b=0 ; b<thread->_blockCount ; b++ )
		{
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* solvers = thread->_solvers[b];
			char id[512];
			SetThisThreadID( id );
			{
				IOServer::StdoutLock lock;
				if( thread->_verbose )
				{
					printf( "Thread Restriction [%s]:    %f\n" , id , t ) , fflush( stdout );
					for( int i=thread->_blockData[b].pData.depths()-1 ; i>=0 ; i-- )
						printf( "\tError[%d x %d] %g -> %g\n" , solvers[i].size() , solvers[i].minor , sqrt( solvers[i].bSquareNorm ) , sqrt( solvers[i].rSquareNorm ) ) , fflush( stdout );
				}
			}
		}

#if DEBUG_SOCKETS
		testCommunication( "Sending solver info" );
#endif // DEBUG_SOCKETS
		Pointer( SolverInfo< PixelChannels > ) solverInfo = thread->_solverInfo;
		for( int b=0 ; b<thread->_blockCount ; b++ )
		{
			SendOnSocket( thread->_blockData[b].serverSocket , ( ConstPointer( SolverInfo< PixelChannels > ) )solverInfo , sizeof( SolverInfo<PixelChannels> )*thread->_blockData[b].pData.depths() , "Failed so send restriction info to server" );
			solverInfo += thread->_blockData[b].pData.depths();
		}
#if DEBUG_SOCKETS
		testCommunication( "Sending average" );
#endif // DEBUG_SOCKETS

		if( !ii && thread->_sRestriction )
		{
			SendOnSocket( thread->_blockData[0].serverSocket , ( ConstPointer( AverageColor< PixelChannels > ) )GetPointer( thread->_sDivergence->average ) , sizeof(thread->_sDivergence->average) , "Failed to send average to server" );
			thread->_sDivergence->UnSetRestriction();
		}
		else highSolver->UnSetRestriction();

#if DEBUG_SOCKETS
		testCommunication( "Finished restriction" );
#endif // DEBUG_SOCKETS

		//////////////////
		// PROLONGATION //
		//////////////////
#if DEBUG_SOCKETS
		testCommunication( "Starting prolongation" );
#endif // DEBUG_SOCKETS
		for( int b=0 ; b<thread->_blockCount ; b++ )
		{
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* solvers = thread->_solvers[b];
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >& lSolver = solvers[0];
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >& hSolver = solvers[thread->_blockData[b].pData.depths()-1];

			for( int i=0 ; i<thread->_blockData[b].pData.depths() ; i++ )
			{
				solvers[i].bSquareNorm = solvers[i].rSquareNorm = solvers[i].xSquareNorm = 0;
				for ( int c=0 ; c<PixelChannels ; c++ )	solvers[i].solutionSum[c] = 0;
				solvers[i].setResidual = thread->_verbose;
			}
			// Clear everything
			lSolver.inX = lSolver.outX = hSolver.inX = hSolver.outX = NULL;
			lSolver.inB = lSolver.outB = hSolver.inB = hSolver.outB = NULL;
			lSolver.outP = hSolver.outP = NULL;

			// If there are multiple blocks the solution for the coarsest level is streamed in
			// and the prolonged solution from the finest level is streamed out
			lSolver.inX = thread->_blockData[b].inLowX;
			lSolver.outP = thread->_blockData[b].outHighP;
		}
		if( thread->_sRestriction )
		{
			if( ii<thread->_vCycles-1 ) highSolver->outX = _X;
			else
			{
				highSolver->outX = thread->_blockData[0].outHighX;
				if( _X ) delete _X , _X = NULL;
				if( _B ) delete _B , _B = NULL;
			}
		}

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		highSolver->InitProlongation();
		lowSolver->SetProlongation();
		// Solve
		// [BADNESS] Why do I have to comment this out?
//		if( ii<_vCycles-1 ) _X->SetServer( &StreamingSolver< Real , Type , Degree , Channels >::server );
		lowSolver->SolveProlongation();

		t = Time() - t;
		idx = 0;
		for( int b=0 ; b<thread->_blockCount ; b++ )
		{
			SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* solvers = thread->_solvers[b];
			for( int i=0 ; i<thread->_blockData[b].pData.depths() ; i++ )
			{
				thread->_solverInfo[idx].bSquareNorm = solvers[i].bSquareNorm;
				thread->_solverInfo[idx].rSquareNorm = solvers[i].rSquareNorm;
				thread->_solverInfo[idx].xSquareNorm = solvers[i].xSquareNorm;
				for( int c=0 ; c<PixelChannels ; c++ )
				{
					thread->_solverInfo[idx].solutionSum[c]  = solvers[i].solutionSum[c];
					thread->_solverInfo[idx].solutionSum[c] /= solvers[i].major;
					thread->_solverInfo[idx].solutionSum[c] /= solvers[i].minor;
				}
				idx++;
			}
		}
		{
			char id[512];
			SetThisThreadID( id );
			{
				IOServer::StdoutLock lock;
				if( thread->_verbose )
				{
					printf( "Thread Prolongation [%s]:    %f\n" , id , t ) , fflush( stdout );
					for( int b=thread->_blockCount-1 ; b>=0 ; b-- )
					{
						SocketedMultiGridStreamingSolver< PixelChannels , StorageType , SyncType >* solvers = thread->_solvers[b];
						for( int i=0 ; i<thread->_blockData[b].pData.depths() ; i++ )
							printf( "\tError[%d x %d] %g -> %g\n" , solvers[i].size() , solvers[i].minor , sqrt( solvers[i].bSquareNorm ) , sqrt( solvers[i].rSquareNorm ) ) , fflush( stdout );
					}
				}
			}
		}
#if DEBUG_SOCKETS
		testCommunication( "Sending solver info" );
#endif // DEBUG_SOCKETS
		solverInfo = thread->_solverInfo;
		for( int b=0 ; b<thread->_blockCount ; b++ )
		{
			SendOnSocket ( thread->_blockData[b].serverSocket , ( ConstPointer( SolverInfo< PixelChannels > ) )solverInfo , sizeof( SolverInfo< PixelChannels > )*thread->_blockData[b].pData.depths() , "Failed so send restriction info to server" );
			solverInfo += thread->_blockData[b].pData.depths();
		}
		lowSolver->UnSetProlongation();
#if DEBUG_SOCKETS
		testCommunication( "Done prolongation" );
#endif // DEBUG_SOCKETS
	}
	return 0;
}
