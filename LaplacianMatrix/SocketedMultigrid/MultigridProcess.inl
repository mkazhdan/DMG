#include "Util/XPlatform.h"
#include "Util/MultiStreamIO.h"
#include "Util/Socket.h"
#include "Util/Half/half.h"
#include "Util/ChannelConverter.h"

#define SOCKET_CONNECT_WAIT 500

/////////////////////////////
// SocketedMultigridClient //
/////////////////////////////
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::SocketedMultigridClient(void)
{
	_threadBlockData = NullPointer< ProcessingBlockData >();
	_pConnectInfo = NullPointer< ProcessConnectionInfo >();
	_processCount = 0;
}
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::~SocketedMultigridClient(void)
{
	if( _threadBlockData )
	{
		for( int i=0 ; i<_processCount ; i++ )
		{
			if( _threadBlockData[i].outHighP ) delete _threadBlockData[i].outHighP , _threadBlockData[i].outHighP = NULL;
			if( _threadBlockData[i].inHighB )  delete _threadBlockData[i].inHighB  , _threadBlockData[i].inHighB  = NULL;
			if( _threadBlockData[i].outLowR )  delete _threadBlockData[i].outLowR  , _threadBlockData[i].outLowR  = NULL;
			if( _threadBlockData[i].inLowX )   delete _threadBlockData[i].inLowX   , _threadBlockData[i].inLowX   = NULL;
		}
		DeletePointer( _threadBlockData );
	}
	DeletePointer( _pConnectInfo );
	_processCount = 0;
}
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
bool SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::shortSync( void ) const { return _globalData.shortSync; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
bool SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::gammaCorrection( void ) const { return _globalData.gammaCorrection; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
bool SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::hdr( void ) const { return _globalData.hdr; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::tileWidth( void ) const { return _globalData.tileWidth; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::tileHeight( void ) const { return _globalData.tileHeight; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
const char* SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::tileExtension( void ) const { return _globalData.tileExt; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::quality( void ) const { return _globalData.quality; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::cSize( void ) const { return _clientData.cEnd - _clientData.start; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::size( void ) const { return _clientData.end - _clientData.start; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::cHeight( void ) const { return _globalData.cHeight; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::height( void ) const { return _globalData.height; }
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
int SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::processCount( void ) const { return _processCount; }

template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
bool SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::SetUp( char* serverAddress , char* prefix , int port , int index , int w , int h , const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& averageMap , bool inCore )
{
	int pid = GetThisProcessID( );
	Socket serverSocket = GetConnectSocket( serverAddress , port , SOCKET_CONNECT_WAIT , false );

	bool haveAddress = false;
	if( prefix )
	{
		if( GetHostEndpointAddress( &_myAddress , prefix ) ) haveAddress = true;
		else fprintf( stderr , "[WARNING]: Failed to find address starting with: %s\n" , prefix );
	}
	if( !haveAddress ) _myAddress = GetLocalSocketEndpointAddress( serverSocket );
	// Send the band information
	_clientData.index  = index;
	_clientData.width  = w;
	_clientData.height = h;
	_inCore = inCore;
	SendOnSocket( serverSocket , ( ConstPointer( ClientData ) )GetPointer( _clientData ) , sizeof( _clientData ) , "Failed to send client data to server" );
	// Send average gradient values here
	{
		int count = (int)( averageMap.size() );
		SendOnSocket( serverSocket , ( ConstPointer( int ) )GetPointer( count ) , sizeof( count ) , "Failed to send gradient count to server" );
		long long size = sizeof( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) * count;
		if( count ) SendOnSocket( serverSocket , ( ConstPointer( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) )GetPointer( averageMap ) , (int)( size ) , "Failed to send average gradients to server" );
	}
	// Now receive back the average gradient values
	{
		int count;
		ReceiveOnSocket( serverSocket , GetPointer( count ) , sizeof( count ) , "Failed to receive gradient count from server" );
		long long size = sizeof( std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > ) * count;
		_averageMap.resize( count );
		if( count ) ReceiveOnSocket( serverSocket , GetPointer( _averageMap ) , (int)( size ) , "Failed to receive average gradients from server" );
	}

	// Get the global information
	ReceiveOnSocket( serverSocket , GetPointer( _globalData ) , sizeof( _globalData ) , "Failed to receive global data from server ");

	// Get the local information
	ReceiveOnSocket( serverSocket , GetPointer( _clientData ) , sizeof( _clientData ) , "Failed to receive client data from server" );

	// Get back the total number of threads this client is responsible for and the information about these threads
	ReceiveOnSocket( serverSocket , GetPointer( _processCount ) , sizeof( _processCount ) , "Failed to receive process count from server" );
	Pointer( ProcessData ) processData = AllocPointer< ProcessData >( _processCount );
	ReceiveOnSocket( serverSocket , processData , sizeof( ProcessData ) * _processCount , "Failed to receive process data from server" );

	// At the point the client should not have to communicate with the server any more
	CloseSocket( serverSocket );

	// Start constructing the threads that will run the different processes
	_threadBlockData = NewPointer< ProcessingBlockData >( _processCount );
	Pointer( IPData ) addressX    = NewPointer< IPData >( _processCount );
	Pointer( IPData ) addressB    = NewPointer< IPData >( _processCount );
	Pointer( IPData ) addressLeft = NewPointer< IPData >( _processCount );
	for( int i=0 ; i<_processCount ; i++ )
	{
		_threadBlockData[i].pData = processData[i];
		memcpy( &addressX[i].address    , &_myAddress , sizeof( _myAddress ) );
		memcpy( &addressB[i].address    , &_myAddress , sizeof( _myAddress ) );
		memcpy( &addressLeft[i].address , &_myAddress , sizeof( _myAddress ) );
	}
	FreePointer( processData );


	_pConnectInfo = NewPointer< ProcessConnectionInfo >( _processCount );
	int listenerCount = 0;
	Pointer( ThreadHandle ) listenerHandles = NewPointer< ThreadHandle >( _processCount*3 );
	Pointer( Socket ) leftSockets = NewPointer< Socket >( _processCount );


	// Establish the connection between each process and the server
	for( int i=0 ; i<_processCount ; i++ )
	{
		_threadBlockData[i].serverSocket = GetConnectSocket( serverAddress , port , SOCKET_CONNECT_WAIT , false );
		SendOnSocket( _threadBlockData[i].serverSocket , ( ConstPointer( int ) )GetPointer( _threadBlockData[i].pData.index  ) , sizeof( _threadBlockData[i].pData.index ) , "Failed to send thread index to server" );
		SendOnSocket( _threadBlockData[i].serverSocket , ( ConstPointer( int ) )GetPointer( _threadBlockData[i].pData.offset ) , sizeof( _threadBlockData[i].pData.offset ) , "Failed to send thread offset to server" );
		if( _threadBlockData[i].pData.index )	// If not at the leaf
		{
			_pConnectInfo[i].init( _threadBlockData[i].pData.children );
			_pConnectInfo[i].listenChildX = GetListenSocket( addressX[i].port );
			if ( _pConnectInfo[i].listenChildX == _INVALID_ACCEPTOR_SOCKET_ ) return false;
			SendOnSocket ( _threadBlockData[i].serverSocket , ( ConstPointer( IPData ) )GetPointer( addressX[i] ) , sizeof( addressX[i] ) , "Failed to send X address to server" );
			listenerHandles[listenerCount] = SpawnListenerSocket( _pConnectInfo[i].listenChildX , _pConnectInfo[i].childX , _threadBlockData[i].pData.children );
			if( !TestThreadHandle( listenerHandles[listenerCount] ) ) return false;
			listenerCount++;

			_pConnectInfo[i].listenChildB = GetListenSocket( addressB[i].port );
			if ( _pConnectInfo[i].listenChildB == _INVALID_ACCEPTOR_SOCKET_ ) return false;
			SendOnSocket ( _threadBlockData[i].serverSocket , ( ConstPointer( IPData ) )GetPointer( addressB[i] ) , sizeof( addressB[i] ) , "Failed to send B address to server" );
			listenerHandles[listenerCount] = SpawnListenerSocket( _pConnectInfo[i].listenChildB , _pConnectInfo[i].childB , _threadBlockData[i].pData.children );
			if( !TestThreadHandle( listenerHandles[listenerCount] ) ) return false;
			listenerCount++;
		}
		else _pConnectInfo[i].listenChildX = _pConnectInfo[i].listenChildB = _INVALID_ACCEPTOR_SOCKET_;

		if( _globalData.periodicType!=NO_PERIODIC || _threadBlockData[i].pData.offset ) // If spherical or not left
		{
			_pConnectInfo[i].listenLeft = GetListenSocket( addressLeft[i].port );
			if ( _pConnectInfo[i].listenLeft == _INVALID_ACCEPTOR_SOCKET_ ) return false;
			SendOnSocket ( _threadBlockData[i].serverSocket , ( ConstPointer( IPData ) )GetPointer( addressLeft[i] ) , sizeof( addressLeft[i] ) , "Failed to send left address to server" );
			listenerHandles[listenerCount] = SpawnListenerSocket( _pConnectInfo[i].listenLeft , &leftSockets[i] , 1 );
			if( !TestThreadHandle( listenerHandles[listenerCount] ) ) return false;
			listenerCount++;
		}
		else leftSockets[i] = _INVALID_SOCKET_;
	}
	
	// Get connection information from the server and make the connection
	for( int i=0 ; i<_processCount ; i++ )
	{
		IPData neighborData;
		ReceiveOnSocket( _threadBlockData[i].serverSocket , GetPointer( neighborData ) , sizeof( neighborData ) , "Failed to get parent X address from server" );
		_pConnectInfo[i].parentX = GetConnectSocket( neighborData.address , neighborData.port , SOCKET_CONNECT_WAIT , false );
		ReceiveOnSocket( _threadBlockData[i].serverSocket , GetPointer( neighborData ) , sizeof( neighborData ) , "Failed to receive parent B address from server" );
		_pConnectInfo[i].parentB = GetConnectSocket( neighborData.address , neighborData.port , SOCKET_CONNECT_WAIT , false );

		if( _globalData.periodicType!=NO_PERIODIC || _threadBlockData[i].pData.offset < _threadBlockData[i].pData.maxOffset )
		{
			ReceiveOnSocket( _threadBlockData[i].serverSocket , GetPointer( neighborData ) , sizeof( IPData ) , "Failed to get right address from server" );
			_pConnectInfo[i].right = new SocketStream( neighborData.address , neighborData.port , SOCKET_CONNECT_WAIT , false );
		}
		else _pConnectInfo[i].right = NULL;
	}
	if( listenerCount ) WaitOnThreads( listenerHandles , listenerCount , 1000 , "SocketedMultigridClient::SetUp" );
	for( int i=0 ; i<_processCount ; i++ ) if( leftSockets[i]!=_INVALID_SOCKET_ ) _pConnectInfo[i].left = new SocketStream( leftSockets[i] );

	if( _globalData.periodicType==SPHERICAL_PERIODIC )
	{
		{
			int i=0;
			int depths = _threadBlockData[i].pData.endDepth - _threadBlockData[i].pData.startDepth + 1;
			Pointer( int ) ports = AllocPointer< int >( depths+1 );
			_threadBlockData[0].syncSockets = NewPointer< Socket >( depths+1 );
			ReceiveOnSocket ( _threadBlockData[i].serverSocket , ports , sizeof(int) , "Failed to get port from server" );
			ReceiveOnSocket ( _threadBlockData[i].serverSocket , ports+1 , sizeof(int)*depths , "Failed to get ports from server" );

			for( int j=0 ; j<depths+1 ; j++ )
			{
				_threadBlockData[i].syncSockets[j] = GetConnectSocket( serverAddress , ports[j] , SOCKET_CONNECT_WAIT , false );
				if( _threadBlockData[i].syncSockets[j] == _INVALID_SOCKET_ )	return false;
			}
			FreePointer( ports );
		}
		for( int i=1 ; i<_processCount ; i++ )
		{
			int depths = _threadBlockData[i].pData.endDepth - _threadBlockData[i].pData.startDepth + 1;
			Pointer( int ) ports = AllocPointer< int >( depths );
			_threadBlockData[i].syncSockets = NewPointer< Socket >( depths );
			ReceiveOnSocket ( _threadBlockData[i].serverSocket , ports , sizeof(int)*depths , "Failed to get ports from server" );
			for( int j=0 ; j<depths ; j++ )
			{
				_threadBlockData[i].syncSockets[j] = GetConnectSocket( serverAddress , ports[j] , SOCKET_CONNECT_WAIT , false );
				if( _threadBlockData[i].syncSockets[j]==_INVALID_SOCKET_ )	return false;
			}
			FreePointer( ports );
		}
	}
	else for( int i=0 ; i<_processCount ; i++ ) _threadBlockData[i].syncSockets = NullPointer< Socket >();
	DeletePointer( listenerHandles );
	DeletePointer( addressX );
	DeletePointer( addressB );
	DeletePointer( addressLeft );
	DeletePointer( leftSockets );
	return true;
}
template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
void SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType >::Run( StreamingGrid* lowPixels , StreamingGrid* pixels , StreamingGrid* labels , StreamingGrid* outGuess , StreamingGrid* out , MultiStreamIOServer* multiStreamIOServer , bool showProgress )
{
	// Initialize the processes
	MultigridThread< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType > thread;
	thread.lowPixels = lowPixels;
	thread.pixels = pixels;
	thread.labels = labels;
	thread.outGuess = outGuess;
	_threadBlockData[0].leftStream  = _pConnectInfo[0].left;
	_threadBlockData[0].rightStream = _pConnectInfo[0].right;
	_threadBlockData[0].outHighX = out;
	_threadBlockData[0].outHighP = NULL;

	for( int i=1 ; i<_processCount ; i++ )
	{
		_threadBlockData[i].leftStream  = _pConnectInfo[i].left;
		_threadBlockData[i].rightStream = _pConnectInfo[i].right;
	}
	int pid = GetThisProcessID( );
	for( int i=0 ; i<_processCount ; i++ )
		if( _threadBlockData[i].rightStream )
		{
			_threadBlockData[i].rightStream->write( ( Pointer( byte ) )GetPointer( _myAddress ) , sizeof( _myAddress ) );
			_threadBlockData[i].rightStream->write( ( Pointer( byte ) )GetPointer( pid ) , sizeof( pid ) );
		}
	for( int i=0 ; i<_processCount ; i++ )
		if( _threadBlockData[i].leftStream )
		{
			EndpointAddress leftAddress;
			int leftPID;
			_threadBlockData[i].leftStream->read( ( Pointer( byte ) )GetPointer( leftAddress ) , sizeof( leftAddress ) );
			_threadBlockData[i].leftStream->read( ( Pointer( byte ) )GetPointer( leftPID ) , sizeof( leftPID ) );

			if( AddressesEqual( _myAddress , leftAddress ) && pid==leftPID )
			{
				SharedMemoryBuffer::StreamPair sPair;
				SharedMemoryBuffer::StreamPair::CreateSharedBufferPair( sPair );
				_threadBlockData[i].leftStream->write( ( Pointer( byte ) )GetPointer( sPair.second ) , sizeof( sPair.second ) );
				delete _threadBlockData[i].leftStream;
				_pConnectInfo[i].left = sPair.first;
				_threadBlockData[i].leftStream = sPair.first;
			}
			else
			{
				SharedMemoryBuffer::SecondStream* sStream = NULL;
				_threadBlockData[i].leftStream->write( ( Pointer( byte ) )GetPointer( sStream ) , sizeof( sStream ) );
			}
		}
	for( int i=0 ; i<_processCount ; i++ )
		if( _threadBlockData[i].rightStream )
		{
			SharedMemoryBuffer::SecondStream* sStream;
			_threadBlockData[i].rightStream->read( ( Pointer( byte ) )GetPointer( sStream ) , sizeof( sStream ) );
			if( sStream )
			{
				delete _threadBlockData[i].rightStream;
				_pConnectInfo[i].right = sStream;
				_threadBlockData[i].rightStream = sStream;
			}
		}

	for( int i=0 ; i<_processCount ; i++ )
	{
		int width  = _threadBlockData[i].pData.stop-_threadBlockData[i].pData.start;
		int height = _threadBlockData[i].pData.height;
		width  >>= _threadBlockData[i].pData.endDepth-_threadBlockData[i].pData.startDepth;
		height >>= _threadBlockData[i].pData.endDepth-_threadBlockData[i].pData.startDepth;

		SendOnSocket( _pConnectInfo[i].parentB , ( ConstPointer( int ) )GetPointer( _threadBlockData[i].pData.offset ) , sizeof( _threadBlockData[i].pData.offset ) , "Failed to send offset to parent B" );
		int width2 = width/2;
		SendOnSocket( _pConnectInfo[i].parentB , ( ConstPointer( int ) )GetPointer( width2 ) , sizeof( int ) , "failed to send width to parent B" );
		SendOnSocket( _pConnectInfo[i].parentX , ( ConstPointer( int ) )GetPointer( _threadBlockData[i].pData.offset ) , sizeof( _threadBlockData[i].pData.offset ) , "Failed to send offset to parent X" );
		SendOnSocket( _pConnectInfo[i].parentX , ( ConstPointer( int ) )GetPointer( width ) , sizeof( int ) , "Failed to send width to parent X" );
		_threadBlockData[i].outLowR = new SocketBackedGrid( _pConnectInfo[i].parentB , width * PixelChannels * sizeof( SyncType ) /2 , height/2 );
		_threadBlockData[i].inLowX  = new SocketBackedGrid( _pConnectInfo[i].parentX , width * PixelChannels * sizeof( SyncType )    , height   );
	}
	for( int i=1 ; i<_processCount ; i++ )
	{
		_threadBlockData[i].inHighB  = GetMultiSocketBackedGrid< SyncType , PixelChannels > ( _pConnectInfo[i].childB , _threadBlockData[i].pData.children , _threadBlockData[i].pData.height    );
		_threadBlockData[i].outHighP = GetMultiSocketBackedGrid< SyncType , PixelChannels > ( _pConnectInfo[i].childX , _threadBlockData[i].pData.children , _threadBlockData[i].pData.height<<1 );
	}
	thread.Initialize( _averageMap , multiStreamIOServer , _threadBlockData , _processCount , _globalData , (_globalData.showProgress && showProgress) , _inCore );
	thread.RunThread( &thread );
	for( int i=0; i<_processCount ; i++ )
	{
		if( _pConnectInfo[i].left ) delete _pConnectInfo[i].left;
		if( _pConnectInfo[i].right ) delete _pConnectInfo[i].right;
		CloseSocket( _threadBlockData[i].serverSocket );
		CloseSocket( _pConnectInfo[i].parentX );
		CloseSocket( _pConnectInfo[i].parentB );
		if( _pConnectInfo[i].childX ) for( int j=0 ; j<_pConnectInfo[i].children ; j++ ) CloseSocket( _pConnectInfo[i].childX[j] );
		if( _pConnectInfo[i].childB ) for( int j=0 ; j<_pConnectInfo[i].children ; j++ ) CloseSocket( _pConnectInfo[i].childB[j] );


		if( _globalData.periodicType==SPHERICAL_PERIODIC )
			for( int i=0 ; i<_processCount ; i++ )
			{
				int depths = _threadBlockData[i].pData.endDepth - _threadBlockData[i].pData.startDepth + 1;
				if( !i ) for( int j=0 ; j<depths+1 ; j++ ) CloseSocket( _threadBlockData[i].syncSockets[j] );
				else	 for( int j=0 ; j<depths   ; j++ ) CloseSocket( _threadBlockData[i].syncSockets[j] );
			}
	}
}
//////////////////////////////////
// SocketedSuperMultigridClient //
//////////////////////////////////
template< class StorageType , class PixelType , class LabelType >
bool SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::SetUp( char* address , char* prefix , int port , int index , int w , int h , int subClients , bool inCore , MultiStreamIOServer* server , bool& getGradientAverage , bool& useGrayImage )
{
	_address = address;
	_prefix = prefix;
	_port = port;
	_clientNum = subClients;
	_inCore = inCore;

	sharedIOServer = server;
	printf( "Connecting to server..." ) , fflush( stdout );
	Socket serverSocket = GetConnectSocket( address , port , SOCKET_CONNECT_WAIT , true );

	// Send the band information
	_clientData.index	= index;
	_clientData.width	= w;
	_clientData.height	= h;
	_clientData.subClients = subClients;
	SendOnSocket( serverSocket , ( ConstPointer( ClientData ) )GetPointer( _clientData ) , sizeof( _clientData ) , "Failed to send client data to server" );

	// Get the local information
	ReceiveOnSocket( serverSocket , GetPointer( _globalData ) , sizeof( _globalData ) , "Failed to get global data from server" );
	ReceiveOnSocket( serverSocket , GetPointer( _clientData ) , sizeof( _clientData ) , "Failed to get client data from server" );
	ReceiveOnSocket( serverSocket , GetPointer( _blockSize  ) , sizeof( _blockSize )  , "Failed to get block size from server" );
	int _getGradientAverage;
	ReceiveOnSocket( serverSocket , GetPointer( _getGradientAverage ) , sizeof( _getGradientAverage ) , "Failed to get gradient average status from server" );
	getGradientAverage = (_getGradientAverage!=0);
	int _useGrayImage;
	ReceiveOnSocket( serverSocket , GetPointer( _useGrayImage ) , sizeof( _useGrayImage ) , "Failed to get use gray status from server" );
	useGrayImage = (_useGrayImage!=0);
	CloseSocket( serverSocket );

	if( _globalData.verbose )
	{
		printf( "\tImage Size: %d x %d\n" , w , h ) , fflush( stdout );
		printf( "\t   Address: %s:%d\n" , address , port ) , fflush( stdout );
	}
	if( (size()/_blockSize)*_blockSize != size() )
		fprintfId( stderr , "Band width is not a multiple of the block size: %d , %d\n", size() , _blockSize );
	return true;
}
template< class StorageType , class PixelType , class LabelType >
template< int PixelChannels , int LabelChannels >
void SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::Run( StreamingGrid* lowPixels , StreamingGrid* pixels , StreamingGrid* labels , StreamingGrid* outGuess , StreamingGrid* out , const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& averageMap )
{
	if( _inCore )
	{
		MemoryBackedGrid *_lowPixels = NULL , *_pixels = NULL , *_labels = NULL , *_outGuess = NULL , *_out = NULL;
		_pixels = new MemoryBackedGrid( pixels->rowSize() , pixels->rows() );
		if( lowPixels ) _lowPixels = new MemoryBackedGrid( lowPixels->rowSize() , lowPixels->rows() );
		if( labels ) _labels = new MemoryBackedGrid( labels->rowSize() , labels->rows() );
		if( outGuess ) _outGuess = new MemoryBackedGrid( outGuess->rowSize() , outGuess->rows() );
		_out = new MemoryBackedGrid( out->rowSize() , out->rows() );
		for( int i=0 ; i<pixels->rows() ; i++ ) memcpy( (*_pixels)[i] , (*pixels)[i] , pixels->rowSize() ) , pixels->advance();
		if( lowPixels ) for( int i=0 ; i<lowPixels->rows() ; i++ ) memcpy( (*_lowPixels)[i] , (*lowPixels)[i] , lowPixels->rowSize() ) , lowPixels->advance();
		if( labels ) for( int i=0 ; i<labels->rows() ; i++ ) memcpy( (*_labels)[i] , (*labels)[i] , labels->rowSize() ) , labels->advance();
		if( outGuess ) for( int i=0 ; i<outGuess->rows() ; i++ ) memcpy( (*_outGuess)[i] , (*outGuess)[i] , outGuess->rowSize() ) , outGuess->advance();
		if( _globalData.shortSync ) _Run< PixelChannels , LabelChannels , half  >( _lowPixels , _pixels , _labels , _outGuess , _out , averageMap );
		else                        _Run< PixelChannels , LabelChannels , float >( _lowPixels , _pixels , _labels , _outGuess , _out , averageMap );
		for( int i=0 ; i<out->rows() ; i++ ) memcpy( (*out)[i] , (*_out)[i] , out->rowSize() ) , out->advance();
		delete _pixels;
		if( _lowPixels ) delete _lowPixels;
		if( _labels ) delete _labels;
		if( _outGuess ) delete _outGuess;
		delete _out;
	}
	else
		if( _globalData.shortSync ) _Run< PixelChannels , LabelChannels , half  >( lowPixels , pixels , labels , outGuess , out , averageMap );
		else                        _Run< PixelChannels , LabelChannels , float >( lowPixels , pixels , labels , outGuess , out , averageMap );
}
template< class StorageType , class PixelType , class LabelType >
template< int PixelChannels , int LabelChannels , class SyncType > 
void SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::_Run( StreamingGrid* lowPixels , StreamingGrid* pixels , StreamingGrid* labels , StreamingGrid* outGuess , StreamingGrid* out , const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& averageMap )
{
	typedef _SocketedMultigridClientData< PixelChannels , LabelChannels , SyncType > SocketedMultigridClientData;
	Pointer( ThreadHandle ) clientHandles = NewPointer< ThreadHandle >( _clientNum  );
	Pointer( SocketedMultigridClientData ) clientData = NewPointer< SocketedMultigridClientData >( _clientNum );

#if 1 // Shouldn't make a difference since size() should be a multiple of _blockSize
	int blockNum = size() / _blockSize;
#else
	int blockNum = (size() + _blockSize - 1 ) / _blockSize;
#endif
	// _blockSize specifies the level of atomicity required by the solver
	// Error with:
	// clientData[_clientNum-1].width = 1280
	// size()					      = 4608
	// cSize()						  = 2800
	// _blockSize					  = 256
	// blockNum						  = 18
	// (Assuming) _clientNum		  = 4
	// clientData[0].width			  = (18/4) * _blockSize			= 1024	->	[   0 , 1024)
	// clientData[1].width			  = (36/4 - 18/4) * _blockSize	= 1280	->	[1024 , 2304)
	// clientData[2].width			  = (54/4 - 36/4) * _blockSize	= 1024	->	[2304 , 3328)
	// clientData[3].width			  = (72/4 - 54/4) * _blockSize	= 1280	->	[3328 , 4608)

	int widthSum = 0;
	for( int i=0 ; i<_clientNum ; i++ )
	{
		clientData[i].address = _address;
		clientData[i].prefix = _prefix;
		clientData[i].port = _port;
		clientData[i].index = _clientData.index+i;
		clientData[i].height = height();
		clientData[i].cHeight = cHeight();
		clientData[i].inCore = _inCore;
		clientData[i].cWidth = clientData[i].width = ( ( ( blockNum * (i+1) ) / _clientNum ) - ( ( blockNum * i ) / _clientNum ) ) * _blockSize;
		widthSum += clientData[i].width;
		if( widthSum>cSize() )
		{
			clientData[i].cWidth = clientData[i].width - ( widthSum - cSize() );
			if( clientData[i].cWidth <= 0 )
			{
				fprintf( stderr , "Warning effective client width set to zero!\n" );
				clientData[i].cWidth = 0;
			}

		}
		clientData[i].averageMapPtr = &averageMap;
		clientHandles[i] = clientData[i].spawnSetUp();
	}
	WaitOnThreads( clientHandles , _clientNum , 10000 , "SocketedSuperMultigridClient::_Run (Set-Up)" );
	int pCount = 0;
	for( int i=0 ; i<_clientNum ; i++ ) pCount += clientData[i].client.processCount();
	for( int i=0 ; i<_clientNum ; i++ ) clientData[i].ioServer = sharedIOServer;
	for( int i=0 ; i<_clientNum ; i++ ) clientData[i].showProgress = false;
	clientData[0].showProgress = true;
	SharedMemoryGrid *lowPixelReader = NULL;
	SharedMemoryGrid *pixelReader = NULL;
	SharedMemoryGrid *labelReader = NULL;
	SharedMemoryGrid *outGuessReader = NULL;
	SharedMemoryGrid *outWriter   = NULL;
	Pointer( PixelType ) lowPixelRow = NullPointer< PixelType >( );
	Pointer( PixelType ) pixelRow = NullPointer< PixelType >( );
	Pointer( LabelType ) labelRow = NullPointer< LabelType >( );
	// These are both "SyncType" as they interface later into the solver
	Pointer( SyncType ) outGuessRow = NullPointer< SyncType >();
	Pointer( SyncType )    outRow = NullPointer<  SyncType >( );
	Pointer( Signal ) lowPixelSignals = NullPointer< Signal >();
	Pointer( Signal ) outGuessSignals = NullPointer< Signal >();
	Pointer( Signal )    pixelSignals = NullPointer< Signal >();
	Pointer( Signal )    labelSignals = NullPointer< Signal >();
	Pointer( Signal )      outSignals = NullPointer< Signal >();
	int offset;
	
	pixelSignals = NewPointer< Signal >( _clientNum );
	pixelRow = AllocPointer< PixelType >( PixelChannels * size() );
	offset = 0;
	for( int i=0 ; i<_clientNum ; i++ )
	{
		clientData[i].pixels = new SharedMemoryGrid( pixelSignals+i , 1 , false , ( Pointer( byte ) )(pixelRow+offset) , clientData[i].width * PixelChannels * sizeof( PixelType ) , height() );
		offset += clientData[i].width * PixelChannels;
	}

	pixelReader = new SharedMemoryGrid( pixelSignals , _clientNum , true , ( Pointer( byte ) )pixelRow , size() * PixelChannels * sizeof( PixelType ) , height() );
	pixelReader->reset( false , 1 );
//pixels->SetServer( sharedIOServer );
	if( lowPixels )
	{
		lowPixelSignals = NewPointer< Signal >( _clientNum );
		lowPixelRow = AllocPointer< PixelType >( PixelChannels * size() );
		offset = 0;
		for( int i=0 ; i<_clientNum ; i++ )
		{
			clientData[i].lowPixels = new SharedMemoryGrid( lowPixelSignals+i , 1 , false , ( Pointer( byte ) )(lowPixelRow+offset) , clientData[i].width*PixelChannels*sizeof(PixelType) , height() );
			offset += clientData[i].width*PixelChannels;
		}
		lowPixelReader = new SharedMemoryGrid( lowPixelSignals , _clientNum , true , ( Pointer( byte ) )lowPixelRow , size()*PixelChannels*sizeof(PixelType) , height() );
		lowPixelReader->reset( false , 1 );
//lowPixel->SetServer( sharedIOServer );
	}
	if( outGuess )
	{
		outGuessSignals = NewPointer< Signal >( _clientNum );
		outGuessRow = AllocPointer< SyncType >( PixelChannels * size() );
		offset = 0;
		for( int i=0 ; i<_clientNum ; i++ )
		{
			clientData[i].outGuess = new SharedMemoryGrid( outGuessSignals+i , 1 , false , ( Pointer( byte ) )(outGuessRow+offset) , clientData[i].width*PixelChannels*sizeof(SyncType) , height() );
			offset += clientData[i].width*PixelChannels;
		}
		outGuessReader = new SharedMemoryGrid( outGuessSignals , _clientNum , true , ( Pointer( byte ) )outGuessRow , cSize()*PixelChannels*sizeof(SyncType) , cHeight() );
		outGuessReader->reset( false , 1 );
//outGuess->SetServer( sharedIOServer );
	}
	if( labels )
	{
		labelSignals = NewPointer< Signal >( _clientNum );
		labelRow = AllocPointer< LabelType >( LabelChannels * size() );
		offset = 0;
		for( int i=0 ; i<_clientNum ; i++ )
		{
			clientData[i].labels = new SharedMemoryGrid( labelSignals+i , 1 , false , ( Pointer( byte ) )(labelRow+offset) , clientData[i].width*LabelChannels*sizeof(LabelType) , height() );
			offset += clientData[i].width*LabelChannels;
		}
		labelReader = new SharedMemoryGrid( labelSignals , _clientNum , true , ( Pointer( byte ) )labelRow , size()*LabelChannels*sizeof(LabelType) , height() );
		labelReader->reset( false , 1 );
//labels->SetServer( sharedIOServer );
	}
	outSignals = NewPointer< Signal >( _clientNum );
	outRow = AllocPointer< SyncType >( PixelChannels * size() );
	offset = 0;
	for( int i=0 ; i<_clientNum ; i++ )
	{
		clientData[i].out = new SharedMemoryGrid( outSignals+i , 1 , true , ( Pointer( byte ) )(outRow+offset) , clientData[i].width*PixelChannels*sizeof( SyncType ) , height() );
		offset += clientData[i].width*PixelChannels;
	}
	outWriter = new SharedMemoryGrid( outSignals , _clientNum , false , ( Pointer( byte ) )outRow , size()*PixelChannels*sizeof( SyncType ) , height() );
	outWriter->reset( true , 1 );

	for( int i=0 ; i<_clientNum ; i++ ) clientHandles[i] = clientData[i].spawnProcess();
	for( int j=0 ; j<height() ; j++ )
	{
		int jj = std::min< int >( j , cHeight()-1 );

		{
			Pointer( PixelType ) pOut = ( Pointer( PixelType ) )( (*pixelReader)[j] );
			Pointer( PixelType ) pIn = ( Pointer( PixelType ) )(*pixels)[jj];
			memcpy( pOut , pIn , sizeof( PixelType ) * PixelChannels * cSize() );
			for( int i=cSize() ; i<size() ; i++ ) for( int c=0 ; c<PixelChannels ; c++ ) pOut[i*PixelChannels+c] = pIn[ (cSize()-1)*PixelChannels+c ];
			if( jj!=cHeight()-1 ) pixels->advance();
			pixelReader->advance();
		}
		if( lowPixels )
		{
			Pointer( PixelType ) lpOut = ( Pointer( PixelType ) )(*lowPixelReader)[j];
			Pointer( PixelType ) lpIn = ( Pointer( PixelType ) )(*lowPixels)[jj];
			memcpy( lpOut , lpIn , sizeof( PixelType ) * PixelChannels * cSize() );
			for( int i=cSize() ; i<size() ; i++ ) for( int c=0 ; c<PixelChannels ; c++ ) lpOut[i*PixelChannels+c] = lpIn[ (cSize()-1)*PixelChannels+c ];
			if( jj!=cHeight()-1 ) lowPixels->advance();
			lowPixelReader->advance();
		}
		if( labels )
		{
			Pointer( LabelType ) lOut = ( Pointer( LabelType ) )(*labelReader)[j];
			Pointer( LabelType ) lIn = ( Pointer( LabelType ) )(*labels)[jj];
			memcpy( lOut , lIn , sizeof( LabelType ) * LabelChannels * cSize() );
			for( int i=cSize() ; i<size() ; i++ ) for( int c=0 ; c<LabelChannels ; c++ ) lOut[i*LabelChannels+c] = lIn[ (cSize()-1)*LabelChannels+c ];
			if( jj!=cHeight()-1 ) labels->advance();
			labelReader->advance();
		}
		if( outGuess )
		{
			Pointer( SyncType ) gOut = ( Pointer( SyncType ) )(*outGuessReader)[j];
			Pointer( PixelType ) gIn = ( Pointer( PixelType ) )(*outGuess)[jj];
			ConvertRow< PixelType , SyncType >( gIn , gOut , cSize() , PixelChannels , PixelChannels );
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				SyncType v = ConvertChannel< PixelType , SyncType >( gIn[ (cSize()-1)*PixelChannels+c ] );
				for( int i=cSize() ; i<size() ; i++ ) gOut[i*PixelChannels+c] = v;
			}
			if( jj!=cHeight()-1 ) outGuess->advance();
			outGuessReader->advance();
		}
	}
	pixels->advance();
	if( lowPixels ) lowPixels->advance();
	if( labels ) labels->advance();
	if( outGuess ) outGuess->advance();
//out->SetServer( sharedIOServer );
	for( int j=0 ; j<height() ; j++ )
	{
		Pointer( SyncType ) oIn = ( Pointer( SyncType ) )(*outWriter)[j];
		if( j<cHeight() )
		{
			Pointer( PixelType ) oOut = ( Pointer( PixelType ) )(*out)[j];
			ConvertRow< SyncType , PixelType >( oIn , oOut , cSize() , PixelChannels , PixelChannels );
			if( j<cHeight()-1 ) out->advance();
		}
		outWriter->advance();
	}
	out->advance();
	WaitOnThreads( clientHandles , _clientNum , 10000 , "SocketedSuperMultigridClient::_Run (Process)" );
	DeletePointer( clientHandles );
	DeletePointer( clientData );
	// Delete client data before deleting all the arrays it depends on
	FreePointer( lowPixelRow );
	FreePointer( outGuessRow );
	FreePointer( pixelRow );
	FreePointer( labelRow );
	FreePointer(   outRow );
	DeletePointer( lowPixelSignals );
	DeletePointer( outGuessSignals );
	DeletePointer(    pixelSignals );
	DeletePointer(    labelSignals );
	DeletePointer(      outSignals );
	if( lowPixelReader ) delete lowPixelReader;
	if( outGuessReader ) delete outGuessReader;
	if( pixelReader ) delete pixelReader;
	if( labelReader ) delete labelReader;
	if(   outWriter ) delete   outWriter;
}

template< class StorageType , class PixelType , class LabelType >
bool SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::gammaCorrection( void ) const { return _globalData.gammaCorrection; }
template< class StorageType , class PixelType , class LabelType >
bool SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::hdr( void ) const { return _globalData.hdr; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::tileWidth( void ) const { return _globalData.tileWidth; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::tileHeight( void ) const { return _globalData.tileHeight; }
template< class StorageType , class PixelType , class LabelType >
const char* SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::tileExtension( void ) const { return _globalData.tileExt; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::quality( void ) const { return _globalData.quality; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::cSize( void ) const { return _clientData.cEnd - _clientData.start; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::size( void ) const { return _clientData.end - _clientData.start; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::cHeight( void ) const { return _globalData.cHeight; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::height( void ) const { return _globalData.height; }
template< class StorageType , class PixelType , class LabelType >
bool SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::showProgress( void ) const { return _globalData.showProgress; }
template< class StorageType , class PixelType , class LabelType >
bool SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::shortSync( void ) const { return _globalData.shortSync; }
template< class StorageType , class PixelType , class LabelType >
int SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::periodicType( void ) const { return _globalData.periodicType; }
template< class StorageType , class PixelType , class LabelType >
bool SocketedSuperMultigridClient< StorageType , PixelType , LabelType >::verbose( void ) const { return _globalData.verbose; }
