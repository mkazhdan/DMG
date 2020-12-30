/*
Copyright (c) 2008, Michael Kazhdan
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

#include "Socket.h"

void printfId( const char* format , ... )
{
	va_list args;
	va_start(args,format);
	char id[512] ; SetThisThreadID( id );
	printf( "%s] " , id );
	vprintf(format,args);
	va_end(args);
}
void fprintfId(FILE* fp , const char* format,...)
{
	va_list args;
	va_start( args , format );
	char id[512] ; SetThisThreadID( id );
	fprintf( fp , "%s] " , id );
	vfprintf( fp , format , args );
	va_end( args);
}

bool GetHostEndpointAddress( EndpointAddress* address , const char* prefix )
{
	boost::asio::ip::tcp::resolver resolver( io_service );
	boost::asio::ip::tcp::resolver::query query( boost::asio::ip::host_name() , std::string( "" ) , boost::asio::ip::resolver_query_base::numeric_service );
	boost::asio::ip::tcp::resolver::iterator iterator = resolver.resolve( query ) , end;
	for( int count=0 ; iterator!=end ; )
	{
		if( (*iterator).endpoint().address().is_v4() )
		{
			std::string str = (*iterator).endpoint().address().to_string();
			const char* _address = str.c_str();
			if( !prefix || strstr( _address , prefix ) )
			{
				*address = (*iterator).endpoint().address();
				return true;
			}
		}
		iterator++;
	}
	return false;
}
bool GetHostAddress( char* address , const char* prefix )
{
	EndpointAddress _address;
	if( !GetHostEndpointAddress( &_address , prefix ) ) return false;
		strcpy( address , _address.to_string().c_str() );
	return true;
}
int GetLocalSocketPort( Socket& s )
{
	return s->local_endpoint().port();
}
EndpointAddress GetLocalSocketEndpointAddress( Socket& s )
{
	return s->local_endpoint().address();
}
int GetPeerSocketPort( Socket& s )
{
	return s->remote_endpoint().port();
}
EndpointAddress GetPeerSocketEndpointAddress( Socket& s )
{
	return s->remote_endpoint().address();
}
Socket GetConnectSocket( const char* address , int port , int ms , bool progress )
{
	char _port[128];
	sprintf( _port , "%d" , port );
	boost::asio::ip::tcp::resolver resolver( io_service );
	boost::asio::ip::tcp::resolver::query query( address , _port );
	boost::asio::ip::tcp::resolver::iterator iterator = resolver.resolve( query );
	Socket s = new boost::asio::ip::tcp::socket( io_service );
	boost::system::error_code ec;
	long long sleepCount = 0;
	do
	{
		boost::asio::connect( *s , resolver.resolve(query) , ec );
		sleepCount++;
		SleepThisThread( 1 );
		if( progress && !(sleepCount%ms) ) printf( "." );
	}
	while( ec );
	if( progress ) printf( "\n" ) , fflush( stdout );
	return s;
}
Socket GetConnectSocket( EndpointAddress address , int port , int ms , bool progress )
{
	char _port[128];
	sprintf( _port , "%d" , port );
	boost::asio::ip::tcp::resolver resolver( io_service );
	boost::asio::ip::tcp::resolver::query query( address.to_string().c_str() , _port );
	boost::asio::ip::tcp::resolver::iterator iterator = resolver.resolve( query );
	Socket s = new boost::asio::ip::tcp::socket( io_service );
	boost::system::error_code ec;
	long long sleepCount = 0;
	do
	{
		boost::asio::connect( *s , resolver.resolve(query) , ec );
		sleepCount++;
		SleepThisThread( 1 );
		if( progress && !(sleepCount%ms) ) printf( "." );
	}
	while( ec );
	if( progress ) printf( "\n" ) , fflush( stdout );
	return s;
}
Socket AcceptSocket( AcceptorSocket listen )
{
	Socket s = new boost::asio::ip::tcp::socket( io_service );
	listen->accept( *s );
	return s;
}

AcceptorSocket GetListenSocket( int& port )
{
	AcceptorSocket s = new boost::asio::ip::tcp::acceptor( io_service , boost::asio::ip::tcp::endpoint( boost::asio::ip::tcp::v4() , port ) );
	port = s->local_endpoint().port();
	return s;
}
void CloseSocket( Socket& s )
{
	delete s;
	s = _INVALID_SOCKET_;
}
void CloseAcceptorSocket( AcceptorSocket& s )
{
	delete s;
	s = _INVALID_ACCEPTOR_SOCKET_;
}

///////////////////////////
// DataStreamConstructor //
///////////////////////////
DataStreamConstructor::DataStreamConstructor( void )
{
	_sock = _INVALID_SOCKET_;
	_stream = NULL;
	GetHostEndpointAddress( &_myAddr );
}
void DataStreamConstructor::init( Socket sock , bool master , bool cleanUp )
{
	_sock = sock;
	_master = master;
	_cleanUp = cleanUp;
}
DataStream* DataStreamConstructor::getDataStream( void ) { return _stream; }
void DataStreamConstructor::doStep( int sNum )
{
	switch( sNum )
	{
	case 0:
		if( !_master )
		{
			SendOnSocket( _sock , ( ConstPointer( EndpointAddress ) )GetPointer(_myAddr) , sizeof(_myAddr) );
			SendOnSocket( _sock , ( ConstPointer( int ) )GetPointer(_myPID ) , sizeof(_myPID) );
		}
		break;
	case 1:
		if( _master )
		{
			EndpointAddress addr;
			int pid;
			ReceiveOnSocket( _sock , GetPointer(addr) , sizeof(addr) );
			ReceiveOnSocket( _sock , GetPointer(pid ) , sizeof(pid) );

			if( _myAddr.to_string()==addr.to_string() && _myPID == pid )
			{
				SharedMemoryBuffer::StreamPair sPair;
				SharedMemoryBuffer::StreamPair::CreateSharedBufferPair( sPair );
				SendOnSocket( _sock , ( ConstPointer( SharedMemoryBuffer::SecondStream* ) )GetPointer(sPair.second) , sizeof( sPair.second ) );
				_stream = sPair.first;
				if( _cleanUp ) CloseSocket( _sock );
			}
			else
			{
				SharedMemoryBuffer::SecondStream* sStream = NULL;
				SendOnSocket( _sock , ( ConstPointer( SharedMemoryBuffer::SecondStream* ) )GetPointer(sStream) , sizeof( sStream ) );
				_stream = new SocketStream( _sock );
			}
		}
		break;
	case 2:
		if( !_master )
		{
			SharedMemoryBuffer::SecondStream* sStream;
			ReceiveOnSocket( _sock , GetPointer(sStream) , sizeof( sStream ) );
			if( sStream )
			{
				if( _cleanUp ) CloseSocket( _sock );
				_stream = sStream;
			}
			else _stream = new SocketStream( _sock  );
		}
		break;
	}
}


////////////////
// DataStream //
////////////////
DataStream* DataStream::GetDataStream( Socket sock , bool master , bool cleanUp )
{
	char address[512];
	EndpointAddress myAddr , addr;
	int pid , myPID = GetThisProcessID( );
	GetHostAddress( address );
	myAddr.from_string( address );

	if( master )
	{
		ReceiveOnSocket( sock , GetPointer(addr) , sizeof(addr) );
		ReceiveOnSocket( sock , GetPointer(pid ) , sizeof(pid) );
		if( myAddr==addr && myPID==pid )
		{
			SharedMemoryBuffer::StreamPair sPair;
			SharedMemoryBuffer::StreamPair::CreateSharedBufferPair( sPair );
			SendOnSocket( sock , ( ConstPointer( SharedMemoryBuffer::SecondStream* ) )GetPointer(sPair.second) , sizeof( sPair.second ) );
			if( cleanUp ) CloseSocket( sock );
			return sPair.first;
		}
		else
		{
			SharedMemoryBuffer::SecondStream* sStream = NULL;
			SendOnSocket( sock , ( ConstPointer( SharedMemoryBuffer::SecondStream* ) )GetPointer(sStream) , sizeof( sStream ) );
			return new SocketStream( sock );
		}
	}
	else
	{
		SendOnSocket( sock , ( ConstPointer( EndpointAddress ) )GetPointer(myAddr) , sizeof(myAddr) );
		SendOnSocket( sock , ( ConstPointer( int ) )GetPointer(myPID) , sizeof(myPID) );
		SharedMemoryBuffer::SecondStream* sStream;
		ReceiveOnSocket( sock , GetPointer(sStream) , sizeof( sStream ) );
		if( sStream )
		{
			if( cleanUp ) CloseSocket( sock );
			return sStream;
		}
		else return new SocketStream( sock  );
	}
}
////////////
// Socket //
////////////
SocketStream::SocketStream( Socket sock ) { _sock = sock; _mySocket = false; }
SocketStream::SocketStream( const char* address , int port , int ms , bool progress )
{
	_sock = GetConnectSocket( address , port , ms , progress );
	_mySocket = true;
}
SocketStream::SocketStream( EndpointAddress address , int port , int ms , bool progress )
{
	_sock = GetConnectSocket( address , port , ms , progress );
	_mySocket = true;
}
SocketStream::~SocketStream( void ) { if( _mySocket ) CloseSocket( _sock); }

bool SocketStream::write( ConstPointer( byte ) buf , int len ) { return SendOnSocket   ( _sock , buf , len ); }
bool SocketStream::read ( Pointer(       byte ) buf , int len ){ return ReceiveOnSocket( _sock , buf , len ); }

////////////////////////
// SharedMemoryBuffer //
////////////////////////
SharedMemoryBuffer::SharedMemoryBuffer( void )
{
	_buf1 = NullPointer< byte >( );
	_buf2 = NullPointer< byte >( );
	_bufSize1 = _bufSize2 = 0;
}
SharedMemoryBuffer::~SharedMemoryBuffer( void )
{
	FreePointer( _buf1 ) ; _bufSize1 = 0;
	FreePointer( _buf2 ) ; _bufSize2 = 0;
}
SharedMemoryBuffer::FirstStream::FirstStream  ( SharedMemoryBuffer* smb ) { _smb = smb; }
SharedMemoryBuffer::SecondStream::SecondStream( SharedMemoryBuffer* smb ) { _smb = smb; }
SharedMemoryBuffer::FirstStream::~FirstStream  ( void ) { if( _smb ) delete _smb , _smb = NULL; }
SharedMemoryBuffer::SecondStream::~SecondStream( void ) {                          _smb = NULL; }
// The first stream reads on buf1 and writes on buf2
bool SharedMemoryBuffer::FirstStream::read( Pointer( byte ) buf , int len ) 
{
	Signaller signaller( _smb->_readyForWriting1 , false );
	if( len>_smb->_bufSize1 )
	{
		printf( "Uh oh 1\n" ) , fflush( stdout );
		return false;
	}
	memcpy( buf , _smb->_buf1 , len );
	return true;
}
bool SharedMemoryBuffer::FirstStream::write( ConstPointer( byte ) buf , int len )
{
	Signaller signaller( _smb->_readyForWriting2 , true );
	if( len>_smb->_bufSize2 )
	{
		FreePointer( _smb->_buf2 );
		_smb->_bufSize2 = 0;
		_smb->_buf2 = AllocPointer< byte >( len );
		if( !_smb->_buf2 )
		{
			printf( "Uh oh 2\n" ) , fflush( stdout );
			return false;
		}
		_smb->_bufSize2 = len;
	}
	memcpy( _smb->_buf2 , buf , len );
	return true;
}
bool SharedMemoryBuffer::SecondStream::read( Pointer( byte ) buf , int len )
{
	Signaller signaller( _smb->_readyForWriting2 , false );
	if( len>_smb->_bufSize2 )
	{
		printf( "Uh oh 3\n" ) , fflush( stdout );
		return false;
	}
	memcpy( buf , _smb->_buf2 , len );
	return true;
}
bool SharedMemoryBuffer::SecondStream::write( ConstPointer( byte ) buf , int len )
{
	Signaller signaller( _smb->_readyForWriting1 , true );
	if( len>_smb->_bufSize1 )
	{
		FreePointer( _smb->_buf1 );
		_smb->_bufSize1 = 0;
		_smb->_buf1 = AllocPointer< byte >( len );
		if( !_smb->_buf1 )
		{
			printf( "Uh oh 4\n" ) , fflush( stdout );
			return false;
		}
		_smb->_bufSize1 = len;
	}
	memcpy( _smb->_buf1 , buf , len );
	return true;
}
bool SharedMemoryBuffer::StreamPair::CreateSharedBufferPair( StreamPair& pair )
{
	SharedMemoryBuffer* smb = new SharedMemoryBuffer( );
	pair.first  = new SharedMemoryBuffer::FirstStream ( smb );
	pair.second = new SharedMemoryBuffer::SecondStream( smb );
	if( !pair.first || !pair.second )
	{
		fprintf( stderr , "Failed to create shared buffer pair\n" );
		return false;
	}
	return true;
}

SharedMemoryBuffer::StreamPair::StreamPair( void )
{
	first = NULL;
	second = NULL;
}
