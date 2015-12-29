/*
Copyright (c) 2010, Michael Kazhdan
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
#ifndef XPLATFORM_INCLUDED
#define XPLATFORM_INCLUDED

#define DEBUG_MUTEX 1

#ifdef _WIN32
#define ALIGN( declaration , bytes ) __declspec(align(bytes)) declaration
const char FileSeparator = '\\';
#else // !_WIN32
#define ALIGN( declaration , bytes ) declaration __attribute__((aligned(bytes)))
const char FileSeparator = '/';
#endif // _WIN32

#ifdef _WIN32
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0501
#endif // !_WIN32_WINNT
#endif // _WIN32

#include <Util/Array.h>

#include <stdint.h>
typedef uint8_t byte;

// locking
#include <boost/thread/recursive_mutex.hpp>
#if DEBUG_MUTEX
typedef boost::recursive_timed_mutex CriticalSectionLock;
#else // !DEBUG_MUTEX
typedef boost::recursive_mutex CriticalSectionLock;
#endif
inline void InitializeCriticalSection( CriticalSectionLock* ){ ; }
inline void     DeleteCriticalSection( CriticalSectionLock* ){ ; }
#if DEBUG_MUTEX
inline void      EnterCriticalSection( CriticalSectionLock* l , long long ms , const char* message=NULL )
{
	if( ms<=0 ) l->lock();
	else
		while( 1 )
			if( l->timed_lock( boost::posix_time::milliseconds( ms ) ) ) return;
			else
				if( message ) fprintf( stderr , "[WARNING] Failed to acquire lock: %s\n" , message );
				else          fprintf( stderr , "[WARNING] Failed to acquire lock\n" );
}
#else // !DEBUG_MUTEX
inline void      EnterCriticalSection( CriticalSectionLock* l ){ l->lock(); }
#endif // DEBUG_MUTEX
inline void      LeaveCriticalSection( CriticalSectionLock* l ){ l->unlock(); }

#include <boost/asio.hpp>
class IOServer
{
	static int _count;
	static CriticalSectionLock _stdoutLock , _stderrLock , _systemLock;
public:
	static boost::asio::io_service io_service;
	static void Load( void );
	static void UnLoad( void );
	class StdoutLock
	{
	public:
		StdoutLock ( void );
		~StdoutLock( void );
	};
	class StderrLock
	{
	public:
		StderrLock ( void );
		~StderrLock( void );
	};
	class SystemLock
	{
	public:
		SystemLock ( void );
		~SystemLock( void );
	};
	static void  printfID(            const char* format , ... );
	static void fprintfID( FILE* fp , const char* format , ... );
};

// sockets
#include <boost/asio.hpp>
typedef boost::asio::ip::tcp::socket* Socket;
typedef boost::asio::ip::tcp::acceptor* AcceptorSocket;
typedef boost::asio::ip::address EndpointAddress;
const Socket _INVALID_SOCKET_ = (Socket)NULL;
const AcceptorSocket _INVALID_ACCEPTOR_SOCKET_ = (AcceptorSocket)NULL;
static boost::asio::io_service io_service;
template< class C > int socket_receive( Socket& s , C* destination , size_t len )
{
	boost::system::error_code ec;
	int ret = (int)( boost::asio::read( *s , boost::asio::buffer( destination , len ) , ec ) );
	if( ec ) return  -1;
	else     return ret;
}
template< class C > int socket_send( Socket& s , const C* source , size_t len )
{
	boost::system::error_code ec;
	int ret = (int)( boost::asio::write( *s , boost::asio::buffer( source , len ) , ec ) );
	if( ec ) return  -1;
	else     return ret;
}
inline bool AddressesEqual( const EndpointAddress& a1 ,  const EndpointAddress& a2 ){ return a1.to_string()==a2.to_string(); }
inline const char *LastSocketError( void ){ return ""; }
void PrintHostAddresses( FILE* fp=stdout );

// signals
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
struct Signal
{
	boost::mutex mutex;
	bool state;
	boost::condition_variable condition;
	Signal( void ){ state=true; }
	Signal( const Signal& s )
	{
		state = s.state;
		fprintf( stderr , "[WARNING] Could not copy signal mutex/condition" );
	}
};
// Used to ping-pong between two clients.
// The client whose state is asigned "true" goes first.
struct Signaller
{
	Signal& signal;
	bool state;
	static void Start( Signal& signal , bool state , long long ms=100 , bool verbose=false )
	{
		boost::unique_lock< boost::mutex > lock( signal.mutex );
		int count = 0;
		while( signal.state!=state )
		{
			signal.condition.timed_wait( lock , boost::posix_time::milliseconds(ms) );
			if( count && verbose ) IOServer::printfID( "waiting on signal (%d)\n" , count );
			count++;
		}
	}
	static void Stop( Signal& signal , bool state )
	{
		{
			boost::unique_lock< boost::mutex > lock( signal.mutex );
			if( signal.state!=state ) fprintf( stderr , "[ERROR] Trying to release other signal\n" ) , exit( 0 );
			signal.state = !state;
		}
		signal.condition.notify_one();
	}
	Signaller( Signal& si , bool st ) : signal( si ) , state( st ) { Start( signal , state ); }
	~Signaller( void ){ Stop( signal , state ); }
};

// threads
#include <boost/thread.hpp>
typedef boost::thread* ThreadHandle;
inline void SetThisThreadID( char* id )
{
	std::stringstream stream;
	stream << boost::this_thread::get_id();
	sprintf( id , "%s" , stream.str().c_str() );
}
inline void CloseThreadHandle( ThreadHandle& t ){ t->join() ; delete t ; t=NULL; }
inline bool TestThreadHandle( ThreadHandle& t ) { return t!=NULL; }
inline ThreadHandle RunThread( int ThreadFunction( void* ) , void* params ){ return new boost::thread( ThreadFunction , params ); }
inline bool InterruptThread( ThreadHandle& t ){ t->interrupt() ; CloseThreadHandle( t ) ; return true; }
inline bool WaitOnThread( ThreadHandle& t , int ms=-1 , const char* message=NULL )
{
	if( !t->joinable() ) return true;
	if( ms<0 ) { CloseThreadHandle( t ) ; return true; }
	else
		if( !t->timed_join( boost::posix_time::milliseconds(ms) ) )
		{
			if( message ) fprintf( stderr , "[WARNING] Failing to wait on thread: %s\n" , message );
			else          fprintf( stderr , "[WARNING] Failing to wait on thread\n" );
			return false;
		}
		else return true;
}
inline bool WaitOnThreads( ThreadHandle* t , int tCount , int ms=-1 , const char* message=NULL )
{
	for( int i=0 ; i<tCount ; i++ )
	{
		if( t[i]->joinable() )
		{
			if( ms<0 ) CloseThreadHandle( t[i] );
			while( t[i]->joinable() && !t[i]->timed_join( boost::posix_time::milliseconds(ms) ) )
			{
				if( message ) fprintf( stderr , "[WARNING] Failing to wait on thread[%d/%d]: %s\n" , i , tCount , message );
				else          fprintf( stderr , "[WARNING] Failing to wait on thread[%d/%d]\n" , i  , tCount );
			}
		}
	}
	return true;
}
inline bool IsActiveThread( ThreadHandle& t ) { return t->joinable(); }
#ifdef OLD_BOOST
inline void SleepThisThread( int ms ){ if(ms<0) boost::this_thread::yield() ; else boost::this_thread::sleep( boost::posix_time::milliseconds( ms ) ); }
#else // !OLD_BOOST
inline void SleepThisThread( int ms ){ if(ms<0) boost::this_thread::yield() ; else boost::this_thread::sleep_for( boost::chrono::milliseconds(ms) ); }
#endif // OLD_BOOST

#ifdef _WIN32
inline int GetThisProcessID( void ){ return GetCurrentProcessId(); }
#ifndef strcasecmp
#define strcasecmp _stricmp
#endif // strcasecmp
#else // !_WIN32
inline int GetThisProcessID( void ){ return getpid(); }
#endif // _WIN32


#ifdef _WIN32
#include <io.h>
#endif // _WIN32
struct _FileHandle
{
	_FileHandle( void ) : _fp( NULL ) , _deleteOnClose( false ){ _fileName[0]=0; }
	~_FileHandle( void )
	{
		if( _fp )
		{
			fclose( _fp );
			if( _fp && _deleteOnClose && strlen( _fileName ) ) remove( _fileName );
		}
		_fp = NULL;
		_fileName[0] = 0;
	}
	const char* fileName( void ) const { return _fileName; }
	bool deleteOnClose( void ) const { return _deleteOnClose; }
	_FileHandle( const char* fileName )
	{
		strcpy( _fileName , fileName );
#ifdef _WIN32
		fopen_s( &_fp , _fileName , "w+b" );
#else // !_WIN32
		_fp = fopen( fileName , "w+b" );
#endif // _WIN32
		if( !_fp ) fprintf( stderr , "[ERROR] Failed to open file: %s\n" , _fileName ) , exit( 0 );
		_deleteOnClose = false;
		setbuf( _fp , NULL );
	}
	_FileHandle( const char* dir , const char* prefix , bool deleteOnClose )
	{
		const char scratch[] = "scratch";
		const char current[] = ".";
		if( !prefix ) prefix = scratch;
		{
			IOServer::SystemLock lock;
			if( !dir || !strlen(dir) ) dir = getenv( "TMP" );	// Try to get the scratch directory from the environment variable
			if( !dir || !strlen(dir) ) dir = current;			// Set it to the current directory
			if( dir[strlen(dir)-1]!=FileSeparator ) sprintf( _fileName , "%s%c%s_%d_XXXXXX" , dir , FileSeparator , prefix , GetThisProcessID() );
			else                                    sprintf( _fileName ,   "%s%s_%d_XXXXXX" , dir ,                 prefix , GetThisProcessID() );
		}
#ifdef _WIN32
		_mktemp( _fileName );
		fopen_s( &_fp , _fileName , "w+b" );
#else // !_WIN32
		mktemp( _fileName );
		_fp = fopen( _fileName , "w+b" );
#endif // _WIN32
		if( !_fp ) fprintf( stderr , "[ERROR] Failed to open file: %s\n" , _fileName ) , exit( 0 );
		_deleteOnClose = deleteOnClose;
		setbuf( _fp , NULL );
	}
	bool seek( long long offset )
	{
#ifdef _WIN32
		return _fseeki64( _fp , offset , SEEK_SET )==0;
#else // !_WIN32
		return  fseek   ( _fp , offset , SEEK_SET )==0;
#endif // _WIN32
	}
	long long read ( Pointer( byte ) buffer , long long bytesToRead  ){ return fread ( buffer , sizeof(byte) , bytesToRead  , _fp ); }
	long long write( Pointer( byte ) buffer , long long bytesToWrite ){ return fwrite( buffer , sizeof(byte) , bytesToWrite , _fp ); }
private:
	FILE* _fp;
	bool _deleteOnClose;
	char _fileName[1024];
};
typedef _FileHandle* FileHandle;
inline FileHandle CreateFileHandle( const char* dir , const char* prefix , bool deleteOnClose ){ return new _FileHandle( dir , prefix , deleteOnClose ); }
inline FileHandle CreateFileHandle( const char* fileName ){ return new _FileHandle( fileName ); }
inline void CloseFileHandle( FileHandle& hFile ){ if( hFile ) delete hFile ; hFile = NULL; }
inline bool SeekFileHandle( FileHandle hFile , long long  distanceToMove ){ return hFile->seek( distanceToMove ); }
inline bool SetEndOfFileHandle( FileHandle hFile ){ return true; }
inline long long ReadFileHandle( FileHandle hFile , Pointer( byte ) buffer , long long bytesToRead ){ return hFile->read( buffer , bytesToRead ); }
inline long long WriteFileHandle( FileHandle hFile , Pointer( byte ) buffer , long long bytesToWrite ){ return hFile->write( buffer , bytesToWrite ); }




#if ARRAY_DEBUG
template< class C >
int socket_receive( Socket& s , Array< C > destination , size_t len )
{
	if( len>destination.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of socket_receive exceeds destination maximum: %zd > %zd\n" , len , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return socket_receive( s , (char*)&destination[0] , len );
}
template< class C >
int socket_send( Socket s , ConstArray< C > source , size_t len )
{
	if( len>source.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of socket_send exceeds source maximum: %zd > %zd\n" , len , source.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return socket_send( s , (char*)&source[0] , len );
}
#endif // ARRAY_DEBUG
#endif // XPLATFORM_INCLUDED
