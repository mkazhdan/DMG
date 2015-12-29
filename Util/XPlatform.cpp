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
#include <stdio.h>
#include <stdarg.h>
#include "XPlatform.h"

int IOServer::_count = 0;
CriticalSectionLock IOServer::_stdoutLock , IOServer::_stderrLock , IOServer::_systemLock;

#if DEBUG_MUTEX
IOServer::StdoutLock::StdoutLock ( void ){ EnterCriticalSection( &_stdoutLock , 100 , "StdoutLock::StdoutLock" ); }
#else // !DEBUG_MUTEX
IOServer::StdoutLock::StdoutLock ( void ){ EnterCriticalSection( &_stdoutLock ); }
#endif // DEBUG_MUTEX
IOServer::StdoutLock::~StdoutLock( void ){ LeaveCriticalSection( &_stdoutLock ); }
#if DEBUG_MUTEX
IOServer::StderrLock::StderrLock ( void ){ EnterCriticalSection( &_stderrLock , 100 , "StderrLock::StderrLock" ); }
#else // !DEBUG_MUTEX
IOServer::StderrLock::StderrLock ( void ){ EnterCriticalSection( &_stderrLock ); }
#endif // DEBUG_MUTEX
IOServer::StderrLock::~StderrLock( void ){ LeaveCriticalSection( &_stderrLock ); }
#if DEBUG_MUTEX
IOServer::SystemLock::SystemLock ( void ){ EnterCriticalSection( &_systemLock , 100 , "SystemLock::SystemLock" ); }
#else // !DEBUG_MUTEX
IOServer::SystemLock::SystemLock ( void ){ EnterCriticalSection( &_systemLock ); }
#endif // DEBUG_MUTEX
IOServer::SystemLock::~SystemLock( void ){ LeaveCriticalSection( &_systemLock ); }

boost::asio::io_service IOServer::io_service;
void IOServer::Load( void )
{
	if( !_count )
	{
		InitializeCriticalSection( &_stdoutLock );
		InitializeCriticalSection( &_stderrLock );
		InitializeCriticalSection( &_systemLock );
	}
	_count++;
}
void IOServer::UnLoad( void )
{
	_count--;
	if( !_count )
	{
		DeleteCriticalSection( &_stdoutLock );
		DeleteCriticalSection( &_stderrLock );
		DeleteCriticalSection( &_systemLock );
	}
}

void IOServer::fprintfID( FILE* fp , const char* format , ... )
{
	StdoutLock lock;
	va_list args;
	va_start( args , format );
	char id[512] ; SetThisThreadID( id );
	fprintf( fp , "%s] " , id );
	vfprintf( fp , format , args );
	va_end( args );
	fflush( fp );
}

void IOServer::printfID( const char* format , ... )
{
	StdoutLock lock;
	va_list args;
	va_start( args , format );
	char id[512] ; SetThisThreadID( id );
	printf( "%s] " , id );
	vprintf( format , args );
	va_end( args );
	fflush( stdout );
}


void PrintHostAddresses( FILE* fp )
{
	boost::asio::ip::tcp::resolver resolver( io_service );
	boost::asio::ip::tcp::resolver::query query( boost::asio::ip::host_name() , std::string( "" ) , boost::asio::ip::resolver_query_base::numeric_service );
	boost::asio::ip::tcp::resolver::iterator iterator = resolver.resolve( query ) , end;
	for( int count=0 ; iterator!=end ; )
	{
		if( (*iterator).endpoint().address().is_v4() ) fprintf( fp , "%d]  %s\n" , count++ , (*iterator).endpoint().address().to_string().c_str() );
//		else                                           fprintf( fp , "%d]* %s\n" , count++ , (*iterator).endpoint().address().to_string().c_str() );
		iterator++;
	}
}
