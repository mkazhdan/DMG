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
#include "XPlatform.h"
#include "BaseMultiStreamIO.h"
#include "Socket.h"

//////////////
// IOClient //
//////////////
const int IOClient::BYTES_PER_SECTOR = 1<<9;					// 512B Sector sizes
//const int IOClient::IO_BLOCK_SIZE = BYTES_PER_SECTOR<<13;		// 4MB IO Chunks
const int IOClient::IO_BLOCK_SIZE = BYTES_PER_SECTOR<<12;		// 2MB IO Chunks
long long IOClient::ReadBytes  = 0;
long long IOClient::WriteBytes = 0;

IOClient::IOClient( void )
{
	InitializeCriticalSection( &lock );
	server = NULL;
}
IOClient::~IOClient( void )
{
	DeleteCriticalSection( &lock );
}
void IOClient::SetServer( class MultiStreamIOServer* server )
{
	this->server = server;
	if( server ) server->AddClient( this );
}
/////////////////////////
// MultiStreamIOServer //
/////////////////////////
MultiStreamIOServer::MultiStreamIOServer( void )
{
	InitializeCriticalSection( &_pendingLock );
	InitializeCriticalSection( &_clientLock );

	_allowNewClients = true;
	_pendingClient = NULL;
	_ioThread = RunThread( _IOThread , this );
	if( !TestThreadHandle( _ioThread ) ) fprintf( stderr , "[ERROR] MultiStreamIOServer::MultiStreamIOServer: Failed to create I/O thread\n" ) , exit( 0 );
}
MultiStreamIOServer::~MultiStreamIOServer( void )
{
	// First get _ioThread to terminate, if it's still active
	if( IsActiveThread( _ioThread ) )
	{
		while( 1 )
		{
			// As long as the server has work to do, it cannot be terminated.
#if DEBUG_MUTEX
			EnterCriticalSection( &_clientLock , 100 , "MultiStreamIOServer::~MultiStreamIOServer" );
#else // !DEBUG_MUTEX
			EnterCriticalSection( &_clientLock );
#endif // DEBUG_MUTEX
			if( !_clients.size() )
			{
				// [WARNING] Dead-lock potential. Trying to interrupt while holding onto _clientLock.
				// _IOThread may block waiting for the lock before reaching the interrupt point.
				_allowNewClients = false;
				LeaveCriticalSection( &_clientLock );
				if( !InterruptThread( _ioThread ) ) fprintf( stderr , "[WARNING] MultiStreamIOServer::~MultiStreamIOServer: Failed to interrupt MultiStreamIOServer thread\n" );
				break;
			}
			LeaveCriticalSection( &_clientLock );
		}
	}

	DeleteCriticalSection( &_pendingLock );
	DeleteCriticalSection( &_clientLock );
}
bool MultiStreamIOServer::SetPending( IOClient* client )
{
#if DEBUG_MUTEX
	EnterCriticalSection( &_pendingLock , 100 , "MultiStreamIOServer::SetPending" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &_pendingLock );
#endif // DEBUG_MUTEX
	if( !_pendingClient && client )
	{
		_pendingClient = client;
		LeaveCriticalSection( &_pendingLock );
		return true;
	}
	else
	{
		LeaveCriticalSection( &_pendingLock );
		return false;
	}
}
void MultiStreamIOServer::AddClient( IOClient* client )
{
#if DEBUG_MUTEX
	EnterCriticalSection( &_clientLock , 100 , "MultiStreamIOServer::AddClient" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &_clientLock );
#endif // DEBUG_MUTEX
	if( _allowNewClients ) _clients.push_back( client );
	else fprintf( stderr , "[ERROR] MultiStreamIOServer::AddClient: Adding new clients disallowed\n" );
	LeaveCriticalSection( &_clientLock );
}
int MultiStreamIOServer::clientNum( void )
{
#if DEBUG_MUTEX
	EnterCriticalSection( &_clientLock , 100 , "MultiStreamIOServer::clientNum" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &_clientLock );
#endif // DEBUG_MUTEX
	int sz = int( _clients.size() );
	LeaveCriticalSection( &_clientLock );
	return sz;
}
int MultiStreamIOServer::_IOThread( void* vparams )
{
static int threadID=0;
int myID = threadID++;
	MultiStreamIOServer* server = (MultiStreamIOServer*)vparams;
	std::vector< IOClient* >& clients = server->_clients;
	int idx = 0;
	{
		while( 1 )
		{
			if( boost::this_thread::interruption_requested() ) break;
			bool hasPending = false;
			{
#if DEBUG_MUTEX
				EnterCriticalSection( &server->_pendingLock , 100 , "MultiStreamIOServer::_IOThread (_pendingLock #1)" );
#else // !DEBUG_MUTEX
				EnterCriticalSection( &server->_pendingLock );
#endif // DEBUG_MUTEX
				if( server->_pendingClient )
				{
					if( server->_pendingClient->Service()==IOClient::NONE ) SleepThisThread( 0 );
					else server->_pendingClient = NULL;
					hasPending = true;
				}
				LeaveCriticalSection( &server->_pendingLock );
			}
			if( !hasPending )
			{
#if DEBUG_MUTEX
				EnterCriticalSection( &server->_clientLock , 100 , "MultiStreamIOServer::_IOThread (_clientLock #1)" );
#else // !DEBUG_MUTEX
				EnterCriticalSection( &server->_clientLock );
#endif // DEBUG_MUTEX
				bool ioDone = false;
				for( int i=0 ; i<clients.size() && !ioDone ; i++ )
				{
					idx = (idx+1)%clients.size();
					switch( clients[idx]->Service() )
					{
					case IOClient::COMPLETE:
						clients[idx]->SetServer( NULL );
						clients[idx] = clients[clients.size()-1];
						clients.pop_back();
					case IOClient::SUCCESS:
						ioDone = true;
						break;
					}
				}
				LeaveCriticalSection( &server->_clientLock );
				if( !ioDone ) SleepThisThread( 1 );
				else SleepThisThread( 0 );
			}
		}
	}
	return 0;
}
