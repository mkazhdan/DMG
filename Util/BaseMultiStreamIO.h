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
#ifndef BASE_MULTI_STREAM_IO
#define BASE_MULTI_STREAM_IO
#include "XPlatform.h"
#include <vector>

// [WARNING] Apparently bad things can happen if the bufferMultiplier is set to one. Possibly related to the fact that rows can overflow buffers? 

class IOClient
{
public:
	static const int BYTES_PER_SECTOR;			// The atomic size of a read/write operation (for some IO applications)
	static const int IO_BLOCK_SIZE;				// A nice size for reads
	static long long ReadBytes , WriteBytes;	// Total bytes read/written
	enum
	{
		NONE,		// If there was no I/O that the client could do
		SUCCESS,	// If the client succeeded in performing the I/O
		COMPLETE,	// If the client won't need to do any more I/O
	};

	friend class MultiStreamIOServer;
	CriticalSectionLock lock;				// A locking devices so that the server doesn't step on the client's toes
	class MultiStreamIOServer* server;		// The server responsible for processing the job requests (this should only be set by the server)!

	IOClient ( void );
	~IOClient( void );

	void SetServer( class MultiStreamIOServer* server );
	virtual int	Service( void ) = 0;		// The task the server asks the client to do
};

class MultiStreamIOServer
{
	ThreadHandle _ioThread;
	static int _IOThread( void* );
	std::vector< IOClient* > _clients;
	CriticalSectionLock _pendingLock , _clientLock;
	IOClient* _pendingClient;
	bool _allowNewClients;
public:
	MultiStreamIOServer ( void );
	~MultiStreamIOServer( void );

	int clientNum( void );
	virtual bool SetPending( IOClient* client );
	virtual void AddClient ( IOClient* client );
};
#endif // BASE_MULTI_STREAM_IO