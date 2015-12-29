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

#ifndef SOCKET_INCLUDED
#define SOCKET_INCLUDED

#include <stdarg.h>
#include "XPlatform.h"
#include <stdio.h>
#include "Util/Time.h"
#include "Util/Array.h"

void printfId( const char* format , ... );
void fprintfId( FILE* fp , const char* format , ... );

class ConnectionData
{
public:
	EndpointAddress localAddr , peerAddr;
	int localPort , peerPort;
};


template< class C > bool ReceiveOnSocket( Socket& s ,      Pointer( C ) data , size_t dataSize );
template< class C > bool SendOnSocket   ( Socket& s , ConstPointer( C ) data , size_t dataSize );
template< class C > void ReceiveOnSocket( Socket& s ,      Pointer( C ) data , size_t dataSize , const char* errorMessage , ... );
template< class C > void SendOnSocket   ( Socket& s , ConstPointer( C ) data , size_t dataSize , const char* errorMessage , ... );

AcceptorSocket GetListenSocket( int& port );
Socket AcceptSocket( AcceptorSocket listen );
Socket GetConnectSocket( const char* address , int port , int ms=5 , bool progress=false );
Socket GetConnectSocket( EndpointAddress , int port , int ms=5 , bool progress=false );
void CloseSocket( Socket& s );
void CloseAcceptorSocket( AcceptorSocket& s );
EndpointAddress GetLocalSocketEndpointAddress( Socket& s );
int             GetLocalSocketPort           ( Socket& s );
EndpointAddress GetLocalSocketEndpointAddress( Socket& s );
int             GetPeerSocketPort            ( Socket& s );
bool GetHostAddress( char* address , const char* prefix = NULL );
bool GetHostEndpointAddress( EndpointAddress* address , const char* prefix=NULL );
void PrintHostAddress( void );


class DataStream
{
public:
	virtual ~DataStream( void ){ }
	virtual bool write( ConstPointer( byte ) buf , int len ) = 0;
	virtual bool read ( Pointer(      byte ) buf , int len ) = 0;

	static DataStream* GetDataStream( Socket sock , bool master , bool cleanUp = true );
};

class DataStreamConstructor
{
	Socket _sock;
	bool _master;
	bool _cleanUp;
	DataStream* _stream;
	EndpointAddress _myAddr;
	int _myPID;

public:
	static const int STEPS = 3;
	DataStreamConstructor( void );
	void init( Socket sock , bool master , bool cleanUp = true );
	void doStep( int sNum );
	DataStream* getDataStream( void );
};
class SocketStream : public DataStream
{
	Socket _sock;
	bool _mySocket;
public:
	SocketStream( Socket sock );
	SocketStream( EndpointAddress address , int port , int ms=5 , bool progress=false );
	SocketStream( const char* address , int port , int ms=5 , bool progress=false );
	~SocketStream( void );
	bool write( ConstPointer( byte ) buf , int len );
	bool read ( Pointer(      byte ) buf , int len );
};

class SharedMemoryBuffer
{
protected:
	Pointer( byte ) _buf1;
	Pointer( byte ) _buf2;
	int _bufSize1 , _bufSize2;
	Signal _readyForWriting1 , _readyForWriting2;
	SharedMemoryBuffer( void );
	~SharedMemoryBuffer( void );
public:
	class FirstStream : public DataStream
	{
		friend class SharedMemoryBuffer;
	protected:
		SharedMemoryBuffer* _smb;
		FirstStream( SharedMemoryBuffer* smb );
	public:
		~FirstStream( void );
		bool write( ConstPointer( byte ) buf , int len );
		bool read ( Pointer(      byte ) buf , int len );
	};
	class SecondStream : public DataStream
	{
		friend class SharedMemoryBuffer;
	protected:
		SharedMemoryBuffer* _smb;
		SecondStream( SharedMemoryBuffer* smb );
	public:
		~SecondStream( void );
		bool write( ConstPointer( byte ) buf , int len );
		bool read (      Pointer( byte ) buf , int len );
	};
	class StreamPair
	{
	public:
		StreamPair( void );
		FirstStream* first;
		SecondStream* second;
		static bool CreateSharedBufferPair( StreamPair& pair );
	};
};

#include "Socket.inl"
#endif // SOCKET_INCLUDED
