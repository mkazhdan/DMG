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
#ifndef MULTI_STREAM_IO
#define MULTI_STREAM_IO
#include <vector>
#include "BaseMultiStreamIO.h"
#include "GridStream.h"

// [WARNING] Apparently bad things can happen if the bufferMultiplier is set to one. Possibly related to the fact that rows can overflow buffers? 



// When only the  maximum size of the requested read/write is known in advance.
class VariableIOClientStream : public IOClient
{
	int		_blockSize;							// The actual size of a read/write operation
	Pointer( byte ) _data;
	int		_bufferMultiplier;					// The number of blocks to store in memory at any given time (must be at least 2)
	FILE*	_fp;								// Where to read/write data
	int		_current;							// The index of the buffer from which data can be requested
	int		_frontAndBack;						// For reading, [ _current , _frontAndBack ] is in the buffer.
												// For writing, [ _frontAndBack , _current ]  is in the buffer
	bool	_read;								// Are we reading data?
	bool	_endIO;								// Have we performed all the requisite I/O?
	long long _pseudoHead;						// Where are we believed to be in the file
	long long _head;							// Where are we actually in the file
public:
	VariableIOClientStream( void );
	void Initialize( const char* fileName , int maxIOSize , int bufferMultiplier , bool read );
	~VariableIOClientStream( void );
	int		Service	( void );
	size_t  read (      Pointer( byte ) buffer , size_t sz );
	size_t  write( ConstPointer( byte ) buffer , size_t sz );
	int     close( void );
};

class FixedIOClientStream : public IOClient
{
public:
	int					blockSize;
	Pointer( int )		off;
	Pointer( byte )		data;
	int					r , rs , win , b;
	FileHandle			hFile;
	bool				read;
	int					current,front,back;

	FixedIOClientStream ( void );
	~FixedIOClientStream( void );
	void	Init	( FileHandle hFile , int rowSize , int rows , int bufferMultiplier );
	void	Reset	( bool read , int minWindowSize );
	void	Unset	( void );
	int		Service	( void );
	bool	Advance	( void );
	Pointer( byte ) operator[] ( int idx );
};
class MultiStreamIOClient : public StreamingGrid
{
	FileHandle hFile;
	FixedIOClientStream stream;
public:
	MultiStreamIOClient	( const char* fileName , int rs , int r , int bufferMultiplier , bool writeOnly=false );
	MultiStreamIOClient	( int rs , int r , int bufferMultiplier , const char* prefix , bool deleteOnClose );
	MultiStreamIOClient	( int rs , int r , int bufferMultiplier , const char* dir , const char* prefix , bool deleteOnClose);
	~MultiStreamIOClient(void);
	void	finalize	(void);

	int		rows		( void ) const;
	int		rowSize		( void ) const;
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
	void	unset		( void );
	Pointer( byte ) operator[] ( int idx );
	void	SetServer	( MultiStreamIOServer* server );
	bool	Advance		( void );
	int		Service		( void );
};

#endif // MULTI_STREAM_IO