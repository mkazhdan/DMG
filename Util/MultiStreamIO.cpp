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
#include "MultiStreamIO.h"
#ifdef _WIN32
#include <atlstr.h>
#endif // _WIN32

////////////////////////////
// VariableIOClientStream //
////////////////////////////
VariableIOClientStream::~VariableIOClientStream( void )
{
	if( _fp ) close( );
	FreePointer( _data );
}
VariableIOClientStream::VariableIOClientStream( void ) : IOClient( )
{
	_fp = NULL;
}
int VariableIOClientStream::close( void )
{
	if( !_fp )
	{
		fprintf( stderr , "Closing null file pointer!\n" );
		return 0;
	}
	{
#if DEBUG_MUTEX
		EnterCriticalSection( &lock , 100 , "VariableIOClientStream::close" );
#else // !DEBUG_MUTEX
		EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
		_endIO = true;
		LeaveCriticalSection( &lock );
	}
	while( server ) SleepThisThread(1);
	if( !_read )
	{
		while( _frontAndBack<_current )
		{
			size_t ioBytes = fwrite( _data + ( _frontAndBack % _bufferMultiplier ) * _blockSize , 1 ,_blockSize , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
			_frontAndBack++;
		}
		if( _head<_pseudoHead )
		{
			size_t ioBytes = fwrite( _data + ( _current % _bufferMultiplier ) * _blockSize , 1 , _pseudoHead-_head , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
		}
	}
	int ret = fclose( _fp );
	_fp = NULL;
	return ret;
}
void VariableIOClientStream::Initialize( const char* fileName , int maxIOSize , int bufferMultiplier , bool read )
{
#if DEBUG_MUTEX
	EnterCriticalSection( &lock , 100 , "VariableIOClientStream::Initialize" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX

	if( read ) _fp = fopen( fileName , "rb" );
	else       _fp = fopen( fileName , "wb" );

	_bufferMultiplier = bufferMultiplier;

	_head = _pseudoHead = 0;
	_read = read;
	_blockSize = maxIOSize;
	while( _blockSize<IO_BLOCK_SIZE ) _blockSize += maxIOSize;
	_data = AllocPointer< byte >( _blockSize*_bufferMultiplier );
	if( !_data )
	{
		fprintf(stderr , "Failed to allocate memory for BufferedIOState buffer: %d x %d\n" , _blockSize , _bufferMultiplier );
		exit(0);
	}
	_endIO = false;
	_current      = 0;
	_frontAndBack = 0;
	if( _read )
	{
		int ioBytes = int( fread( _data , 1 , _blockSize , _fp ) );
		ReadBytes += ioBytes;
		_head      = ioBytes;
		if( ioBytes<_blockSize ) _endIO = true;
		_frontAndBack++;
	}
	LeaveCriticalSection( &lock );
}

// Assumes that whatever data wants to be read has already been brought into the buffer.
size_t VariableIOClientStream::read( Pointer( byte ) buffer , size_t sz )
{
	if( sz > _blockSize ) fprintf( stderr , "Requesting a read that is too large: %zd\n" , sz );
#if DEBUG_MUTEX
	EnterCriticalSection( &lock , 100 , "VariableIOClientStream::read" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
	long long start , end;
	int startBlock , endBlock;
	start = _pseudoHead;
	end   = start + sz;
	if( end>_head && _endIO ) end = _head;
	if( (start/_blockSize) != (end/_blockSize) && !server ) Service();		// If we need data from the next block and there is no server
																			// go ahead and read it.

	startBlock = int(  start  / _blockSize );
	endBlock   = int( (end-1) / _blockSize );
	int startIdx = startBlock % _bufferMultiplier;
	int   endIdx =   endBlock % _bufferMultiplier;
	int startOffset =    start % _blockSize;
	int   endOffset =( (end-1) % _blockSize ) + 1;
	if( startIdx == endIdx ) memcpy( buffer , _data + startIdx*_blockSize + startOffset , endOffset-startOffset );
	else
	{
		memcpy( buffer , _data + startIdx*_blockSize + startOffset , _blockSize-startOffset );
		memcpy( (void*)( size_t( GetAddress( buffer ) ) + _blockSize-startOffset ) , _data + endIdx*_blockSize , endOffset );
	}
	_pseudoHead += end-start;
	if( (start/_blockSize) != (end/_blockSize) ) _current++;
	LeaveCriticalSection( &lock );
	return end-start;
}
// Assumes that there is room in the buffer for storing the data.
size_t VariableIOClientStream::write( ConstPointer( byte ) buffer , size_t sz )
{
	if( sz > _blockSize ) fprintf( stderr , "Requesting a write that is too large: %zd\n" , sz );
#if DEBUG_MUTEX
	EnterCriticalSection( &lock , 100 , "VariableIOClientStream::write" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
	long long start , end;
	int startBlock , endBlock;
	start = _pseudoHead;
	end   = start + sz;
	if( (start/_blockSize) != (end/_blockSize) && !server ) Service();		// If we need data from the next block and there is no server
																			// go ahead and read it.

	startBlock = int(  start  / _blockSize );
	endBlock   = int( (end-1) / _blockSize );
	int startIdx = startBlock % _bufferMultiplier;
	int   endIdx =   endBlock % _bufferMultiplier;
	int startOffset = start % _blockSize;
	int   endOffset =( (end-1) % _blockSize ) + 1;
	if( startIdx == endIdx )
	{
		memcpy( _data + startIdx*_blockSize + startOffset , buffer , endOffset-startOffset );
	}
	else
	{
		memcpy( _data + startIdx*_blockSize + startOffset , buffer , _blockSize-startOffset );
		memcpy( _data + endIdx*_blockSize , (void*)( size_t(buffer) + _blockSize-startOffset ) , endOffset );
	}
	_pseudoHead += end-start;
	if( (start/_blockSize) != (end/_blockSize) ) _current++;
	LeaveCriticalSection( &lock );
	return end-start;
}

int VariableIOClientStream::Service( void )
{
	size_t ioBytes;
	// Try and grab the lock:
	//   [W] _head
	//   [W] _frontAndBack
	//   [W] _endIO 
	//   [R] _data
#if DEBUG_MUTEX
	EnterCriticalSection( &lock , 100 , "VariableIOClientStream::Service" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX

	// If we are reading, see if we can advance the front pointer
	if( _read )
		if( !_endIO && _frontAndBack-_current<_bufferMultiplier )
		{
			ioBytes = fread( _data + ( _frontAndBack % _bufferMultiplier ) * _blockSize , 1 ,_blockSize , _fp );
			if( ioBytes != _blockSize ) _endIO = true;
			ReadBytes += ioBytes;
			_head += ioBytes;
			_frontAndBack++;
			LeaveCriticalSection( &lock );
			if( _endIO ) return COMPLETE;
			else		 return SUCCESS;
		}
	// If we are writing, see if we can flush the back pointer
	if( !_read )
		if( _frontAndBack<_current )
		{
			ioBytes = fwrite( _data + ( _frontAndBack % _bufferMultiplier ) * _blockSize , 1 ,_blockSize , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
			_frontAndBack++;
			LeaveCriticalSection( &lock );
			return SUCCESS;
		}
		else if( _endIO && _head<_pseudoHead )
		{
			ioBytes = fwrite( _data + ( _current % _bufferMultiplier ) * _blockSize , 1 ,_pseudoHead-_head , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
			LeaveCriticalSection( &lock );
			return COMPLETE;
		}
	if( _endIO )
	{
		LeaveCriticalSection( &lock );
		return COMPLETE;
	}
	else
	{
		LeaveCriticalSection( &lock );
		return NONE;
	}
}
/////////////////////////
// FixedIOClientStream //
/////////////////////////

FixedIOClientStream::FixedIOClientStream( void ) : IOClient( )
{
	hFile=NULL;
	off = NullPointer< int >( );
	data = NullPointer< byte >( );
	r=rs=win=b=0;
	blockSize=0;
}
FixedIOClientStream::~FixedIOClientStream( void )
{
	FreePointer( off );
	FreePointer( data );
}
void FixedIOClientStream::Init( FileHandle hFile , int rowSize , int rows , int bufferMultiplier )
{
	FreePointer( off );
	FreePointer( data );
	this->hFile = hFile;
	r = rows;
	rs = rowSize;
	b = bufferMultiplier;
	off = AllocPointer< int >( b );
	if( !off )
	{
		fprintf( stderr , "Failed to allocate memory for offsets\n" );
		exit(0);
	}
}
void FixedIOClientStream::Reset( bool read , int minWindowSize )
{
	this->read=read;
	win=minWindowSize<r ? minWindowSize : r;
	while( win*rs<IO_BLOCK_SIZE && win<r ) win++;

	blockSize=((win*rs+(BYTES_PER_SECTOR-1))/BYTES_PER_SECTOR+1)*BYTES_PER_SECTOR;
	if( !((win*rs)%BYTES_PER_SECTOR) ) blockSize -= BYTES_PER_SECTOR;
	FreePointer( data );
	data = AllocPointer< byte >( blockSize * b );
	if( !data )
	{
		fprintf( stderr , "Failed to allocate memory for StreamState buffer: %d x %d\n" ,  blockSize , b );
		exit(0);
	}
	if( read )
	{
		long long ioBytes;
		SeekFileHandle( hFile , 0 );
		ioBytes = ReadFileHandle( hFile , data , blockSize );
		ReadBytes += ioBytes;
	}
	current = 0;
	back    = 0;
	front   = win;
	off[0]  = 0;
}
void FixedIOClientStream::Unset(void)
{
	read=false;
	win=0;
	blockSize=0;
	FreePointer( data );
	current=0;
	back=0;
	front=win;
	off[0]=0;
}
Pointer( byte ) FixedIOClientStream::operator[]	(int idx)
{
	Pointer( byte ) rowData;
#if DEBUG_MUTEX
	EnterCriticalSection( &lock , 100 , "FixedIOClientStream::operator[]" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
	int bIndex = (idx/win)%b;		// Which block is it in?
	int wIndex =  idx%win;			// Which row of the block?
	rowData = data + bIndex*blockSize + off[bIndex] + wIndex*rs;
	LeaveCriticalSection( &lock );
	return rowData;
}

bool FixedIOClientStream::Advance( void )
{
#if DEBUG_MUTEX
	EnterCriticalSection( &lock , 100 , "FixedIOClientStream::Advance" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
	if( current+1>=back && current+1<front )
	{
		current++;
		LeaveCriticalSection( &lock );
		return true;
	}
	else if( current+1>=r )
	{
		current = r;
		LeaveCriticalSection( &lock );
		return true;
	}
	else
	{
		LeaveCriticalSection( &lock );
		return false;
	}
}
int FixedIOClientStream::Service(void)
{
	long long ioBytes;
	int ioRows=win;
	int ioBlockSize=blockSize;
	// Try and grab the lock
#if DEBUG_MUTEX
	EnterCriticalSection( &lock , 100 , "FixedIOClientStream::Service (1)" );
#else // !DEBUG_MUTEX
	EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
	// First see if we can advance the front pointer
	if(front<r && front+win-back <= win*b)
	{
		long long locationOnDisk=(long long)(front)*rs;
		long long readStart=(locationOnDisk/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
		int offset=int( locationOnDisk-readStart );
		int bIndex=(front/win)%b;
		if(read)
		{
			LeaveCriticalSection( &lock );
			if(front+win<=r)	ioRows=win;
			else				ioRows=r-front;
			ioBlockSize=((offset+ioRows*rs+BYTES_PER_SECTOR-1)/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
			SeekFileHandle( hFile , readStart );
			ioBytes = ReadFileHandle( hFile , data + bIndex*blockSize , ioBlockSize );
			ReadBytes+=ioBytes;
#if DEBUG_MUTEX
			EnterCriticalSection( &lock , 100 , "FixedIOClientStream::Service (2)" );
#else // !DEBUG_MUTEX
			EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
		}
		off[bIndex]=offset;
		front+=ioRows;
		LeaveCriticalSection( &lock );
		return SUCCESS;
	}
	// Now try to free up trailing memory
	else if( (back+win<=current || current>=r) && back<r)	// If we won't write out needed data and there is what to write
	{
		long long locationOnDisk=(long long)(back)*rs;
		long long writeStart=(locationOnDisk/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
		int offset=int( locationOnDisk-writeStart );
		int bIndex=(back/win)%b;
		if(!read)											// If we are doing a write, write out the data
		{
			LeaveCriticalSection( &lock );
			if(back+win<=r)	ioRows=win;
			else			ioRows=r-back;
			ioBlockSize=((offset+ioRows*rs+BYTES_PER_SECTOR-1)/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
			SeekFileHandle( hFile , writeStart );
			ioBytes = WriteFileHandle( hFile , data + bIndex*blockSize , ioBlockSize );
			WriteBytes+=ioBytes;
#if DEBUG_MUTEX
			EnterCriticalSection( &lock , 100 , "FixedIOClientStream::Service (3)" );
#else // !DEBUG_MUTEX
			EnterCriticalSection( &lock );
#endif // DEBUG_MUTEX
		}
		back+=ioRows;
		if(!read && back<r)
		{
			long long locationOnDisk=(long long)(back)*rs;
			long long writeStart=(locationOnDisk/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
			int offset=int( locationOnDisk-writeStart );		// The number of bytes that need to be copied over from the previous buffer
			int bIndex=(back/win)%b;
			int oldBIndex=(bIndex+b-1)%b;
			if( offset ) memcpy( data + bIndex*blockSize , data + oldBIndex*blockSize + ioBlockSize-BYTES_PER_SECTOR , offset );
		}
		LeaveCriticalSection( &lock );
		return SUCCESS;
	}
	// Check if we are done
	else if( back>=r )
	{
		LeaveCriticalSection( &lock );
		return COMPLETE;
	}
	else
	{
		LeaveCriticalSection( &lock );
		return NONE;
	}
}
/////////////////////////
// MultiStreamIOClient //
/////////////////////////
MultiStreamIOClient::MultiStreamIOClient( const char* fileName , int rs , int r , int bufferMultiplier , bool writeOnly )
{
	hFile = CreateFileHandle( fileName );
	if( !hFile ) fprintf( stderr , "[ERROR] Failed to create file handle\n" ) , exit(0);
	if( writeOnly )
	{
		// Pre-allocate file space
		long long fileSize;
		fileSize=(long long)(r)*rs;
		fileSize=((fileSize+IOClient::BYTES_PER_SECTOR-1)/IOClient::BYTES_PER_SECTOR)*IOClient::BYTES_PER_SECTOR;
		if( !SeekFileHandle( hFile , fileSize ) ) fprintf( stderr , "[ERROR] Failed to seek file\n" ) , exit( 0 );
		if( !SetEndOfFileHandle( hFile ) ) fprintf( stderr , "[WARNING] Failed to set end of file\n"  );
	}
	stream.Init( hFile , rs , r , bufferMultiplier );
	server = NULL;
}
MultiStreamIOClient::MultiStreamIOClient( int rs , int r , int bufferMultiplier , const char* prefix , bool deleteOnClose )
{
	hFile = CreateFileHandle( NULL , prefix , deleteOnClose );
	if( !hFile ) fprintf( stderr , "[ERROR] Failed to create file handle\n" ) , exit( 0 );

	// Pre-allocate file space
	long long fileSize;
	fileSize = (long long)(r)*rs;
	fileSize = ((fileSize+IOClient::BYTES_PER_SECTOR-1)/IOClient::BYTES_PER_SECTOR)*IOClient::BYTES_PER_SECTOR;
	if( !SeekFileHandle( hFile , fileSize ) ) fprintf( stderr , "[ERROR] Failed to seek file\n" ) , exit( 0 );
	if( !SetEndOfFileHandle( hFile ) ) fprintf( stderr , "[WARNING] Failed to set end of file\n"  );
	stream.Init( hFile , rs , r , bufferMultiplier );
	server=NULL;
}
MultiStreamIOClient::MultiStreamIOClient( int rs , int r , int bufferMultiplier , const char* dir , const char* prefix , bool deleteOnClose )
{
	hFile = CreateFileHandle( dir , prefix , deleteOnClose );
	if( !hFile ) fprintf( stderr , "[ERROR] Failed to create file handle\n" ) , exit( 0 );

	// Pre-allocate file space
	long long fileSize;
	fileSize=(long long)(r)*rs;
	fileSize=((fileSize+IOClient::BYTES_PER_SECTOR-1)/IOClient::BYTES_PER_SECTOR)*IOClient::BYTES_PER_SECTOR;
	if( !SeekFileHandle( hFile , fileSize ) ) fprintf( stderr , "[ERROR] Failed to seek file\n" ) , exit( 0 );
	if( !SetEndOfFileHandle( hFile ) ) fprintf( stderr , "[WARNING] Failed to set end of file\n"  );
	stream.Init( hFile , rs , r , bufferMultiplier );
	server=NULL;
}
MultiStreamIOClient::~MultiStreamIOClient(void)
{
	finalize();
}
void MultiStreamIOClient::finalize(void)
{
	CloseFileHandle( hFile );
	hFile=NULL;
}
int MultiStreamIOClient::Service		( void )								{ return stream.Service(); }
int		MultiStreamIOClient::rows		( void )						const	{ return stream.r;  }
int		MultiStreamIOClient::rowSize	( void )						const	{ return stream.rs; }
void	MultiStreamIOClient::reset		( bool r , int minWindowSize )			{ stream.Reset( r , minWindowSize ); }
void	MultiStreamIOClient::unset		( void )								{ stream.Unset(); }
Pointer( byte ) MultiStreamIOClient::operator []( int idx ) { return stream[idx]; }
bool	MultiStreamIOClient::Advance	( void )								{ return stream.Advance(); }
void	MultiStreamIOClient::advance	( void )
{
	// BADNESS!!! Server may not be NULL if it was set for an earlier server and not for the newer one.
	if( server )
	{
		while(1)
			if( Advance() )	return;
			else server->SetPending( this );
	}
	else
		if( !Advance() )
		{
			int serviceState = Service();
			while( serviceState==IOClient::SUCCESS ) serviceState = Service();
			if( !Advance() ) fprintf( stderr , "[ERROR] MultiStreamIOClient::advance: failed to advance: %d %d %d\n" , stream.current , stream.front , stream.r ) , exit(0);
		}
		else Service();	// To make sure that the last row gets written
}
void MultiStreamIOClient::SetServer( MultiStreamIOServer* server ) { IOClient::SetServer( server ); }