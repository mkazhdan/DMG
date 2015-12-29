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
#ifndef GRID_STREAM_INCLUDED
#define GRID_STREAM_INCLUDED

#include <Util/Array.h>
#include "BaseMultiStreamIO.h"
#include "Socket.h"

class StreamingGrid : public IOClient
{
public:
	virtual			~StreamingGrid	( void )							{;}
	virtual int		rows			( void ) const						= 0;
	virtual int		rowSize			( void ) const						= 0;
	virtual Pointer( byte ) operator[] ( int idx )						= 0;
	virtual void	advance			( void )							{;}
	virtual void	reset			( bool read , int minWindowSize )	{;}
	virtual void	unset			( void )							{;}
	virtual void	SetServer		( MultiStreamIOServer* server )		{;}
	virtual int		Service			( void )							{ return IOClient::NONE; }
};
class NULLStreamingGrid : public StreamingGrid
{
	Pointer( byte ) _data;
	int _rows , _rowSize;
public:
	NULLStreamingGrid	( int rowSize , int rows );
	~NULLStreamingGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
};

class BufferedStreamingGrid : public StreamingGrid
{
	bool read;
	Pointer( byte ) data;
	int current , win;
	StreamingGrid* sg;
public:
	BufferedStreamingGrid	(StreamingGrid* sg);
	~BufferedStreamingGrid	(void);
	int		rows			(void) const;
	int		rowSize			(void) const;
	Pointer( byte )			operator[] ( int idx );
	void	advance			(void);
	void	reset			(bool read,int minWindowSize);
};
template< class Real , int Channels >
class MultiSocketBackedGrid : public StreamingGrid
{
	bool _read , _readComplete;
	Pointer( Real ) _data;
	Pointer( Pointer( Real ) ) _subData;
	int		_rows , *_rowSizes , _rowSize , _current , _sockCount;
	Socket*	_socks;
public:
	MultiSocketBackedGrid( Socket* socks , int *rowSizes , int sockCount , int rows );
	~MultiSocketBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
};

class SocketBackedGrid : public StreamingGrid
{
	bool _read , _readComplete;
	Pointer( byte ) _data;
	int		_r , _rs , _c;
	Socket	_sock;
public:
	SocketBackedGrid	( Socket sock , int rowSize , int rows );
	~SocketBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer ( byte ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
};
class SharedMemoryGrid : public StreamingGrid
{
	bool    _read , _dataReady;
	Pointer( byte ) _data;
	int     _r , _rs , _c;
	Signal *_signals;
	bool _signalState;
	int     _sCount;
public:
	SharedMemoryGrid	( Signal* signals , int signalCount , bool signalState , Pointer( byte ) row , int rowSize , int rows );
	~SharedMemoryGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
};
class MemoryBackedGrid : public StreamingGrid
{
	Pointer( byte ) _data;
	int		_r,_rs;
	bool	_del;
public:
	MemoryBackedGrid	( Pointer( byte ) data , int rowSize , int rows );
	MemoryBackedGrid	( int rowSize , int rows );
	~MemoryBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
};
class FileBackedGrid : public StreamingGrid
{
	long long _rSize,_wSize;
	FileHandle _fh;
	Pointer( byte ) _data;
	int		_r,_rs,_c,_w;
	bool	_read;
	void	_readNext(void);
	void	_writeNext(void);
public:
	FileBackedGrid		( int rowSize , int rows , FileHandle fh );
	FileBackedGrid		( int rowSize , int rows );
	~FileBackedGrid		(void);
	int		rows		(void) const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
	void	reset		(bool read,int minWindowSize);
};
class CompressedFileBackedGrid : public StreamingGrid
{
	int		_sSize;
	long long _rSize,_wSize;
	FileHandle _fh;
	Pointer( byte ) _data;
	Pointer( byte ) _scratch;
	int		_r,_rs,_c,_w;
	bool	_read;
	void	_readNext(void);
	void	_writeNext(void);
public:
	CompressedFileBackedGrid	( int rowSize , int rows , FileHandle fh );
	CompressedFileBackedGrid	( int rowSize , int rows );
	~CompressedFileBackedGrid	(void);
	int		rows				(void) const;
	int		rowSize				(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance				(void);
	void	reset				(bool read,int minWindowSize);
};

#include "GridStream.inl"
#endif // GRID_STREAM_INCLUDED