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
#include "GridStream.h"
#include "Time.h"

///////////////////////
// NULLStreamingGrid //
///////////////////////
NULLStreamingGrid::NULLStreamingGrid( int rowSize , int rows )
{
	_rows = rows;
	_rowSize = rowSize;
	_data = AllocPointer< byte >( rowSize );
	memset( _data , 0 , rowSize );
}
NULLStreamingGrid::~NULLStreamingGrid( void ){ FreePointer( _data ); }

int NULLStreamingGrid::rows( void ) const { return _rows; }
int NULLStreamingGrid::rowSize( void ) const { return _rowSize; }
Pointer( byte ) NULLStreamingGrid::operator [] ( int idx ){ return _data; }

///////////////////////////
// BufferedStreamingGrid //
///////////////////////////
BufferedStreamingGrid::BufferedStreamingGrid(StreamingGrid* sg)
{
	data = NullPointer< byte >( );
	current=0;
	this->sg=sg;
}
BufferedStreamingGrid::~BufferedStreamingGrid(void)
{
	FreePointer( data );
	sg=NULL;
}
int BufferedStreamingGrid::rows		(void) const	{return sg->rows();}
int BufferedStreamingGrid::rowSize	(void) const	{return sg->rowSize();}
void BufferedStreamingGrid::reset	(bool read,int minWindowSize)
{
	this->read=read;
	win=minWindowSize<rows() ? minWindowSize : rows();

	FreePointer( data );
	data = AllocPointer< byte >( rowSize() * win );
	if( !data ) fprintf( stderr , "[ERROR] Failed to allocate memory for BufferedStreamingGrid\n" ) , exit( 0 );
	sg->reset( read , 1 );
	if( read )
	{
		for(int w=0;w<win;w++)
		{
			Pointer( byte ) row = (*sg)[w];
			memcpy( data+w*rowSize() , row , rowSize() );
			sg->advance();
		}
	}
	current=0;
}
void BufferedStreamingGrid::advance(void)
{
	if(read)
	{
		current++;
		if(current+win-1<rows())
		{
			sg->advance();
			memcpy( data + ((current+win-1)%win)*rowSize() , (*sg)[current+win-1] , rowSize() );
		}
	}
	else
	{
		if(current-win+1>=0)
		{
			memcpy( (*sg)[current-win+1] , data+((current-win+1)%win)*rowSize() , rowSize() );
			sg->advance();
		}
		current++;
	}
}
Pointer( byte ) BufferedStreamingGrid::operator[]	(int idx)
{
	return data+(idx%win)*rowSize();
}
//////////////////////
// MemoryBackedGrid //
//////////////////////
MemoryBackedGrid::MemoryBackedGrid( Pointer( byte ) data , int rs , int r )
{
	_del=false;
	_r=r, _rs=rs;
	_data=data;
}
MemoryBackedGrid::MemoryBackedGrid( int rs , int r )
{
	_del=true;
	_r=r, _rs=rs;
	_data = AllocPointer< byte >( r * rs );
	if( !_data )	fprintf(stderr,"Failed to allocate memory backed grid of size %d x %d\n" , rs , r ) , exit(0);
}
MemoryBackedGrid::~MemoryBackedGrid(void)
{
	if( _del ) FreePointer( _data );
	_r=_rs=0;
}
int	MemoryBackedGrid::rows(void) const
{
	return _r;
}
int MemoryBackedGrid::rowSize(void) const
{
	return _rs;
}
Pointer( byte ) MemoryBackedGrid::operator[]( int idx )
{
	return _data + _rs * idx;
}
//////////////////////
// SocketBackedGrid //
//////////////////////
SocketBackedGrid::SocketBackedGrid( Socket sock , int rs , int r )
{
	_sock = sock;
	_read = true;
	_r = r;
	_rs = rs;
	_data = AllocPointer< byte >( rs );
	if( !_data ) fprintfId( stderr , "Failed to allocate\n" ) , exit(0);
}
SocketBackedGrid::~SocketBackedGrid( void )
{
	FreePointer( _data );
	_data = NullPointer< byte >( );
	_r = _rs = 0;
}
int	SocketBackedGrid::rows( void ) const { return _r; }
int SocketBackedGrid::rowSize( void ) const { return _rs; }
Pointer( byte ) SocketBackedGrid::operator[]( int idx )
{
	if( _read && !_readComplete ) ReceiveOnSocket( _sock , _data , _rs , "SocketBackedGrid::operator[]" ) , _readComplete = true;
	return _data;
}
void SocketBackedGrid::advance( void )
{
	if( !_read )
	{
		SendOnSocket( _sock , ( ConstPointer( byte ) )_data , _rs , "SocketBackedGrid::advance" );
	}
	else
	{
		if( !_readComplete ) ReceiveOnSocket( _sock , _data , _rs , "SocketBackedGrid::advance" );
		_readComplete = false;
	}
	_c++;
}
void SocketBackedGrid::reset( bool read , int minWindowSize )
{
	_read = read;
	_readComplete = !_read;
	_c = 0;
}

//////////////////////
// SharedMemoryGrid // 
//////////////////////
SharedMemoryGrid::SharedMemoryGrid( Signal* signals , int signalCount , bool signalState , Pointer( byte ) row , int rowSize , int rows )
{
	_read = true;
	_r = rows;
	_rs = rowSize;
	_data = row;
	_sCount = signalCount;
	_signals = signals;
	_signalState = signalState;
}
SharedMemoryGrid::~SharedMemoryGrid(void)
{
	_data = NullPointer< byte >( );
	_r = _rs = 0;
}
int	SharedMemoryGrid::rows(void) const { return _r; }
int SharedMemoryGrid::rowSize(void) const { return _rs; }
Pointer( byte ) SharedMemoryGrid::operator[]( int idx )
{
	if( !_dataReady ) for( int i=0 ; i<_sCount ; i++ ) Signaller::Start( _signals[i] , _signalState );
	_dataReady = true;
	return _data;
}
void SharedMemoryGrid::advance( void )
{
	for( int i=0 ; i<_sCount ; i++ ) Signaller::Stop( _signals[i] , _signalState );

	_dataReady = false;
	_c++;
}
void SharedMemoryGrid::reset( bool read , int minWindowSize )
{
	_read = read;
	_dataReady = false;
	_c = 0;
}


////////////////////
// FileBackedGrid //
////////////////////
FileBackedGrid::FileBackedGrid( int rs , int r )
{
	_rSize=_wSize=0;
	_fh = CreateFileHandle( "." , NULL , true );
	_read=true;
	_r=r, _rs=rs;
	_data = NullPointer< byte >( );
	_c=_w=0;
}
FileBackedGrid::FileBackedGrid( int rs , int r , FileHandle fh )
{
	_rSize=_wSize=0;
	_fh = fh;
	_read=true;
	_r=r, _rs=rs;
	_data = NullPointer< byte >( );
	_c=_w=0;
}
FileBackedGrid::~FileBackedGrid(void)
{
	if( !_read ) for( int i=0 ; i<_w ; i++ ) _writeNext();
	FreePointer( _data );
	_r=_rs=_c=_w=0;
	CloseFileHandle( _fh );
	printf( "Read/Write Size: %lld MB/ %lld MB\n" , _rSize>>20 , _wSize>>20 ) , fflush( stdout );
}
int	FileBackedGrid::rows(void) const	{	return _r;	}
int FileBackedGrid::rowSize(void) const	{	return _rs;	}
Pointer( byte ) FileBackedGrid::operator[]( int idx )
{
	return _data + _rs * (idx%_w);
}
void FileBackedGrid::reset(bool read,int minWindowSize)
{
	if(!_read)	for(int i=0;i<_w;i++)	_writeNext();
	FreePointer( _data );

	SeekFileHandle( _fh , 0 );
	_read=read;
	_w=minWindowSize;
	if(_w>_r)	_w=_r;
	_data = AllocPointer< byte >( _rs * _w );
	if( !_data ) fprintf(stderr,"Failed to allocate in FileBackedGrid::reset\n") , exit(0);
	if(_read)
	{
		_c = -_w;
		for(int i=0;i<_w;i++)	_readNext();
	}
	else	_c = 0;
}
void FileBackedGrid::_readNext(void)
{
	if( _c+_w<_r ) ReadFileHandle( _fh , _data + _rs*((_c+_w)%_w) , _rs ) , _rSize+=_rs;
	_c++;
}
void FileBackedGrid::_writeNext(void)
{
	_c++;
	if(_c-_w>=0)	WriteFileHandle( _fh , _data + _rs*((_c+_w)%_w) , _rs ) ,	_wSize+=_rs;
}
void FileBackedGrid::advance(void)
{
	if(_read)	_readNext();
	else		_writeNext();
}

//////////////////////////////
// CompressedFileBackedGrid //
//////////////////////////////
CompressedFileBackedGrid::CompressedFileBackedGrid(int rs,int r)
{
	_rSize=_wSize=0;
	_fh = CreateFileHandle( "." , NULL , true );
	_read=true;
	_r=r, _rs=rs;
	_c=_w=0;
	_sSize=int( (_rs+12)*1.002 );
	_data = NullPointer< byte >( );
	_scratch = AllocPointer< byte >( _sSize );
	if( !_scratch ) fprintf( stderr , "[ERROR] CompressedFileBackedGrid::CompressedFileBackedGrid failed\n") , exit(0);
}
CompressedFileBackedGrid::CompressedFileBackedGrid( int rs , int r , FileHandle fh )
{
	_rSize=_wSize=0;
	_fh = fh;
	_read=true;
	_r=r, _rs=rs;
	_c=_w=0;
	_sSize=int( (_rs+12)*1.002 );
	_data = NullPointer< byte >( );
	_scratch = AllocPointer< byte >( _sSize );
	if( !_scratch ) fprintf( stderr , "[ERROR] CompressedFileBackedGrid::CompressedFileBackedGrid failed\n") , exit(0);
}
CompressedFileBackedGrid::~CompressedFileBackedGrid(void)
{
	if(!_read)	for(int i=0;i<_w;i++)	_writeNext();
	FreePointer( _data );
	_r=_rs=_c=_w=0;
	CloseFileHandle( _fh );
	printf( "Read/Write Size: %lld MB / %lld MB\n" , _rSize>>20 , _wSize>>20 ) , fflush( stdout );
}
int	CompressedFileBackedGrid::rows(void) const { return _r;	}
int CompressedFileBackedGrid::rowSize(void) const { return _rs;	}
Pointer( byte ) CompressedFileBackedGrid::operator[](int idx)
{
	return _data + _rs*(idx%_w);
}
void CompressedFileBackedGrid::reset(bool read,int minWindowSize)
{
	if(!_read)	for(int i=0;i<_w;i++)	_writeNext();

	FreePointer( _data );

	SeekFileHandle( _fh , 0 );
	_read=read;
	_w=minWindowSize;
	if(_w>_r)	_w=_r;
	_data = AllocPointer< byte >( _rs * _w );
	if( !_data ) fprintf( stderr , "[ERROR] CompressedFileBackedGrid::reset failed\n" ) , exit(0);
	if( _read ) for( _c = -_w ; _c < 0 ; _readNext() )	;
	else _c = 0;
}
#ifdef _WIN32
#include <Util/ZLIB/ZLIB.h>
#else // !_WIN32
#include <zlib.h>
#endif // _WIN32
void CompressedFileBackedGrid::_readNext( void )
{
	if( _c+_w<_r )
	{
		uLong destLen=_rs,sourceLen;
		Bytef *dest = (Bytef*)(long long)( (_data)+_rs*((_c+_w)%_w) );
		const Bytef *source = (Bytef*)GetAddress( _scratch );
		if( ReadFileHandle( _fh , ( Pointer( byte ) )( GetPointer( sourceLen ) ) , sizeof(uLong) )!=1 ) fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_readNext failed\n" ) , exit(0);
		if( ReadFileHandle( _fh , _scratch , sourceLen )!=sourceLen ) fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_readNext failed\n") , exit(0);
		_rSize += sizeof(uLong)+sourceLen;

		switch(uncompress(dest,&destLen,source,sourceLen))
		{
			case Z_OK:         break;
			case Z_MEM_ERROR:  fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_readNext -- Z_MEM_ERROR\n" ) , exit(0);
			case Z_BUF_ERROR:  fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_readNext -- Z_BUF_ERROR\n" ) , exit(0);
			case Z_DATA_ERROR: fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_readNext -- Z_DATA_ERROR\n") , exit(0);
		}
	}
	_c++;
}
void CompressedFileBackedGrid::_writeNext(void)
{
	_c++;
	if( _c-_w>=0 && _c-_w<_r )
	{
		uLong destLen=_sSize,sourceLen=_rs;
		Bytef *dest = (Bytef*)GetAddress( _scratch );
		Bytef *source = (Bytef*)( (long long)(_data)+_rs*((_c+_w)%_w) );
		switch( compress( dest , &destLen , source , sourceLen ) )
		{
			case Z_OK:        break;
			case Z_MEM_ERROR: fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_writeNext -- Z_MEM_ERROR\n"   ) , exit(0);
			case Z_BUF_ERROR: fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_writeNext -- Z_BUF_ERROR\n"   ) , exit(0);
			default:          fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_writeNext -- unknown error\n" ) , exit(0);
		}
		if( WriteFileHandle( _fh , ( Pointer( byte ) )( GetPointer( destLen ) ) , sizeof(uLong) )!=1 ) fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_writeNext failed\n") , exit(0);
		if( WriteFileHandle( _fh , _scratch , destLen )!=destLen ) fprintf( stderr , "[ERROR] CompressedFileBackedGrid::_writeNext failed\n" ) , exit(0);
		_wSize += sizeof(uLong)+destLen;
	}
}
void CompressedFileBackedGrid::advance(void)
{
	if(_read)	_readNext();
	else		_writeNext();
}
