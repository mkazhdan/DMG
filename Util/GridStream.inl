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
#include "Time.h"

///////////////////////////
// MultiSocketBackedGrid //
///////////////////////////
template< class Real , int Channels >
MultiSocketBackedGrid< Real , Channels >::MultiSocketBackedGrid( Socket* socks , int* rowSizes , int sockCount , int rows )
{
	_read = true;
	_sockCount = sockCount;
	_rows = rows;
	_rowSize = 0;
	for( int i=0 ; i<_sockCount ; i++ ) _rowSize += rowSizes[i];

	_rowSizes = ( int* )malloc( sizeof(int) * _sockCount );
	_socks = ( Socket* )malloc( sizeof(Socket) * _sockCount );
	_subData = AllocPointer< Pointer( Real ) >( _sockCount );
	if( !_rowSizes || !_socks || !_subData ) fprintf( stderr , "Failed to allocate\n" ) ,  exit(0);
	for( int i=0 ; i<_sockCount ; i++ )
	{
		_rowSizes[i] = rowSizes[i];
		_socks[i] = socks[i];
		_subData[i] = AllocPointer< Real >( Channels * rowSizes[i] );
		if( !_subData[i] ) fprintf( stderr , "Failed to allocate\n" ) , exit(0);
	}

	_data = AllocPointer< Real >( _rowSize * Channels );
	if( !_data ) fprintf( stderr , "Failed to allocate\n" ) , exit(0);
}
template< class Real , int Channels >
MultiSocketBackedGrid< Real , Channels >::~MultiSocketBackedGrid(void)
{
	for( int i=0 ; i<_sockCount ; i++ ) FreePointer( _subData[i] );
	FreePointer( _data );
	FreePointer( _subData );
	free( _socks );
	free( _rowSizes );

	_rows = _rowSize = 0;
	_rowSizes = NULL;
	_socks = NULL;
}
template< class Real , int Channels > int MultiSocketBackedGrid< Real , Channels >::rows( void ) const { return _rows; }
template< class Real , int Channels > int MultiSocketBackedGrid< Real , Channels >::rowSize( void ) const { return _rowSize * Channels * sizeof(Real); }
template< class Real , int Channels >
Pointer( byte ) MultiSocketBackedGrid< Real , Channels >::operator[]( int idx )
{
	if( _read && !_readComplete )
	{
		int offset = 0;
		for( int i=0 ; i<_sockCount ; i++ )
		{
			ReceiveOnSocket( _socks[i] , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) , "MultiSocketBackedGrid::[]" );
			memcpy( _data+offset*Channels , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) );
			offset += _rowSizes[i];
		}
		_readComplete = true;
	}
	return ( Pointer( byte ) )_data;
}
template< class Real , int Channels >
void MultiSocketBackedGrid< Real , Channels >::advance( void )
{
	if( !_read )
	{
		int offset = 0;
		for( int i=0 ; i<_sockCount ; i++ )
		{
			memcpy( _subData[i] , _data+offset*Channels , _rowSizes[i]*Channels*sizeof(Real) );
			SendOnSocket( _socks[i] , ( ConstPointer( Real ) )_subData[i] , _rowSizes[i]*Channels*sizeof(Real) , "MultiSocketBackedGrid::advance" );
			offset += _rowSizes[i];
		}
	}
	else
	{
		if( !_readComplete )
		{
			int offset = 0;
			for( int i=0 ; i<_sockCount ; i++ )
			{
				memcpy( _subData[i] , _data+offset*Channels , _rowSizes[i]*Channels*sizeof(Real) );
				ReceiveOnSocket( _socks[i] , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) , "MultiSocketBackedGrid::advance"  );
				offset += _rowSizes[i];
			}
		}
		_readComplete = false;
	}
	_current++;
}
template< class Real , int Channels >
void MultiSocketBackedGrid< Real , Channels >::reset( bool read , int minWindowSize )
{
	_read = read;
	_readComplete = !_read;
	_current = 0;
}

