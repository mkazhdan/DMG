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
#ifndef KRO_STREAM_INCLUDED
#define KRO_STREAM_INCLUDED
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <Util/BaseMultiStreamIO.h>
#include "ChannelConverter.h"

// Code courtesy of: http://www.codeguru.com/forum/showthread.php?t=292902
inline void endian_swap(unsigned short& x)
{
	x = (x>>8) | 
		(x<<8);
}

inline void endian_swap(unsigned int& x)
{
	x = (x>>24) | 
		((x<<8) & 0x00FF0000) |
		((x>>8) & 0x0000FF00) |
		(x<<24);
}

// __int64 for MSVC, "long long" for gcc
inline void endian_swap(uint64_t& x)
{
	x = (x>>56) | 
		((x<<40) & 0x00FF000000000000) |
		((x<<24) & 0x0000FF0000000000) |
		((x<<8)  & 0x000000FF00000000) |
		((x>>8)  & 0x00000000FF000000) |
		((x>>24) & 0x0000000000FF0000) |
		((x>>40) & 0x000000000000FF00) |
		(x<<56);
}

/*
Header is 20 bytes long :
3 bytes : "KRO" signature in hex 0x4B 0x52 0x4F
1 byte  : 0x01 version
unsigned long : Width
unsigned long : Height
unsigned long : depth = > 8 bits, 16 bits, 32 bits
unsigned long : ncomp => number of compoment, 4 by default, RGB + Alpha
*/

const char KRO_HEADER[] = { 0x4B , 0x52 , 0x4F };
const char KRO_VERSION  = 0x01;
struct KROHeader
{
	char header[3];
	char version;
	unsigned int width;
	unsigned int height;
	unsigned int depth;
	unsigned int ncomp;
};
struct KROInfo
{
	FILE* fp;
	int width , height;
	int bytesPerChannel , channelsPerPixel;
	Pointer( byte ) data;
};

inline void KROGetImageInfo( char* fn , int& width , int& height , int& channels , int& bytesPerChannel )
{
	KROInfo info;
	info.fp = fopen( fn , "rb" );
	if( !info.fp ) fprintf( stderr , "Failed to open: %s\n" , fn ) , exit(0);

	KROHeader header;
	fread( &header , 1 , sizeof( KROHeader ) , info.fp );
	endian_swap( header.width ) , endian_swap( header.height ) , endian_swap( header.depth ) , endian_swap( header.ncomp );
	if( header.header[0]!=KRO_HEADER[0] || header.header[1]!=KRO_HEADER[1] || header.header[2]!=KRO_HEADER[2] )
		fprintf( stderr , "Invalid header in %s\n" , fn ) , exit( 0 );

	width           = header.width;
	height          = header.height;
	channels        = header.ncomp;
	bytesPerChannel = header.depth / 8;
	fclose( info.fp );
}
inline void* KROInitRead( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	KROInfo *info = ( KROInfo* )malloc( sizeof( KROInfo ) );
	info->fp = fopen( fileName , "rb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	KROHeader header;
	fread( &header , 1 , sizeof( KROHeader ) , info->fp );
	endian_swap( header.width ) , endian_swap( header.height ) , endian_swap( header.depth ) , endian_swap( header.ncomp );
	if( header.header[0]!=KRO_HEADER[0] || header.header[1]!=KRO_HEADER[1] || header.header[2]!=KRO_HEADER[2] )
		fprintf( stderr , "Invalid header in %s\n" , fileName ) , exit( 0 );
	if( header.depth!=8 && header.depth!=16 && header.depth!=32 )
		fprintf( stderr , "Only 8 , 16, and 32 bits-per-channel supported in KRO file-format: %d\n" , header.ncomp ) , exit( 0 );
	if( header.ncomp!=1 && header.ncomp!=3 && header.ncomp!=4 ) fprintf( stderr , "[ERROR] KROInitRead: Only gray, RGB, and RGBA supported for kro: %d\n" , header.ncomp ) , exit( 0 );
	info->width  = header.width;
	info->height = header.height;
	info->bytesPerChannel  = header.depth / 8;
	info->channelsPerPixel = header.ncomp;

	info->data = AllocPointer< byte >( info->width * info->bytesPerChannel * info->channelsPerPixel );
	if( !info->data ) fprintf( stderr , "[ERROR] KROInitRead: Failed to allocate KRO row of size: %d\n" , header.width ) , exit( 0 );

	width  = info->width;
	height = info->height;
	return info;
}
template< int Channels , bool HDR >
void* KROInitWrite( char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	KROInfo *info = ( KROInfo* )malloc( sizeof( KROInfo ) );
	info->fp = fopen( fileName , "wb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	KROHeader header;
	header.header[0] = KRO_HEADER[0] , header.header[1] = KRO_HEADER[1] , header.header[2] = KRO_HEADER[2];
	header.version = KRO_VERSION;
	header.width   = width;
	header.height  = height;
	header.depth   = HDR ? 16 : 8;
	header.ncomp   = Channels;

	info->bytesPerChannel = header.depth / 8;
	info->data = AllocPointer< byte >( info->width * info->bytesPerChannel * info->channelsPerPixel );
	if( !info->data ) fprintf( stderr , "Failed to allocate KRO row of size: %d\n" , header.width ) , exit( 0 );

	endian_swap( header.width ) , endian_swap( header.height ) , endian_swap( header.depth ) , endian_swap( header.ncomp );
	fwrite( &header , 1 , sizeof( KROHeader ) , info->fp );

	info->width  = header.width;
	info->height = header.height;
	info->bytesPerChannel  = header.depth / 8;
	info->channelsPerPixel = header.ncomp;

	return info;
}
template< int Channels , class ChannelType >
void KROReadRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	KROInfo* info = ( KROInfo* )v;
	fread( info->data , 1 , info->bytesPerChannel * info->channelsPerPixel * info->width , info->fp );
	switch( info->bytesPerChannel )
	{
	case 1: ConvertRow< unsigned char , ChannelType >( ( ConstPointer(unsigned char) )info->data , pixels , info->width , info->channelsPerPixel , Channels ) ; break;
	case 2: ConvertRow< uint16_t      , ChannelType >( ( ConstPointer(uint16_t     ) )info->data , pixels , info->width , info->channelsPerPixel , Channels ) ; break;
	case 4: ConvertRow< float         , ChannelType >( ( ConstPointer(float        ) )info->data , pixels , info->width , info->channelsPerPixel , Channels ) ; break;
	default: fprintf( stderr , "[ERROR] KROReadRow: Invalid number of bytes per channel: %d\n" , info->bytesPerChannel ) , exit( 0 );
	}
}
template< int Channels , class ChannelType >
void KROWriteRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	KROInfo* info = ( KROInfo* )v;
	if( Channels!=info->channelsPerPixel ) fprintf( stderr , "[ERROR] KROWriteRow: channels don't match: %d != %d\n" , Channels , info->channelsPerPixel ) , exit( 0 );
	switch( info->bytesPerChannel )
	{
	case 1: ConvertRow< ChannelType , unsigned char >( pixels , ( Pointer(unsigned char) )info->data , info->width , Channels , info->channelsPerPixel ) ; break;
	case 2: ConvertRow< ChannelType , uint16_t      >( pixels , ( Pointer(uint16_t     ) )info->data , info->width , Channels , info->channelsPerPixel ) ; break;
	case 4: ConvertRow< ChannelType , float         >( pixels , ( Pointer(float        ) )info->data , info->width , Channels , info->channelsPerPixel ) ; break;
	}
	fwrite( info->data , 1 , info->bytesPerChannel * info->channelsPerPixel * info->width , info->fp );
}

inline void KROFinalizeRead( void* v )
{
	KROInfo* info = (KROInfo*)v;
	fclose( info->fp );
	FreePointer( info->data );
	free( info );
}
inline void KROFinalizeWrite( void* v )
{
	KROInfo* info = (KROInfo*)v;
	fclose( info->fp );
	FreePointer( info->data );
	free( info );
}
#endif // KRO_STREAM_INCLUDED