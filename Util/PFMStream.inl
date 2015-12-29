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
#ifndef PFM_STREAM_INCLUDED
#define PFM_STREAM_INCLUDED
#include <stdio.h>
#include <stdlib.h>
#include <Util/BaseMultiStreamIO.h>
#include "ChannelConverter.h"

struct PFMReadInfo
{
	FILE* fp;
	int width , height;
	Pointer( byte ) data;
};

inline void PFMGetImageInfo( char* fn , int& width , int& height , int& channels , int& bytesPerChannel )
{
	PFMReadInfo info;
	info.fp = fopen( fn , "rb" );
	if( !info.fp ) fprintf( stderr , "Failed to open: %s\n" , fn ) , exit(0);

	char buf[1024];
	fgets( buf , 1024 , info.fp );
	fgets( buf , 1024 , info.fp );
	while (buf[0] == '#') fgets( buf , 1024 , info.fp );
	sscanf(buf, "%d %d" , &info.width , &info.height);
	if( info.width<0 || info.height<0 )
	{
		fprintf(stderr, "Invalid header or image not valid PFM.\n");
		exit( 0 );
	}
	width  = info.width;
	height = info.height;
	channels = 3;
	bytesPerChannel = sizeof( float );

	fclose( info.fp );
}
template< int Channels , bool HDR >
void* PFMInitWrite( char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	fprintf( stderr , "[ERROR] PFMInitWrite: Function not supported\n" ) , exit( 0 );
}
inline void PFMFinalizeWrite( void* v )
{
	fprintf( stderr , "[ERROR] PFMFinalizeWrite: Function not supported\n" ) , exit( 0 );
}
template< int Channels , class ChannelType >
void PFMWriteRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	fprintf( stderr , "[ERROR] PFMWriteRow: Function not supported\n" ) , exit( 0 );
}

inline void* PFMInitRead( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	PFMReadInfo *info=(PFMReadInfo*)malloc(sizeof(PFMReadInfo));
	info->fp = fopen( fileName , "rb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	char buf[1024];
	fgets( buf , 1024 , info->fp );
	fgets( buf , 1024 , info->fp );
	while (buf[0] == '#') fgets( buf , 1024 , info->fp );
	sscanf(buf, "%d %d", &info->width, &info->height);
	if( info->width<0 || info->height<0 )
	{
		fprintf(stderr, "Invalid header or image not valid PFM.\n");
		exit( 0 );
	}
	fgets( buf, 1024 , info->fp );
	while ( buf[0] == '#') fgets( buf , 1024 , info->fp );


	width  = info->width;
	height = info->height;
	info->data = AllocPointer< byte >( sizeof(float)*3*width );

	return info;
}
template< int Channels , class ChannelType >
void PFMReadRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	PFMReadInfo* info = ( PFMReadInfo* )v;
	if( !fread( info->data , sizeof(float) , info->width * 3 , info->fp ) ) fprintf( stderr , "Failed to read PFM row of size %d\n" , info->width ) , exit( 0 );
	ConvertRow< float , ChannelType >( ( ConstPointer( float ) )info->data , pixels , info->width , 3 , Channels );
}

inline void PFMFinalizeRead( void* v )
{
	PFMReadInfo* info = (PFMReadInfo*)v;
	fclose( info->fp );
	FreePointer( info->data );
	free( info );
}
#endif // PFM_STREAM_INCLUDED