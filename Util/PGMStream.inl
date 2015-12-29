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
#ifndef PGM_STREAM_INCLUDED
#define PGM_STREAM_INCLUDED
#include "XPlatform.h"
#include <stdio.h>
#include <Util/BaseMultiStreamIO.h>
#include "ChannelConverter.h"


struct PGMReadInfo
{
	FILE* fp;
	int width;
	int ascii , pgm;
	Pointer( unsigned char ) data;
};

inline void PGMGetImageInfo( char* fileName , int& width , int& height , int& channels , int& bytesPerChannel )
{
	FILE* fp = fopen( fileName , "rb" );
	if( !fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	char header[512];
	int maxValue;
	if( fscanf( fp , "%s %d %d %d" , header , &width , &height , &maxValue )!=4 ) 
		fprintf( stderr , "Failed to parse header\n" ) , exit( 0 );
	while (fgetc( fp ) != '\n');
	if( ( strcasecmp( header , "P5" ) ) && ( strcasecmp( header , "P2" ) ) &&  // binary and ascii pgm
		( strcasecmp( header , "P6" ) ) && ( strcasecmp( header , "P3" ) ) )   // binary and ascii ppm
		fprintf( stderr , "Error: binary/ascii PPM (P6/P3) or binary/ascii PGM only (P5/P2): %s\n" , header ) , exit( 0 );
	int ascii = ( !strcasecmp( header , "P2" ) ) || ( !strcasecmp( header , "P3" ) );
	int pgm   = ( !strcasecmp( header , "P5" ) ) || ( !strcasecmp( header , "P2" ) );

	channels = pgm ? 1 : 3;
	bytesPerChannel = 1;
	
	fclose( fp );
}

template< int Channels , bool HDR >
void* PGMInitWrite( char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	fprintf( stderr , "[ERROR] PGMInitWrite: Function not supported\n" ) , exit( 0 );
}
inline void PGMFinalizeWrite( void* v )
{
	fprintf( stderr , "[ERROR] PGMFinalizeWrite: Function not supported\n" ) , exit( 0 );
}
template< int Channels , class ChannelType >
void PGMWriteRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	fprintf( stderr , "[ERROR] PGMWriteRow: Function not supported\n" ) , exit( 0 );
}

inline void* PGMInitRead( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	PGMReadInfo* info = (PGMReadInfo*)malloc( sizeof( PGMReadInfo ) );
	info->fp = fopen( fileName , "rb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	char header[512];
	int maxValue;
	if( fscanf( info->fp , "%s %d %d %d" , header , &width , &height , &maxValue )!=4 ) 
		fprintf( stderr , "Failed to parse header\n" ) , exit( 0 );
	while (fgetc( info->fp ) != '\n');
	if( ( strcasecmp( header , "P5" ) ) && ( strcasecmp( header , "P2" ) ) &&  // binary and ascii pgm
		( strcasecmp( header , "P6" ) ) && ( strcasecmp( header , "P3" ) ) )   // binary and ascii ppm
		fprintf( stderr , "Error: binary/ascii PPM (P6/P3) or binary/ascii PGM only (P5/P2): %s\n" , header ) , exit( 0 );
	info->ascii = ( !strcasecmp( header , "P2" ) ) || ( !strcasecmp( header , "P3" ) );
	info->pgm   = ( !strcasecmp( header , "P5" ) ) || ( !strcasecmp( header , "P2" ) );

	info->width = width;
	info->data = AllocPointer< unsigned char >( width * 3 );
	if( !info->data ) fprintf( stderr , "Could not allocate memory for bitmap data\n" ) , exit( 0 );

	return info;
}

template< int Channels , class ChannelType >
void PGMReadRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	PGMReadInfo* info = ( PGMReadInfo* )v;
	unsigned char pixel_byte;
	if( info->pgm )
		if( info->ascii )
			for( int i=0 ; i<info->width ; i++ )
			{
				int pixel_int;
				fscanf( info->fp , " %d " , &pixel_int );
				pixel_byte = (unsigned char) pixel_int;
				info->data[3*i] = info->data[3*i+1] = info->data[3*i+2] = pixel_byte;
			}
		else
			for( int i=0 ; i<info->width ; i++ )
			{
				fscanf( info->fp , "%c" , &pixel_byte );
				info->data[3*i] = info->data[3*i+1] = info->data[3*i+2] = pixel_byte;
			}
	else
		if( info->ascii )
		{
			int pixel_int[3];
			for( int i=0 ; i<info->width ; i++ )
			{
				fscanf( info->fp , " %d %d %d " , pixel_int+0 , pixel_int+1 , pixel_int+2 );
				info->data[3*i+0] = (unsigned char) pixel_int[0];
				info->data[3*i+1] = (unsigned char) pixel_int[1];
				info->data[3*i+2] = (unsigned char) pixel_int[2];
			}
		}
		else
			fread( info->data , sizeof(unsigned char) , 3*info->width , info->fp );

	ConvertRow< unsigned char , ChannelType >( ( ConstPointer(unsigned char) )info->data , pixels , info->width , info->pgm ? 1 : 3 , Channels );
}

inline void PGMFinalizeRead( void* v  )
{
	PGMReadInfo* info = ( PGMReadInfo* )v;
	fclose( info->fp );
	FreePointer( info->data );
	free( info );
}
#endif // PGM_STREAM_INCLUDED