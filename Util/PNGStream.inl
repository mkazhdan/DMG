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
#ifndef PNG_STREAM_INCLUDED
#define PNG_STREAM_INCLUDED

#include <stdlib.h>
#include <stdint.h>
#define PNG_SETJMP_NOT_SUPPORTED // Trying to get around the double inclusion of setjmp.h
#ifdef _WIN32
#include <Util/PNG/png.h>
#include <Util/ZLIB/zlib.h>
#else // !_WIN32
#include <png.h>
#include <zlib.h>
#endif // _WIN32
#include <Util/MultiStreamIO.h>


#define MONOCHROME_FIX 1

struct PNGWriteInfo
{
	VariableIOClientStream* outStream;
	png_structp png_ptr;
	png_infop info_ptr;
	Pointer( byte ) data;
};
struct PNGReadInfo
{
	VariableIOClientStream* inStream;
	png_structp png_ptr;
	png_infop info_ptr , end_info;
	int width;
	Pointer( byte ) data;
};

inline void myPngWriteFunction( png_structp png_ptr , png_bytep buffer , png_size_t byteNum )
{
	VariableIOClientStream* outStream = (VariableIOClientStream*)png_get_io_ptr( png_ptr );
	if( outStream->write( GetPointer< byte >( buffer , byteNum ) , byteNum ) != byteNum ) fprintf( stderr , "PNG Failed to write %zd bytes\n" , byteNum );
}
inline void myPngReadFunction( png_structp png_ptr , png_bytep buffer , png_size_t byteNum )
{
	VariableIOClientStream* inStream = (VariableIOClientStream*)png_get_io_ptr( png_ptr );
	if( inStream->read( GetPointer< byte >( buffer , byteNum ) , byteNum ) != byteNum ) fprintf( stderr , "PNG Failed to read %zd bytes\n" , byteNum );
}
inline void myPngFlushFunction( png_structp png_ptr )
{
}


inline void PNGGetImageInfo( char* fileName , int& width , int& height , int& channels , int& bytesPerChannel )
{
	FILE* fp;
	png_structp png_ptr;
	png_infop info_ptr , end_info;

	fp = fopen( fileName , "rb" );
	if( !fp ) fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);

	png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING , 0 , 0 , 0 );
	if( !png_ptr ) fprintf(stderr,"Failed to create PNG read structure\n") , exit( 0 );
	info_ptr = png_create_info_struct( png_ptr );
	if( !info_ptr ) fprintf(stderr,"Failed to create PNG info structure 1\n") , exit( 0 );
	end_info = png_create_info_struct(png_ptr);
	if( !end_info ) fprintf(stderr,"Failed to create PNG info structure 2\n") , exit(0);
	png_init_io( png_ptr, fp );
	png_read_info( png_ptr, info_ptr );
#if MONOCHROME_FIX
	png_set_expand( png_ptr );
	png_read_update_info( png_ptr , info_ptr );
#endif // MONOCHROME_FIX
	width           = png_get_image_width ( png_ptr , info_ptr );
	height          = png_get_image_height( png_ptr , info_ptr );
	channels        = png_get_channels    ( png_ptr , info_ptr );
#if MONOCHROME_FIX
	bytesPerChannel = std::max< int >( 1 , png_get_bit_depth( png_ptr , info_ptr ) / 8 );
#else // !MONOCHROME_FIX
	bytesPerChannel = png_get_bit_depth( png_ptr , info_ptr ) / 8;
#endif // MONOCHROME_FIX

	png_destroy_read_struct( &png_ptr , &info_ptr , &end_info );
	fclose( fp );
}

template< int Channels , bool HDR >
void* PNGInitWrite( char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	PNGWriteInfo* info = (PNGWriteInfo*)malloc(sizeof(PNGWriteInfo));
	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );

	info->png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING , 0 , 0 , 0 );
	if(!info->png_ptr)	fprintf(stderr,"Failed to create png write struct\n")	,	exit(0);
	info->info_ptr = png_create_info_struct(info->png_ptr);
	if(!info->info_ptr)	fprintf(stderr,"Failed to create png info struct\n")	,	exit(0);

	png_set_write_fn( info->png_ptr , info->outStream , myPngWriteFunction , myPngFlushFunction );

	png_set_compression_level(info->png_ptr,Z_BEST_SPEED);

	int pngColorType;
	switch( Channels )
	{
	case 1: pngColorType = PNG_COLOR_TYPE_GRAY ; break;
	case 3: pngColorType = PNG_COLOR_TYPE_RGB  ; break;
	case 4: pngColorType = PNG_COLOR_TYPE_RGBA ; break;
	};
	png_set_IHDR(info->png_ptr, info->info_ptr, width, height,
		HDR ? 16 : 8 , pngColorType ,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(info->png_ptr, info->info_ptr);

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if(swap)	png_set_swap(info->png_ptr);
	}
	info->data = AllocPointer< byte >( ( HDR ? sizeof( uint16_t ) : sizeof( unsigned char ) ) * Channels * width );
	return info;
}
template< int Channels , class ChannelType >
void PNGWriteRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	PNGWriteInfo* info = (PNGWriteInfo*)v;
	int width           = png_get_image_width( info->png_ptr , info->info_ptr );
	int channels        = png_get_channels   ( info->png_ptr , info->info_ptr );
	int bytesPerChannel = png_get_bit_depth  ( info->png_ptr , info->info_ptr ) / 8;
	if( bytesPerChannel==1 ) ConvertRow< ChannelType , unsigned char >( pixels , ( Pointer(unsigned char) )info->data , width , Channels , channels );
	else                     ConvertRow< ChannelType , uint16_t      >( pixels , ( Pointer(uint16_t     ) )info->data , width , Channels , channels );
	png_bytep row_pointer = png_bytep( info->data );
	png_write_row( info->png_ptr , row_pointer );
}
inline void PNGFinalizeWrite( void* v )
{
	PNGWriteInfo* info = (PNGWriteInfo*)v;
	png_write_end(info->png_ptr, NULL);
	png_destroy_write_struct(&info->png_ptr, &info->info_ptr);
	delete info->outStream;
	FreePointer( info->data );
	free( info );
}

inline void* PNGInitRead( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	PNGReadInfo* info = (PNGReadInfo*)malloc( sizeof( PNGReadInfo ) );
	info->inStream = new VariableIOClientStream( );
	info->inStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , true );
	info->inStream->SetServer( ioServer );

	info->png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING , 0 , 0 , 0 );
	if( !info->png_ptr ) fprintf( stderr , "Failed to create PNG read structure\n" ) , exit(0);

	info->info_ptr = png_create_info_struct( info->png_ptr );
	if( !info->info_ptr ) fprintf( stderr , "Failed to create PNG info structure 1\n" ) , exit(0);

	info->end_info = png_create_info_struct(info->png_ptr);
	if( !info->end_info ) fprintf( stderr , "Failed to create PNG info structure 2\n" ) , exit(0);

	png_set_read_fn( info->png_ptr , info->inStream , myPngReadFunction );

	png_read_info( info->png_ptr, info->info_ptr );

#if MONOCHROME_FIX
	png_set_expand( info->png_ptr );
	png_read_update_info( info->png_ptr , info->info_ptr );
	png_uint_32 _width , _height;
	int bit_depth , color_type , ncomp;
	png_get_IHDR( info->png_ptr , info->info_ptr , &_width , &_height , &bit_depth , &color_type , NULL , NULL , NULL );
	ncomp = png_get_channels( info->png_ptr , info->info_ptr );
	info->width = width = (int)_width , height = (int)_height;
	if( ncomp<1 || ncomp>4 )            fprintf( stderr , "[ERROR] Invalid number of components: %d\n" , ncomp ) , exit( 0 );
	if( bit_depth!=8 && bit_depth!=16 ) fprintf( stderr , "[ERROR] Invalid bit depth: %d\n" , bit_depth ) , exit(0);
	info->data = AllocPointer< byte >( (bit_depth/8) * width * ncomp );
#else // !MONOCHROME_FIX
	info->width = width = png_get_image_width( info->png_ptr , info->info_ptr );
	height = png_get_image_height( info->png_ptr , info->info_ptr );
	int ncomp = png_get_channels( info->png_ptr , info->info_ptr );
	int bit_depth = png_get_bit_depth( info->png_ptr , info->info_ptr );
	int color_type = png_get_color_type( info->png_ptr , info->info_ptr );
	if( width<=0 || height<=0 ) fprintf( stderr , "[ERROR] Invalid resolution [ %d x %d ]\n" , width , height ) , exit(0);
	if( ncomp<1 || ncomp>4 ) fprintf( stderr , "[ERRROR] Invalid number of components: %d\n" , ncomp ) , exit( 0 );
	if( bit_depth!=8 && bit_depth!=16 ) fprintf( stderr , "[ERROR] Invalid bit depth: %d\n" , bit_depth ) , exit(0);
	if( color_type==PNG_COLOR_TYPE_PALETTE ) png_set_expand( info->png_ptr ) , printf( "Expanding PNG color pallette\n" );
	info->data = AllocPointer< byte >( (bit_depth/8) * width * ncomp );
#endif // MONOCHROME_FIX

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( info->png_ptr );
	}
	return info;
}
template< int Channels , class ChannelType >
void PNGReadRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	PNGReadInfo* info = (PNGReadInfo*)v;
	int width     = png_get_image_width( info->png_ptr , info->info_ptr );
	int ncomp     = png_get_channels   ( info->png_ptr , info->info_ptr );
	int bit_depth = png_get_bit_depth  ( info->png_ptr , info->info_ptr );

	png_bytep row = (png_bytep)GetAddress( info->data );
	png_read_row( info->png_ptr , row , NULL );

	if( bit_depth==8 ) ConvertRow< unsigned char , ChannelType >( ( ConstPointer(unsigned char) )info->data , pixels , width , ncomp , Channels );
	else               ConvertRow< uint16_t      , ChannelType >( ( ConstPointer(uint16_t     ) )info->data , pixels , width , ncomp , Channels );
}

inline void PNGFinalizeRead( void* v )
{
	PNGReadInfo* info = (PNGReadInfo*)v;
	png_destroy_read_struct(&info->png_ptr, &info->info_ptr, &info->end_info);
	delete info->inStream;
	FreePointer( info->data );
	free( info );
}
#endif // PNG_STREAM_INCLUDED