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
#ifndef JPEG_STREAM_INCLUDED
#define JPEG_STREAM_INCLUDED
#include "XPlatform.h"
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#ifdef _WIN32
#include "JPEG/jpeglib.h"
#include "JPEG/jerror.h"
#include "JPEG/jmorecfg.h"
#else // !_WIN32
#include <jpeglib.h>
#include <jerror.h>
#include <jmorecfg.h>
#endif // _WIN32

#include <Util/MultiStreamIO.h>
#include "ChannelConverter.h"
#include "Array.h"

struct misha_destination_mgr
{
	struct jpeg_destination_mgr pub;	/* public fields */
	VariableIOClientStream* outStream;	/* target stream */
	JOCTET * buffer;					/* start of buffer */
};
typedef misha_destination_mgr * misha_dest_ptr;

struct misha_source_mgr
{
	struct jpeg_source_mgr pub;			/* public fields */
	VariableIOClientStream* inStream; 	/* source stream */
	JOCTET * buffer;					/* start of buffer */
	bool start_of_file;					/* have we gotten any data yet? */
};
typedef misha_source_mgr * misha_src_ptr;

////////////////////////////////////////////////////////////////////
// Code adapted from example.c provided with libjpeg distribution //
////////////////////////////////////////////////////////////////////

struct my_error_mgr
{
	struct jpeg_error_mgr pub;    // "public" fields
	jmp_buf setjmp_buffer;        // for return to caller
};
typedef struct my_error_mgr * my_error_ptr;

struct JPEGReadInfo
{
	VariableIOClientStream* inStream;
	struct jpeg_decompress_struct cInfo;
	struct my_error_mgr	jErr;
	Pointer( unsigned char ) data;
};
struct JPEGWriteInfo
{
	VariableIOClientStream* outStream;
	struct jpeg_compress_struct cInfo;
	struct jpeg_error_mgr	jErr;
	Pointer( unsigned char ) data;
};

#ifndef SIZEOF
#define SIZEOF sizeof
#endif // !SIZEOF


////////////////////////////////////////////////////////////////////
// Code adapted from example.c provided with libjpeg distribution //
////////////////////////////////////////////////////////////////////

////////////////////
// Info for input //
////////////////////

#define INPUT_BUF_SIZE  4096

inline
METHODDEF(void)
init_source (j_decompress_ptr cinfo)
{
	misha_src_ptr src = (misha_src_ptr) cinfo->src;
	src->start_of_file = true;
}
inline
METHODDEF(boolean)
fill_input_buffer (j_decompress_ptr cinfo)
{
	misha_src_ptr src = (misha_src_ptr) cinfo->src;
	size_t nbytes;

	nbytes = src->inStream->read( GetPointer< byte >( src->buffer , INPUT_BUF_SIZE ) , INPUT_BUF_SIZE );

	if (nbytes <= 0)
	{
		if (src->start_of_file)	/* Treat empty input file as fatal error */
			ERREXIT(cinfo, JERR_INPUT_EMPTY);
		WARNMS(cinfo, JWRN_JPEG_EOF);
		/* Insert a fake EOI marker */
		src->buffer[0] = (JOCTET) 0xFF;
		src->buffer[1] = (JOCTET) JPEG_EOI;
		nbytes = 2;
	}

	src->pub.next_input_byte = src->buffer;
	src->pub.bytes_in_buffer = nbytes;
	src->start_of_file = false;
	return TRUE;
}



inline
METHODDEF(void)
skip_input_data (j_decompress_ptr cinfo, long num_bytes)
{
	misha_src_ptr src = (misha_src_ptr) cinfo->src;

	if (num_bytes > 0)
	{
		while (num_bytes > (long) src->pub.bytes_in_buffer)
		{
			num_bytes -= (long) src->pub.bytes_in_buffer;
			(void) fill_input_buffer(cinfo);
		}
		src->pub.next_input_byte += (size_t) num_bytes;
		src->pub.bytes_in_buffer -= (size_t) num_bytes;
	}
}

inline
METHODDEF(void)
term_source (j_decompress_ptr cinfo)
{
}

inline
GLOBAL(void)
jpeg_stdio_src (j_decompress_ptr cinfo, VariableIOClientStream * inStream )
{
	misha_src_ptr src;
	if (cinfo->src == NULL)
	{

		cinfo->src = (struct jpeg_source_mgr *)	(*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT, SIZEOF(misha_source_mgr));
		src = (misha_src_ptr) cinfo->src;
		src->buffer = (JOCTET *) (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,	INPUT_BUF_SIZE * SIZEOF(JOCTET));
	}

	src = (misha_src_ptr) cinfo->src;
	src->pub.init_source = init_source;
	src->pub.fill_input_buffer = fill_input_buffer;
	src->pub.skip_input_data = skip_input_data;
	src->pub.resync_to_restart = jpeg_resync_to_restart; /* use default method */
	src->pub.term_source = term_source;
	src->inStream = inStream;
	src->pub.bytes_in_buffer = 0; /* forces fill_input_buffer on first read */
	src->pub.next_input_byte = NULL; /* until buffer loaded */
}

/////////////////////
// Info for output //
/////////////////////
#define OUTPUT_BUF_SIZE  4096

inline
METHODDEF(void)
init_destination ( j_compress_ptr cinfo )
{
	misha_dest_ptr dest = (misha_dest_ptr) cinfo->dest;

	/* Allocate the output buffer --- it will be released when done with image */
	dest->buffer = (JOCTET *) (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE, OUTPUT_BUF_SIZE * SIZEOF(JOCTET));
	dest->pub.next_output_byte = dest->buffer;
	dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;
}


inline
METHODDEF(boolean)
empty_output_buffer (j_compress_ptr cinfo)
{
	misha_dest_ptr dest = (misha_dest_ptr) cinfo->dest;

	if( dest->outStream->write( GetPointer< byte >( dest->buffer , OUTPUT_BUF_SIZE ) , OUTPUT_BUF_SIZE ) != (size_t) OUTPUT_BUF_SIZE ) ERREXIT(cinfo, JERR_FILE_WRITE);
	dest->pub.next_output_byte = dest->buffer;
	dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;

	return TRUE;
}

inline
METHODDEF(void)
term_destination (j_compress_ptr cinfo)
{
	misha_dest_ptr dest = (misha_dest_ptr) cinfo->dest;
	size_t datacount = OUTPUT_BUF_SIZE - dest->pub.free_in_buffer;

	/* Write any data remaining in the buffer */
	if (datacount > 0) if( dest->outStream->write( GetPointer< byte >( dest->buffer , datacount ) , datacount ) != datacount ) ERREXIT(cinfo, JERR_FILE_WRITE);
	//  fflush(dest->outfile);
	/* Make sure we wrote the output file OK */
	// if (ferror(dest->outfile))
	//  ERREXIT(cinfo, JERR_FILE_WRITE);
}

inline
GLOBAL(void)
jpeg_stdio_dest (j_compress_ptr cinfo, VariableIOClientStream * outStream )
{
	misha_dest_ptr dest;
	if (cinfo->dest == NULL)/* first time for this JPEG object? */
		cinfo->dest = (struct jpeg_destination_mgr *) (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT, SIZEOF(misha_destination_mgr));

	dest = (misha_dest_ptr) cinfo->dest;
	dest->pub.init_destination = init_destination;
	dest->pub.empty_output_buffer = empty_output_buffer;
	dest->pub.term_destination = term_destination;
	dest->outStream = outStream;
}

inline
METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
	// cinfo->err really points to a my_error_mgr struct, so coerce pointer
	my_error_ptr myerr = (my_error_ptr) cinfo->err;

	// Always display the message.
	// We could postpone this until after returning, if we chose.
	(*cinfo->err->output_message) (cinfo);

	// Return control to the setjmp point
	longjmp(myerr->setjmp_buffer, 1);
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

inline void JPGGetImageInfo( char* fileName , int& width , int& height , int& channels , int& bytesPerChannel )
{
	FILE* fp;
	struct jpeg_decompress_struct cInfo;
	struct my_error_mgr	jErr;

	fp = fopen( fileName , "rb" );
	if(!fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);

	cInfo.err = jpeg_std_error( &jErr.pub );
	jErr.pub.error_exit = my_error_exit;
	if( setjmp( jErr.setjmp_buffer ) )
	{
		jpeg_destroy_decompress( &cInfo );
		fprintf( stderr , "[ERROR] JPEGGetImageSize: JPEG error occured\n" );
		return;
	}
	jpeg_create_decompress( &cInfo );
	jpeg_stdio_src( &cInfo, fp );	
	(void) jpeg_read_header( &cInfo , TRUE );
	width  = cInfo.image_width;
	height = cInfo.image_height;
	channels = cInfo.output_components;
	bytesPerChannel = 1;
	jpeg_destroy_decompress( &cInfo );
	fclose( fp );
}

template< int Channels , bool HDR >
void* JPGInitWrite( char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	if( HDR ) fprintf( stderr , "[WARNING] No HDR support for JPEG\n" );
	JPEGWriteInfo *info = (JPEGWriteInfo*)malloc( sizeof( JPEGWriteInfo ) );

	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , OUTPUT_BUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );

	info->cInfo.err = jpeg_std_error( &info->jErr );
	jpeg_create_compress(&info->cInfo);

	jpeg_stdio_dest( &info->cInfo , info->outStream );

	info->cInfo.image_width = width;    /* image width and height, in pixels */
	info->cInfo.image_height = height;
	info->cInfo.input_components = Channels;           /* # of color components per pixel */
	info->cInfo.in_color_space = Channels==3 ? JCS_RGB : JCS_GRAYSCALE;       /* colorspace of input image */

	jpeg_set_defaults( &info->cInfo );
	jpeg_set_quality( &info->cInfo , quality , TRUE);

	jpeg_start_compress( &info->cInfo , TRUE );
	info->data = AllocPointer< unsigned char >( width * info->cInfo.input_components );
	return info;
}

inline void JPGFinalizeWrite( void* v )
{
	JPEGWriteInfo* info = (JPEGWriteInfo*)v;
	jpeg_finish_compress( &info->cInfo );
	jpeg_destroy_compress( &info->cInfo );
	delete info->outStream;
	FreePointer( info->data );
	free( info );
}

inline void* JPGInitRead( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	JPEGReadInfo *info = (JPEGReadInfo*)malloc( sizeof( JPEGReadInfo ) );

	info->inStream = new VariableIOClientStream( );
	info->inStream->Initialize( fileName , OUTPUT_BUF_SIZE , 2 , true );
	info->inStream->SetServer( ioServer );

	info->cInfo.err = jpeg_std_error(&info->jErr.pub);
	info->jErr.pub.error_exit = my_error_exit;
	if (setjmp(info->jErr.setjmp_buffer))
	{
		jpeg_destroy_decompress(&info->cInfo);
		fprintf(stderr,"JPEG error occured\n");
		return 0;
	}

	jpeg_create_decompress(&info->cInfo);
	jpeg_stdio_src( &info->cInfo , info->inStream );

	(void) jpeg_read_header(&info->cInfo, TRUE);
	(void) jpeg_start_decompress(&info->cInfo);

	if( info->cInfo.output_components!=1 && info->cInfo.output_components!=3 )
	{
		fprintf( stderr , "[ERROR] JPEGInitRead: Only gray and RGB supported for jpeg: %d (%s)\n" , info->cInfo.output_components , fileName );
		return NULL;
	}
	width = info->cInfo.output_width;
	height = info->cInfo.output_height;
	info->data = AllocPointer< unsigned char >( width * info->cInfo.output_components );
	return info;
}

inline void JPGFinalizeRead( void* v )
{
	JPEGReadInfo* info = ( JPEGReadInfo* )v;
	// The JPEG reader will not allow us to close until we've actually read the scan-lines.
	// How annoying!!!
	if( info->cInfo.output_scanline < info->cInfo.output_height )
	{
		unsigned char* pixels = (unsigned char*) malloc( sizeof ( unsigned char) * 3 * info->cInfo.image_width );
		JSAMPROW row_pointers[1];
		row_pointers[0] = pixels;
		while( info->cInfo.output_scanline < info->cInfo.output_height ) jpeg_read_scanlines( &info->cInfo , row_pointers , 1 );
	}
	(void) jpeg_finish_decompress( &info->cInfo );
	jpeg_destroy_decompress( &info->cInfo );
	delete info->inStream;
	FreePointer( info->data );
	free( info );
}

template< int Channels , class ChannelType >
void JPGReadRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	JPEGReadInfo* info = (JPEGReadInfo*)v;
	if(info->cInfo.output_scanline >= info->cInfo.output_height)
		fprintf( stderr , "[ERROR] JPEGReadRow: Trying to read beyond the end of the jpeg file: %d >= %d\n" , info->cInfo.output_scanline , info->cInfo.output_height ) , exit( 0 );
	JSAMPROW row_pointers[1];

	row_pointers[0] = GetAddress( info->data );
	jpeg_read_scanlines( &info->cInfo , row_pointers , 1 );
	ConvertRow< unsigned char , ChannelType >( ( ConstPointer( unsigned char ) )info->data , pixels , info->cInfo.image_width , info->cInfo.output_components , Channels );
}
template< int Channels , class ChannelType >
void JPGWriteRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	JPEGWriteInfo* info = (JPEGWriteInfo*)v;
	JSAMPROW row_pointer[1];

	row_pointer[0] = GetAddress( info->data );
	ConvertRow< ChannelType , unsigned char >( pixels , ( Pointer( unsigned char ) )info->data , info->cInfo.image_width , Channels , info->cInfo.input_components );
	(void) jpeg_write_scanlines( &info->cInfo , row_pointer , 1 );
}
#endif // JPEG_STREAM_INCLUDED