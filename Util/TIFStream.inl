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
#include <stdlib.h>
#include <stdint.h>
#ifdef _WIN32
#include <Util/TIFF/tiffio.h>
#else // !_WIN32
#include <tiffio.h>
#endif // _WIN32
#include <Util/BaseMultiStreamIO.h>
#include "ChannelConverter.h"

// Code taken from:
// http://www.codeguru.com/cpp/g-m/bitmap/otherformats/article.php/c4933/

struct TIFFInfo
{
	TIFF* tiff;
	int width;
	uint16 bitsPerSample , samplesPerPixel;
	Pointer( byte ) data;
};

inline void TIFGetImageInfo( char* fn , int& width , int& height , int& channels , int& bytesPerChannel )
{
	TIFF* tiff;
	tiff = TIFFOpen( fn , "r" );
	if( !tiff ) fprintf( stderr , "[ERROR] TIFGetImageInfo: TIFFOpen failed\n" ) , exit( 0 );

	TIFFGetField( tiff , TIFFTAG_IMAGEWIDTH      , &width  );
	TIFFGetField( tiff , TIFFTAG_IMAGELENGTH     , &height );
	TIFFGetField( tiff , TIFFTAG_SAMPLESPERPIXEL , &channels );
	TIFFGetField( tiff , TIFFTAG_BITSPERSAMPLE   , &bytesPerChannel ) , bytesPerChannel /= 8;
	TIFFClose( tiff );
}
inline void* TIFInitRead( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	uint16 temp16;
	TIFFInfo* info = (TIFFInfo*)malloc(sizeof(TIFFInfo));
	info->tiff = TIFFOpen( fileName , "r" );
	if( !info->tiff ) fprintf( stderr , "TIFFOpen failed\n" ) , exit( 0 );

	TIFFGetField( info->tiff , TIFFTAG_IMAGEWIDTH  , &width  );
	info->width = width;
	TIFFGetField( info->tiff , TIFFTAG_IMAGELENGTH , &height );

	TIFFGetField( info->tiff , TIFFTAG_BITSPERSAMPLE   , &(info->bitsPerSample) );
	if( info->bitsPerSample!=8 && info->bitsPerSample!=16 ) fprintf( stderr , "Unsupported number of bits per sample: %d\n" , info->bitsPerSample ) , exit( 0 );
	TIFFGetField( info->tiff , TIFFTAG_SAMPLESPERPIXEL , &info->samplesPerPixel );
	if( info->samplesPerPixel!=1 && info->samplesPerPixel!=3 && info->samplesPerPixel!=4 ) fprintf( stderr , "Unsupported number of samples per pixel: %d\n" , info->bitsPerSample ) , exit( 0 );
	int scanLineSize1 = TIFFScanlineSize( info->tiff ) , scanLineSize2 = ( info->bitsPerSample / 8 ) * info->samplesPerPixel * info->width;
	if( scanLineSize1!=scanLineSize2 ) fprintf( stderr , "Scanline sizes do not agree: %d!=%d\n" , scanLineSize1 , scanLineSize2 ) , exit( 0 );
	TIFFGetField( info->tiff , TIFFTAG_PLANARCONFIG , &temp16 );
	if( temp16!=PLANARCONFIG_CONTIG ) fprintf( stderr , "Planar configuration must be contiguous: %d != %d\n" , temp16 , PLANARCONFIG_CONTIG );
	info->data = AllocPointer< byte >( scanLineSize1 );
	return info;
}

template< int Channels , class ChannelType >
void TIFReadRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	TIFFInfo* info = (TIFFInfo*)v;
	TIFFReadScanline( info->tiff , GetAddress( info->data ) , j );

	if( info->bitsPerSample==8 ) ConvertRow< unsigned char , ChannelType >( ( ConstPointer(unsigned char) )info->data , pixels , info->width , info->samplesPerPixel , Channels );
	else                         ConvertRow< uint16_t      , ChannelType >( ( ConstPointer(uint16_t     ) )info->data , pixels , info->width , info->samplesPerPixel , Channels );
}
inline void TIFFinalizeRead( void* v )
{
	TIFFInfo* info = (TIFFInfo*)v;
	TIFFClose( info->tiff );
	FreePointer( info->data );
	free(info);
}
template< int Channels , bool HDR >
void* TIFInitWrite( char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	TIFFInfo* info = (TIFFInfo*)malloc(sizeof(TIFFInfo));

	info->tiff=TIFFOpen(fileName,"w");
	if (!info->tiff)																fprintf(stderr,"TIFFOpen failed\n")								,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGEWIDTH, width))						fprintf (stderr, "Can't set ImageWidth tag.\n")					,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGELENGTH, height))						fprintf (stderr, "Can't set ImageLength tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_BITSPERSAMPLE, HDR?16:8))					fprintf (stderr, "Can't set BitsPerSample tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_SAMPLESPERPIXEL, Channels))				fprintf (stderr, "Can't set SamplesPerPixel tag.\n")			,	exit(0);
	if( !quality )
	{
		if (!TIFFSetField(info->tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE))		fprintf (stderr, "Can't set Compression tag.\n")				,	exit(0);
	}
	else
	{
		if (!TIFFSetField(info->tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW))		fprintf (stderr, "Can't set Compression tag.\n")				,	exit(0);
	}
	if (!TIFFSetField(info->tiff, TIFFTAG_ROWSPERSTRIP, 1))							fprintf (stderr, "Can't set RowsPerStrip tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG))		fprintf (stderr, "Can't set PlanarConfiguration tag.\n")		,	exit(0);
	if( Channels==1 )
	{
		if (!TIFFSetField(info->tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK))	fprintf (stderr, "Can't set PhotometricInterpretation tag.\n")	,	exit(0);
	}
	else
	{
		if (!TIFFSetField(info->tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB))		fprintf (stderr, "Can't set PhotometricInterpretation tag.\n")	,	exit(0);
	}
	if (!TIFFSetField(info->tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE))			fprintf (stderr, "Can't set ResolutionUnit tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT))		fprintf (stderr, "Can't set Orientation tag.\n")				,	exit(0);

	info->width = width;
	info->bitsPerSample = HDR ? 16 : 8;
	info->samplesPerPixel = Channels;

	int scanLineSize1 = TIFFScanlineSize( info->tiff ) , scanLineSize2 = ( info->bitsPerSample / 8 ) * info->samplesPerPixel * info->width;
	if( scanLineSize1!=scanLineSize2 ) fprintf( stderr , "Scanline sizes do not agree: %d!=%d\n" , scanLineSize1 , scanLineSize2 ) , exit( 0 );
	info->data = AllocPointer< byte >( scanLineSize1 );

	return info;
}
template< int Channels , class ChannelType >
void TIFWriteRow( Pointer( ChannelType ) pixels , void* v , int j )
{
	TIFFInfo* info = (TIFFInfo*)v;
	if( info->bitsPerSample==8 ) ConvertRow< ChannelType , unsigned char >( pixels , ( Pointer(unsigned char) )info->data , info->width , Channels , info->samplesPerPixel );
	else                         ConvertRow< ChannelType , uint16_t      >( pixels , ( Pointer(uint16_t     ) )info->data , info->width , Channels , info->samplesPerPixel );
	TIFFWriteScanline( info->tiff , GetAddress( info->data ) , j , 0 );
}
inline void TIFFinalizeWrite(void* v)
{
	TIFFInfo* info = (TIFFInfo*)v;
	TIFFClose( info->tiff );
	FreePointer( info->data );
	free( info );
}
