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
#ifndef IMAGE_STREAM_INCLUDED
#define IMAGE_STREAM_INCLUDED

static char JPEGFileExtension[] = "jpg";
static int DefaultOutputTileWidth = 8192;
static int DefaultOutputTileHeight = 8192;
static const char* DefaultOutputTileExtension = JPEGFileExtension;
static double GAMMA = 2.2;

#include <math.h>
#include <stdint.h>
#include "Util/Half/half.h"
#include "Util/GridStream.h"
#include "Util/Color.h"
#include "RAWStream.inl"
#include "BMPStream.inl"
#include "KROStream.inl"
#include "JPGStream.inl"
#include "PNGStream.inl"
#include "PFMStream.inl"
#include "PGMStream.inl"
#ifndef NO_TIFF_SUPPORT
#include "TIFStream.inl"
#endif // NO_TIFF_SUPPORT
#include <omp.h>


template< class ChannelType >
struct ImageReadInfo
{
	void* (*InitRead    ) ( char* , int& , int& , MultiStreamIOServer* );
	void  (*ReadRow     ) ( Pointer( ChannelType ) , void* , int);
	void  (*FinalizeRead) ( void* );
	ImageReadInfo( void* (*I)( char* , int& , int& , MultiStreamIOServer* ) , void (*R)( Pointer( ChannelType ) , void* , int ) , void (*F)( void* ) ){ InitRead = I , ReadRow = R , FinalizeRead = F; }
};
template< class ChannelType >
struct ImageWriteInfo
{
	void* (*InitWrite    ) ( char* , int , int , int , MultiStreamIOServer* );
	void  (*WriteRow     ) ( Pointer( ChannelType ) , void* , int );
	void  (*FinalizeWrite) ( void* );
	ImageWriteInfo( void* (*I)( char* , int , int , int , MultiStreamIOServer* ) , void (*W)( Pointer( ChannelType ) , void* , int ) , void (*F)( void* ) ){ InitWrite = I , WriteRow = W , FinalizeWrite = F; }
};

inline void RAWdoubleGetImageInfo( char* fn , int& w , int& h , int& c , int& bpc ){ return RAWGetImageInfo( fn , w , h , c , bpc ); }
inline void RAWfloatGetImageInfo ( char* fn , int& w , int& h , int& c , int& bpc ){ return RAWGetImageInfo( fn , w , h , c , bpc ); }
inline void RAWhalfGetImageInfo  ( char* fn , int& w , int& h , int& c , int& bpc ){ return RAWGetImageInfo( fn , w , h , c , bpc ); }
inline void RAWint32GetImageInfo ( char* fn , int& w , int& h , int& c , int& bpc ){ return RAWGetImageInfo( fn , w , h , c , bpc ); }
inline void RAWint16GetImageInfo ( char* fn , int& w , int& h , int& c , int& bpc ){ return RAWGetImageInfo( fn , w , h , c , bpc ); }
inline void RAWint8GetImageInfo  ( char* fn , int& w , int& h , int& c , int& bpc ){ return RAWGetImageInfo( fn , w , h , c , bpc ); }

template< int Channels , bool HDR > void* RAWdoubleInitWrite( char* fn , int w , int h , int q , MultiStreamIOServer* io ){ return RAWInitWrite< Channels >( RAWInfo::CHANNEL_DOUBLE  , fn , w , h , q , io ); }
template< int Channels , bool HDR > void* RAWfloatInitWrite ( char* fn , int w , int h , int q , MultiStreamIOServer* io ){ return RAWInitWrite< Channels >( RAWInfo::CHANNEL_FLOAT   , fn , w , h , q , io ); }
template< int Channels , bool HDR > void* RAWhalfInitWrite  ( char* fn , int w , int h , int q , MultiStreamIOServer* io ){ return RAWInitWrite< Channels >( RAWInfo::CHANNEL_HALF    , fn , w , h , q , io ); }
template< int Channels , bool HDR > void* RAWint32InitWrite ( char* fn , int w , int h , int q , MultiStreamIOServer* io ){ return RAWInitWrite< Channels >( RAWInfo::CHANNEL_UINT_32 , fn , w , h , q , io ); }
template< int Channels , bool HDR > void* RAWint16InitWrite ( char* fn , int w , int h , int q , MultiStreamIOServer* io ){ return RAWInitWrite< Channels >( RAWInfo::CHANNEL_UINT_16 , fn , w , h , q , io ); }
template< int Channels , bool HDR > void* RAWint8InitWrite  ( char* fn , int w , int h , int q , MultiStreamIOServer* io ){ return RAWInitWrite< Channels >( RAWInfo::CHANNEL_UINT_8  , fn , w , h , q , io ); }

inline void* RAWdoubleInitRead( char* fn , int& w , int& h , MultiStreamIOServer* io ){ return RAWInitRead( fn , w , h , io ); }
inline void* RAWfloatInitRead ( char* fn , int& w , int& h , MultiStreamIOServer* io ){ return RAWInitRead( fn , w , h , io ); }
inline void* RAWhalfInitRead  ( char* fn , int& w , int& h , MultiStreamIOServer* io ){ return RAWInitRead( fn , w , h , io ); }
inline void* RAWint32InitRead ( char* fn , int& w , int& h , MultiStreamIOServer* io ){ return RAWInitRead( fn , w , h , io ); }
inline void* RAWint16InitRead ( char* fn , int& w , int& h , MultiStreamIOServer* io ){ return RAWInitRead( fn , w , h , io ); }
inline void* RAWint8InitRead  ( char* fn , int& w , int& h , MultiStreamIOServer* io ){ return RAWInitRead( fn , w , h , io ); }

template< int Channels , class ChannelType > void RAWdoubleWriteRow( Pointer(ChannelType) pixels , void* info , int j ){ return RAWWriteRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWfloatWriteRow ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWWriteRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWhalfWriteRow  ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWWriteRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWint32WriteRow ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWWriteRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWint16WriteRow ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWWriteRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWint8WriteRow  ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWWriteRow< Channels , ChannelType >( pixels , info , j ); }

template< int Channels , class ChannelType > void RAWdoubleReadRow( Pointer(ChannelType) pixels , void* info , int j ){ return RAWReadRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWfloatReadRow ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWReadRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWhalfReadRow  ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWReadRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWint32ReadRow ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWReadRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWint16ReadRow ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWReadRow< Channels , ChannelType >( pixels , info , j ); }
template< int Channels , class ChannelType > void RAWint8ReadRow  ( Pointer(ChannelType) pixels , void* info , int j ){ return RAWReadRow< Channels , ChannelType >( pixels , info , j ); }

inline void RAWdoubleFinalizeWrite( void* info ){ return RAWFinalizeWrite( info ); }
inline void RAWfloatFinalizeWrite ( void* info ){ return RAWFinalizeWrite( info ); }
inline void RAWhalfFinalizeWrite  ( void* info ){ return RAWFinalizeWrite( info ); }
inline void RAWint32FinalizeWrite ( void* info ){ return RAWFinalizeWrite( info ); }
inline void RAWint16FinalizeWrite ( void* info ){ return RAWFinalizeWrite( info ); }
inline void RAWint8FinalizeWrite  ( void* info ){ return RAWFinalizeWrite( info ); }

inline void RAWdoubleFinalizeRead( void* info ){ return RAWFinalizeRead( info ); }
inline void RAWfloatFinalizeRead ( void* info ){ return RAWFinalizeRead( info ); }
inline void RAWhalfFinalizeRead  ( void* info ){ return RAWFinalizeRead( info ); }
inline void RAWint32FinalizeRead ( void* info ){ return RAWFinalizeRead( info ); }
inline void RAWint16FinalizeRead ( void* info ){ return RAWFinalizeRead( info ); }
inline void RAWint8FinalizeRead  ( void* info ){ return RAWFinalizeRead( info ); }

#define DECLARE_IMAGE_TYPE( TYPE ) \
	static ImageReadInfo < ChannelType > TYPE ##     ReadInfo( void ); \
	static ImageWriteInfo< ChannelType > TYPE ##    WriteInfo( void ); \
	static ImageWriteInfo< ChannelType > TYPE ## HDRWriteInfo( void );

#define DEFINE_IMAGE_TYPE( TYPE ) \
	template< class ChannelType , int Channels > ImageReadInfo < ChannelType > ImageReadWriteInfo< ChannelType , Channels >::TYPE ##     ReadInfo( void ){ return ImageReadInfo < ChannelType >( TYPE ## InitRead                      , TYPE ##  ReadRow< Channels , ChannelType > , TYPE ##  FinalizeRead ); } \
	template< class ChannelType , int Channels > ImageWriteInfo< ChannelType > ImageReadWriteInfo< ChannelType , Channels >::TYPE ##    WriteInfo( void ){ return ImageWriteInfo< ChannelType >( TYPE ## InitWrite< Channels , false > , TYPE ## WriteRow< Channels , ChannelType > , TYPE ## FinalizeWrite ); } \
	template< class ChannelType , int Channels > ImageWriteInfo< ChannelType > ImageReadWriteInfo< ChannelType , Channels >::TYPE ## HDRWriteInfo( void ){ return ImageWriteInfo< ChannelType >( TYPE ## InitWrite< Channels , true  > , TYPE ## WriteRow< Channels , ChannelType > , TYPE ## FinalizeWrite ); } \
	static void TYPE ## GetImageSize( char* fileName , int& width , int& height )       \
	{                                                                                   \
		int channels , bytesPerChannel;                                                 \
		TYPE ## GetImageInfo( fileName , width , height , channels , bytesPerChannel ); \
	}

template< class ChannelType , int Channels >
class ImageReadWriteInfo
{
public:
	DECLARE_IMAGE_TYPE( RAWdouble );
	DECLARE_IMAGE_TYPE( RAWfloat );
	DECLARE_IMAGE_TYPE( RAWhalf );
	DECLARE_IMAGE_TYPE( RAWint32 );
	DECLARE_IMAGE_TYPE( RAWint16 );
	DECLARE_IMAGE_TYPE( RAWint8 );
	DECLARE_IMAGE_TYPE( BMP );
	DECLARE_IMAGE_TYPE( KRO );
	DECLARE_IMAGE_TYPE( JPG );
	DECLARE_IMAGE_TYPE( PNG );
	DECLARE_IMAGE_TYPE( PFM );
	DECLARE_IMAGE_TYPE( PGM );
#ifndef NO_TIFF_SUPPORT
	DECLARE_IMAGE_TYPE( TIF );
#endif // NO_TIFF_SUPPORT
};
DEFINE_IMAGE_TYPE( RAWdouble );
DEFINE_IMAGE_TYPE( RAWfloat );
DEFINE_IMAGE_TYPE( RAWhalf );
DEFINE_IMAGE_TYPE( RAWint32 );
DEFINE_IMAGE_TYPE( RAWint16 );
DEFINE_IMAGE_TYPE( RAWint8 );
DEFINE_IMAGE_TYPE( BMP );
DEFINE_IMAGE_TYPE( KRO );
DEFINE_IMAGE_TYPE( JPG );
DEFINE_IMAGE_TYPE( PNG );
DEFINE_IMAGE_TYPE( PFM );
DEFINE_IMAGE_TYPE( PGM );
#ifndef NO_TIFF_SUPPORT
DEFINE_IMAGE_TYPE( TIF );
#endif // NO_TIFF_SUPPORT
#undef DECLARE_IMAGE_TYPE
#undef DEFINE_IMAGE_TYPE


template< class ChannelType , int Channels >
class WriteImageStream : public StreamingGrid
{
	Pointer( ChannelType ) _pixels;

	int current;
	bool _clamp , _gammaEncode;
	int _q , _w , _h;
	void _Init( char* fileName , MultiStreamIOServer* ioServer );

	void* (*_InitWrite)		( char* , int , int , int , MultiStreamIOServer* );
	void  (*_WriteRow)		( Pointer(ChannelType) , void* , int );
	void  (*_FinalizeWrite)	( void* );
protected:
	void* info;
public:
	~WriteImageStream( void );
	WriteImageStream ( void );
	WriteImageStream ( ImageWriteInfo< ChannelType > iwInfo , char* fileName , int width , int height , bool gammaEncode , int quality , bool clamp , MultiStreamIOServer* ioServer );
	void Init        ( ImageWriteInfo< ChannelType > iwInfo , char* fileName , int width , int height , bool gammaEncode , int quality , bool clamp , MultiStreamIOServer* ioServer );
	int		rows			(void) const;
	int		rowSize			(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance			(void);
};

template< class ChannelType , int Channels >
class ReadImageStream : public StreamingGrid
{
	Pointer( ChannelType ) _pixels;

	bool _gammaDecode;
	int current;
	int _w , _h;

	void* (*_InitRead)		( char* , int& , int& , MultiStreamIOServer* );
	void  (*_ReadRow)		( Pointer(ChannelType) , void* , int);
	void  (*_FinalizeRead)	( void* );
protected:
	void* info;
public:
	~ReadImageStream( void );
	ReadImageStream ( void );
	ReadImageStream	( ImageReadInfo< ChannelType > irInfo , char* fileName , int &width , int &height , bool gammaDecode , MultiStreamIOServer* ioServer );
	void Init       ( ImageReadInfo< ChannelType > irInfo , char* fileName , int &width , int &height , bool gammaDecode , MultiStreamIOServer* ioServer );
	int		rows			(void) const;
	int		rowSize			(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance			(void);
};

template< class Real , int Channels >
class StreamingInputTile
{
	StreamingGrid* data;
	int w,h,rc;
	int sX,sY;
	char tileName[1024];
	Real** rows;
	double average[Channels];
public:
	StreamingInputTile(void);
	~StreamingInputTile(void);

//	void Init( const char* tName , int rowCount=1 , int startX=0 , int startY=0 , bool gammaDecode=false );
	void Init( const char* tName , int rowCount , int startX , int startY , bool gammaDecode );
	void init( int idx );
	Color< Real , Channels > operator() ( int i , int j );
	Color< Real , Channels > operator() ( int i , int j , Real& a );
	int width(void) const;
	int height(void) const;
	int startX(void) const;
	int startY(void) const;
};

template< class Real , int Channels >
class RGridOfImages : public StreamingGrid
{
	int _c , _r , _w , _h;
	int *_widths , *_heights;
	int _current;
	bool _gammaDecode;
	int _threads;
	StreamingGrid*** _grid;
	char*** _fileNames;
	Pointer( Real ) _buffer;
	MultiStreamIOServer* _server;
public:
	static void GetImageSize( char* fileName , int& width , int& height );
	RGridOfImages		(char* fileName , int& width , int& height , bool gammaDecode , MultiStreamIOServer* ioServer , int threads );
	~RGridOfImages		(void);
	int		rows		(void) const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
};
template< class Real , int Channels >
class WGridOfImages : public StreamingGrid
{
	int _c,_r,_w,_h;
	int *_widths,*_heights;
	int _current,_quality;
	bool _hdr , _gammaEncode;
	int _threads;
	StreamingGrid*** _grid;
	char*** _fileNames;
	Pointer( Real ) _buffer;
	MultiStreamIOServer* _server;
	void _init(char* fileName,const char* header,const char* ext,int width,int height , bool gammaEncode ,int quality,int tileWidth,int tileHeight , MultiStreamIOServer* ioServer , int threads , bool hdr );
public:
	WGridOfImages		( char* fileName ,                                        int width , int height , bool gammaEncode , int quality , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , int threads , bool hdr );
	WGridOfImages		( char* fileName , const char* header , const char* ext , int width , int height , bool gammaEncode , int quality , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , int threads , bool hdr );
	~WGridOfImages		(void);
	int		rows		(void) const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
};

template< class Real , int Channels >
class WImageClipper : public StreamingGrid
{
	bool deleteGrid;
	int _w,_h,_ox,_oy,_ow,_oh;
	int current;
	bool _hdr;
	StreamingGrid* grid;
	Pointer( Real ) buffer;
public:
	WImageClipper		( char* fileName , int width , int height , bool gammaEncode , int outOffX , int outOffY , int outWidth , int outHeight , int quality , MultiStreamIOServer* server , bool hdr , int threads );
	WImageClipper		( StreamingGrid* outGrid , int width , int height , int outOffX , int outOffY , int outWidth , int outHeight , bool hdr );
	~WImageClipper		(void);
	int		rows		(void)	const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
};

template< class Real , int Channels >
class RChannelExtractor : public StreamingGrid
{
	int _c , _w , _h;
	int current;
	StreamingGrid* grid;
	Real* buffer;
public:
	RChannelExtractor	( char* fileName , int& width , int& height , bool gammaDecode , int channel );
	~RChannelExtractor	(void);
	int		rows		(void)	const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
};
template< class Real , int Channels >
class RImageSampler : public StreamingGrid
{
	int _inW , _inH , _outW , _outH , _inCurrent , _outCurrent;
	bool _nearest;
	StreamingGrid *_grid;
	Pointer( Real ) _inBuffers[2];
	Pointer( Real ) _outBuffer;
public:
	RImageSampler( char* fileName , int width , int height , bool nearest , bool gammaDecode , MultiStreamIOServer* ioServer , int threads );
	~RImageSampler( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void advance( void );
};
template< class PixelType , class LabelType , int PChannels , int LChannels >
class RMaskedImageSampler : public StreamingGrid
{
	int _inW , _inH , _outW , _outH , _inCurrent , _outCurrent;
	StreamingGrid *_pixelGrid , *_labelGrid;
	Pointer( PixelType ) _pixelBuffers[2];
	Pointer( LabelType ) _labelBuffers[2];
	Pointer( PixelType ) _outBuffer;
public:
	RMaskedImageSampler( char* pixelName , char* labelName , int width , int height , bool gammaDecode , MultiStreamIOServer* ioServer , int threads );
	~RMaskedImageSampler( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void advance( void );
};
template< class Real , int Channels >
class WImageSampler : public StreamingGrid
{
	int _inW , _inH , _outW , _outH , _inCurrent , _outCurrent;
	bool _nearest;
	StreamingGrid* _grid;
	Pointer( Real ) _inBuffers[2];
public:
	WImageSampler( char* fileName , int inWidth , int inHeight , int outWidth , int outHeight , bool nearest , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , int threads );
	~WImageSampler( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void advance( void );
};
inline void GetReadSize( char* fileName , int& width , int& height );
template< class Real , int Channels > StreamingGrid* GetReadStream ( char* fileName ,       int& width ,       int& height , bool gammaDecode ,                      MultiStreamIOServer* ioServer ,            int threads );
template< class Real , int Channels > StreamingGrid* GetWriteStream( char* fileName , const int& width , const int& height , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , bool hdr , int threads );

#include "ImageStream.inl"
#endif // GRID_STREAM_INCLUDED