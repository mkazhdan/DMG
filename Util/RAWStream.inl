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
#include <Util/BaseMultiStreamIO.h>
#include "ChannelConverter.h"

struct RAWInfo
{
	enum
	{
		CHANNEL_DOUBLE ,
		CHANNEL_FLOAT ,
		CHANNEL_HALF ,
		CHANNEL_UINT_32 ,
		CHANNEL_UINT_16 ,
		CHANNEL_UINT_8
	};
	FILE* fp;
	int channels , channelType , width;
	Pointer( byte ) data;

	static size_t TypeSize( int type )
	{
		switch( type )
		{
			case CHANNEL_DOUBLE:  return sizeof( double        );
			case CHANNEL_FLOAT:   return sizeof( float         );
			case CHANNEL_HALF:    return sizeof( half          );
			case CHANNEL_UINT_32: return sizeof( unsigned int  );
			case CHANNEL_UINT_16: return sizeof( uint16_t      );
			case CHANNEL_UINT_8:  return sizeof( unsigned char );
			default: fprintf( stderr , "[ERROR] RAWInfo::TypeSize: unrecognized channel type: %d\n" , type ) , exit( 0 );
		}
	}
};

inline void RAWGetImageInfo( char* fileName , int& width , int& height , int& channels , int& bytesPerChannel )
{
	FILE* fp;
	fp = fopen( fileName , "rb" );
	if( !fp )	fprintf( stderr , "[ERROR] GetImageSize: Failed to open: %s\n" , fileName ) , exit(0);
	int type;
	fread( &width    , sizeof(int) , 1 , fp );
	fread( &height   , sizeof(int) , 1 , fp );
	fread( &channels , sizeof(int) , 1 , fp );
	fread( &type     , sizeof(int) , 1 , fp );
	fclose( fp );
	bytesPerChannel = (int)RAWInfo::TypeSize( type );
}

template< int Channels >
void* RAWInitWrite( int channelType , char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	RAWInfo* info = ( RAWInfo* )malloc( sizeof( RAWInfo ) );
	info->width = width;
	info->channelType = channelType;
	info->channels = Channels;
	info->fp = fopen( fileName , "wb" );
	if( !info->fp ) fprintf( stderr , "[ERROR] InitWrite: Failed to open: %s\n" , fileName ) , exit(0);

	fwrite( &width       , sizeof(int) , 1 , info->fp );
	fwrite( &height      , sizeof(int) , 1 , info->fp );
	fwrite( &info->channels    , sizeof(int) , 1 , info->fp );
	fwrite( &info->channelType , sizeof(int) , 1 , info->fp );
	info->data = AllocPointer< byte >( RAWInfo::TypeSize( info->channelType ) * Channels * width );
	return info;
}

template< int Channels , class ChannelType >
void RAWWriteRow( Pointer(ChannelType) pixels , void* v , int j )
{
	RAWInfo* info = ( RAWInfo* )v;
	switch( info->channelType )
	{
	case RAWInfo::CHANNEL_DOUBLE:  ConvertRow< ChannelType , double        >( pixels , ( Pointer(double       ) )info->data , info->width , Channels , info->channels ) ; break;
	case RAWInfo::CHANNEL_FLOAT:   ConvertRow< ChannelType , float         >( pixels , ( Pointer(float        ) )info->data , info->width , Channels , info->channels ) ; break;
	case RAWInfo::CHANNEL_HALF:    ConvertRow< ChannelType , half          >( pixels , ( Pointer(half         ) )info->data , info->width , Channels , info->channels ) ; break;
	case RAWInfo::CHANNEL_UINT_32: ConvertRow< ChannelType , unsigned int  >( pixels , ( Pointer(unsigned int ) )info->data , info->width , Channels , info->channels ) ; break;
	case RAWInfo::CHANNEL_UINT_16: ConvertRow< ChannelType , uint16_t      >( pixels , ( Pointer(uint16_t     ) )info->data , info->width , Channels , info->channels ) ; break;
	case RAWInfo::CHANNEL_UINT_8:  ConvertRow< ChannelType , unsigned char >( pixels , ( Pointer(unsigned char) )info->data , info->width , Channels , info->channels ) ; break;
	default: fprintf( stderr , "[ERROR] RAWWriteRow: unrecognized channel type: %d\n" , info->channelType ) , exit( 0 );
	}
	fwrite( info->data , RAWInfo::TypeSize( info->channelType ) , info->channels * info->width , info->fp );
}

inline void RAWFinalizeWrite( void* v )
{
	RAWInfo* info = ( RAWInfo* )v;
	fclose( info->fp );
	FreePointer( info->data );
	free( info );
}
inline void* RAWInitRead( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	RAWInfo* info = ( RAWInfo* )malloc( sizeof( RAWInfo ) );
	info->fp=fopen( fileName , "rb" );
	if( !info->fp ) fprintf( stderr , "[ERROR] InitRead: Failed to open: %s\n" , fileName ) , exit(0);
	fread( &width  , sizeof(int) , 1 , info->fp );
	fread( &height , sizeof(int) , 1 , info->fp );
	fwrite( &info->channels    , sizeof(int) , 1 , info->fp );
	fwrite( &info->channelType , sizeof(int) , 1 , info->fp );
	info->width = width;
	info->data = AllocPointer< byte >( RAWInfo::TypeSize( info->channelType ) * info->channels * width );
	return info;
}
template< int Channels , class ChannelType >
void RAWReadRow( Pointer(ChannelType) pixels , void* v , int j )
{
	RAWInfo* info = ( RAWInfo* )v;

	fread( info->data , RAWInfo::TypeSize( info->channelType ) , info->channels*info->width , info->fp );
	switch( info->channelType )
	{
	case RAWInfo::CHANNEL_DOUBLE:  ConvertRow< double        , ChannelType >( ( ConstPointer(double       ) )info->data , pixels , info->width , info->channels , Channels ) ; break;
	case RAWInfo::CHANNEL_FLOAT:   ConvertRow< float         , ChannelType >( ( ConstPointer(float        ) )info->data , pixels , info->width , info->channels , Channels ) ; break;
	case RAWInfo::CHANNEL_HALF:    ConvertRow< half          , ChannelType >( ( ConstPointer(half         ) )info->data , pixels , info->width , info->channels , Channels ) ; break;
	case RAWInfo::CHANNEL_UINT_32: ConvertRow< unsigned int  , ChannelType >( ( ConstPointer(unsigned int ) )info->data , pixels , info->width , info->channels , Channels ) ; break;
	case RAWInfo::CHANNEL_UINT_16: ConvertRow< uint16_t      , ChannelType >( ( ConstPointer(uint16_t     ) )info->data , pixels , info->width , info->channels , Channels ) ; break;
	case RAWInfo::CHANNEL_UINT_8:  ConvertRow< unsigned char , ChannelType >( ( ConstPointer(unsigned char) )info->data , pixels , info->width , info->channels , Channels ) ; break;
	default: fprintf( stderr , "[ERROR] RAWReadRow: unrecognized channel type: %d\n" , info->channelType ) , exit( 0 );
	}
}

inline void RAWFinalizeRead( void* v )
{
	RAWInfo* info = ( RAWInfo* )v;
	fclose( info->fp );
	FreePointer( info->data );
	free( info );
}
