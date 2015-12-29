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
#ifndef CHANNEL_CONVERTER_INCLUDED
#define CHANNEL_CONVERTER_INCLUDED

#include <stdint.h>
#include "Util/Half/half.h"

template< class Real >	bool IsFloatType( void );
template< > inline		bool IsFloatType< half >	( void ){ return true;  }
template< > inline		bool IsFloatType< float >	( void ){ return true;  }
template< > inline		bool IsFloatType< double >	( void ){ return true;  }
template< class Real >	bool IsFloatType			( void ){ return false; }

template< class Real >
class ConversionFactor
{
public:
	static const double Factor;
	static const double FactorR;
	static const double Offset;
};
template<> const double ConversionFactor<          char >::Factor = (double)( ( (long long)(1)<< 8)-1 );
template<> const double ConversionFactor< unsigned char >::Factor = (double)( ( (long long)(1)<< 8)-1 );
template<> const double ConversionFactor<       int16_t >::Factor = (double)( ( (long long)(1)<<16)-1 );
template<> const double ConversionFactor<      uint16_t >::Factor = (double)( ( (long long)(1)<<16)-1 );
template<> const double ConversionFactor<           int >::Factor = (double)( ( (long long)(1)<<32)-1 );
template<> const double ConversionFactor<  unsigned int >::Factor = (double)( ( (long long)(1)<<32)-1 );
template<> const double ConversionFactor<        double >::Factor = 1.;
template<> const double ConversionFactor<          half >::Factor = 1.;
template<> const double ConversionFactor<         float >::Factor = 1.;
template<> const double ConversionFactor<          char >::FactorR = 1./(double)( ( (long long)(1)<< 8)-1 );
template<> const double ConversionFactor< unsigned char >::FactorR = 1./(double)( ( (long long)(1)<< 8)-1 );
template<> const double ConversionFactor<       int16_t >::FactorR = 1./(double)( ( (long long)(1)<<16)-1 );
template<> const double ConversionFactor<      uint16_t >::FactorR = 1./(double)( ( (long long)(1)<<16)-1 );
template<> const double ConversionFactor<           int >::FactorR = 1./(double)( ( (long long)(1)<<32)-1 );
template<> const double ConversionFactor<  unsigned int >::FactorR = 1./(double)( ( (long long)(1)<<32)-1 );
template<> const double ConversionFactor<        double >::FactorR = 1.;
template<> const double ConversionFactor<          half >::FactorR = 1.;
template<> const double ConversionFactor<         float >::FactorR = 1.;
template<> const double ConversionFactor<          char >::Offset = 0.5;
template<> const double ConversionFactor< unsigned char >::Offset = 0.5;
template<> const double ConversionFactor<       int16_t >::Offset = 0.5;
template<> const double ConversionFactor<      uint16_t >::Offset = 0.5;
template<> const double ConversionFactor<           int >::Offset = 0.5;
template<> const double ConversionFactor<  unsigned int >::Offset = 0.5;
template<> const double ConversionFactor<        double >::Offset = 0.0;
template<> const double ConversionFactor<          half >::Offset = 0.0;
template<> const double ConversionFactor<         float >::Offset = 0.0;


template< class In > inline double ChannelToDouble( const In& in ){ return ( (double)in + ConversionFactor< In >::Offset ) * ConversionFactor< In >::FactorR; }
template< class Out > inline Out ChannelFromDouble( double in ){ return (Out)( in * ConversionFactor< Out >::Factor ); }

template< class In , class Out >
class Converter
{
public:
	static const double ConvertFactor;
	static const double ConvertOffset;
};
template< class In , class Out > const double Converter< In , Out >::ConvertFactor = ConversionFactor< Out >::Factor * ConversionFactor< In >::FactorR;
template< class In , class Out > const double Converter< In , Out >::ConvertOffset = ConversionFactor< Out >::Offset;

// [WARNING] When using integer indices, this definition of a channel converter allocates slightly smaller bins for the first and last buckets. 
// The alternative is to use half open intervals of the same length, but this would require explicitly performing bounds testing.
template< class In , class Out > inline Out ConvertChannel( const In& in ) { return Out( double( in )*Converter< In , Out >::ConvertFactor + Converter< In , Out >::ConvertOffset ); }

template< class InChannelType , class OutChannelType >
void ConvertRow( ConstPointer( InChannelType ) inRow , Pointer( OutChannelType ) outRow , int width , int inChannels , int outChannels )
{
	if(  inChannels!=1 &&  inChannels!=3 &&  inChannels!=4 ) fprintf( stderr , "[ERROR] ConvertRow: only Gray, RGB, and RGBA supported: %d\n" ,  inChannels ) , exit( 0 );
	if( outChannels!=1 && outChannels!=3 && outChannels!=4 ) fprintf( stderr , "[ERROR] ConvertRow: only Gray, RGB, and RGBA supported: %d\n" , outChannels ) , exit( 0 );
	if( inChannels==outChannels ) for( int i=0 ; i<width*inChannels ; i++ ) outRow[i] = ConvertChannel< InChannelType , OutChannelType >( inRow[i] );
	else if( outChannels==1 )
	{
		if( inChannels==3 || inChannels==4 )
		{
			for( int i=0 ; i<width ; i++ ) 
			{
				double sum = ConvertChannel< InChannelType , double >( inRow[i*inChannels] ) + ConvertChannel< InChannelType , double >( inRow[i*inChannels+1] ) + ConvertChannel< InChannelType , double >( inRow[i*inChannels+2] );
				outRow[i] = ConvertChannel< double , OutChannelType >( sum/3. );
			}
		}
	}
	else if( outChannels==3 )
	{
		if     ( inChannels==1 ) for( int i=0 ; i<width ; i++ ) outRow[i*3+0] = outRow[i*3+1] = outRow[i*3+2] = ConvertChannel< InChannelType , OutChannelType >( inRow[i] );
		else if( inChannels==4 ) for( int i=0 ; i<width ; i++ ) for( int c=0 ; c<3 ; c++ ) outRow[i*3+c] = ConvertChannel< InChannelType , OutChannelType >( inRow[i*4+c] );
	}
	else if( outChannels==4 )
	{
		if     ( inChannels==1 ) for( int i=0 ; i<width ; i++ ) outRow[i*4+0] = outRow[i*4+1] = outRow[i*4+2] = ConvertChannel< InChannelType , OutChannelType >( inRow[i] ) , outRow[i*4+3] = ConvertChannel< double , OutChannelType >( 1. );
		else if( inChannels==3 ) for( int i=0 ; i<width ; i++ ) for( int c=0 ; c<3 ; c++ ) outRow[i*4+c] = ConvertChannel< InChannelType , OutChannelType >( inRow[i*3+c] ) , outRow[i*4+3] = ConvertChannel< double , OutChannelType >( 1. );;
	}
}

#endif // CHANNEL_CONVERTER_INCLUDED