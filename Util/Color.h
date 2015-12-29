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
#ifndef COLOR_INCLUDED
#define COLOR_INCLUDED

template< class Real , int Channels >
class Color
{
	Real c[Channels];
public:
	Color( void );
	Color( Real clr );
	Color( const Real* clr );
	template< class Real2 >
	Color( const Color< Real2 , Channels >& clr );

//	Real luminance( void ) const;

	Real& operator[] (int i);
	const Real& operator[] (int i) const;

	Color operator + ( const Color& clr ) const;
	Color operator - ( const Color& clr ) const;
	Color operator - ( void             ) const;
	Color operator * ( const Real& s    ) const;
	Color operator / ( const Real& s    ) const;
	// An essential method for getting the conjugate gradient solver to work
//	Real  operator * ( const Color& clr ) const;
	Color  operator * ( const Color& clr ) const;

	Color& operator += ( const Color& clr );
	Color& operator -= ( const Color& clr );
	Color& operator *= ( const Real& s );
	Color& operator /= ( const Real& s );
	bool operator == ( const Color< Real , Channels >& clr ) const;
	bool operator != ( const Color< Real , Channels >& clr ) const;
	bool operator <  ( const Color< Real , Channels >& clr ) const;
	bool operator <= ( const Color< Real , Channels >& clr ) const;
	bool operator >  ( const Color< Real , Channels >& clr ) const;
	bool operator >= ( const Color< Real , Channels >& clr ) const;
};
template< class Real , int Channels >
inline Color< Real , Channels > pow( const Color< Real , Channels >& in , double e )
{
	Color< Real , Channels > out;
	for( int c=0 ; c<Channels ; c++ ) out[c] = (Real)pow( (double)in[c] , e );
	return out;
}

#include "Color.inl"
#endif // COLOR_INCLUDED
