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
template< class Real , int Channels > Color< Real , Channels >::Color( void ){ memset( c , 0 , sizeof( Real ) * Channels ); }
template< class Real , int Channels > Color< Real , Channels >::Color( Real clr ){ c[0] = c[1] = c[2] = clr; }
template< class Real , int Channels > Color< Real , Channels >::Color( const Real* clr ){ memcpy( c , clr , sizeof(Real)*Channels ); }

//template< class Real , int Channels > Color< Real , Channels >::Color(Real c0,Real c1,Real c2)	{ c[0]=c0;	c[1]=c1;	c[2]=c2; }

template< class Real , int Channels >
template< class Real2 >
Color< Real , Channels >::Color( const Color< Real2 , Channels >& clr ){ for( int c=0 ; c<Channels ; c++ ) this->c[c] = (Real)clr[c]; }

//template<class Real>
//Real Color< Real , Channels >::luminance(void) const		{return Real(c[0]*0.3+c[1]*0.59+c[2]*0.11);}

template< class Real , int Channels >
Real& Color< Real , Channels >::operator[] (int i) {return c[i];}
template< class Real , int Channels >
const Real& Color< Real , Channels >::operator[] (int i) const {return c[i];}

template< class Real , int Channels > Color< Real , Channels > Color< Real , Channels >::operator + ( const Color& clr ) const { Color _clr; for( int c=0 ; c<Channels ; c++ ) _clr[c] = this->c[c]+clr.c[c] ; return _clr; };
template< class Real , int Channels > Color< Real , Channels > Color< Real , Channels >::operator - ( const Color& clr ) const { Color _clr; for( int c=0 ; c<Channels ; c++ ) _clr[c] = this->c[c]-clr.c[c] ; return _clr; };
template< class Real , int Channels > Color< Real , Channels > Color< Real , Channels >::operator - ( void             ) const { Color _clr; for( int c=0 ; c<Channels ; c++ ) _clr[c] = -this->c[c] ; return _clr; };
template< class Real , int Channels > Color< Real , Channels > Color< Real , Channels >::operator * ( const Real& s    ) const { Color _clr; for( int c=0 ; c<Channels ; c++ ) _clr[c] = this->c[c] * s ; return _clr; };
template< class Real , int Channels > Color< Real , Channels > Color< Real , Channels >::operator / ( const Real& s    ) const { Color _clr; for( int c=0 ; c<Channels ; c++ ) _clr[c] = this->c[c] / s ; return _clr; };
// An essential method for getting the conjugate gradient solver to work
//template< class Real , int Channels > Real Color< Real , Channels >::operator * ( const Color& clr ) const { Real r = (Real)0 ; for( int c=0 ; c<Channels ; c++ ) r += this->c[c]*clr[c] ; return r; };
template< class Real , int Channels > Color< Real , Channels > Color< Real , Channels >::operator * ( const Color& clr ) const { Color _clr ; for( int c=0 ; c<Channels ; c++ ) _clr[c] = this->c[c]*clr[c] ; return _clr; };

template< class Real , int Channels > Color< Real , Channels >& Color< Real , Channels >::operator += ( const Color& clr ){ for( int c=0 ; c<Channels ; c++ ) this->c[c] += clr.c[c] ; return (*this); }
template< class Real , int Channels > Color< Real , Channels >& Color< Real , Channels >::operator -= ( const Color& clr ){ for( int c=0 ; c<Channels ; c++ ) this->c[c] -= clr.c[c] ; return (*this); }
template< class Real , int Channels > Color< Real , Channels >& Color< Real , Channels >::operator *= ( const Real& s    ){ for( int c=0 ; c<Channels ; c++ ) this->c[c] *= s ; return (*this); }
template< class Real , int Channels > Color< Real , Channels >& Color< Real , Channels >::operator /= ( const Real& s    ){ for( int c=0 ; c<Channels ; c++ ) this->c[c] /= s ; return (*this); }
template< class Real , int Channels >
bool Color< Real , Channels >::operator == ( const Color< Real , Channels >& clr ) const
{
	for( int c=0 ; c<Channels ; c++ ) if( this->c[c]!=clr.c[c] ) return false;
	return true;
}
template< class Real , int Channels >
bool Color< Real , Channels >::operator != ( const Color< Real , Channels >& clr ) const
{
	for( int c=0 ; c<Channels ; c++ ) if( this->c[c]!=clr.c[c] ) return true;
	return false;
}
template< class Real , int Channels >
bool Color< Real , Channels >::operator <  ( const Color< Real , Channels >& clr ) const
{
	for( int c=0 ; c<Channels ; c++ )
		if     ( this->c[c]<clr.c[c] ) return true;
		else if( this->c[c]>clr.c[c] ) return false;
	return false;
}
template< class Real , int Channels >
bool Color< Real , Channels >::operator <= ( const Color< Real , Channels >& clr ) const
{
	for( int c=0 ; c<Channels ; c++ )
		if     ( this->c[c]<clr.c[c] ) return true;
		else if( this->c[c]>clr.c[c] ) return false;
	return true;
}
template< class Real , int Channels >
bool Color< Real , Channels >::operator >  ( const Color< Real , Channels >& clr ) const
{
	for( int c=0 ; c<Channels ; c++ )
		if     ( clr.c[c]<this->c[c] ) return true;
		else if( clr.c[c]>this->c[c] ) return false;
	return false;
}
template< class Real , int Channels >
bool Color< Real , Channels >::operator >= ( const Color< Real , Channels >& clr ) const
{
	for( int c=0 ; c<Channels ; c++ )
		if     ( clr.c[c]<this->c[c] ) return true;
		else if( clr.c[c]>this->c[c] ) return false;
	return true;
}
