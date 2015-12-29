/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

#ifndef POLYNOMIAL_INCLUDED
#define POLYNOMIAL_INCLUDED



//template< int Degree , class Real=double >
template< int Degree , class Real >
class Polynomial
{
public:
	Real coefficients[Degree+1];

	Polynomial( void );
	template< int Degree2 , class Real2 > Polynomial( const Polynomial< Degree2 , Real2 >& P );

	Real operator()(const Real& t) const;
	Real integral			( const Real& tMin , const Real& tMax ) const;
	Real integralSine		( const Real& tMin , const Real& tMax ) const;
	Real integralCosine		( const Real& tMin , const Real& tMax ) const;
	Real integralCosecant	( const Real& tMin , const Real& tMax , const int& samples ) const;
	Real integralSecant		( const Real& tMin , const Real& tMax , const int& samples ) const;
	Real integralSine		( const Real& a , const Real& tMin , const Real& tMax ) const;
	Real integralCosine		( const Real& a , const Real& tMin , const Real& tMax ) const;
	Real integralCosecant	( const Real& a , const Real& tMin , const Real& tMax , const int& samples ) const;
	Real integralSecant		( const Real& a , const Real& tMin , const Real& tMax , const int& samples ) const;
	Real integralCosecant	( const Real& a , const Real& tMin , const Real& tMax , bool noZero=false ) const;

	int operator == ( const Polynomial& p ) const;
	int operator != ( const Polynomial& p ) const;
	int isZero(void) const;
	void setZero(void);

	template< int Degree2 > Polynomial& operator  = ( const Polynomial< Degree2 , Real > &p);
	Polynomial& operator += (const Polynomial& p);
	Polynomial& operator -= (const Polynomial& p);
	Polynomial  operator -  (void) const;
	Polynomial  operator +  (const Polynomial& p) const;
	Polynomial  operator -  (const Polynomial& p) const;
	template< int Degree2 >
	Polynomial< Degree+Degree2 , Real >  operator * ( const Polynomial< Degree2 , Real >& p) const;

	Polynomial& operator += (const Real& s);
	Polynomial& operator -= (const Real& s);
	Polynomial& operator *= (const Real& s);
	Polynomial& operator /= (const Real& s);
	Polynomial  operator +  (const Real& s) const;
	Polynomial  operator -  (const Real& s) const;
	Polynomial  operator *  (const Real& s) const;
	Polynomial  operator /  (const Real& s) const;

	Polynomial scale(const Real& s) const;
	Polynomial shift(const Real& t) const;

	Polynomial< Degree-1 , Real > derivative(void) const;
	Polynomial< Degree+1 , Real > integral(void) const;

	void printnl(void) const;

	Polynomial& addScaled(const Polynomial& p,const Real& scale);

	static void Negate(const Polynomial& in,Polynomial& out);
	static void Subtract(const Polynomial& p1,const Polynomial& p2,Polynomial& q);
	static void Scale(const Polynomial& p,const Real& w,Polynomial& q);
	static void AddScaled(const Polynomial& p1,const Real& w1,const Polynomial& p2,const Real& w2,Polynomial& q);
	static void AddScaled(const Polynomial& p1,const Polynomial& p2,const Real& w2,Polynomial& q);
	static void AddScaled(const Polynomial& p1,const Real& w1,const Polynomial& p2,Polynomial& q);

	void getSolutions( const Real& c , std::vector<Real>& roots , const Real& EPS ) const;

	static Polynomial Cosecant( Real center );
	static Polynomial XCscX( void );
	static Real IntegratePolyCosecant( const Real& a , const Real& tMin , const Real& tMax , const int& polyDegree );
};
template< class Real > Real Bernoulli( int N );

#include "Polynomial.inl"
#endif // POLYNOMIAL_INCLUDED
