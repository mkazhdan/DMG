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

#ifndef P_POLYNOMIAL_INCLUDED
#define P_POLYNOMIAL_INCLUDED
#include <vector>
#include "Polynomial.h"

#include "Util/Array.h"

template< int Degree , class Real >
class StartingPolynomial
{
public:
	Polynomial< Degree , Real > p;
	Real start;

	template<int Degree2>
	StartingPolynomial< Degree+Degree2 , Real >  operator * ( const StartingPolynomial< Degree2 , Real >& p ) const;
	StartingPolynomial scale(const Real& s) const;
	StartingPolynomial shift(const Real& t) const;
	int operator < (const StartingPolynomial& sp) const;
	static int Compare(const void* v1,const void* v2);
};

template< int Degree , class Real=double >
class PPolynomial
{
public:
	size_t polyCount;
	Pointer( StartingPolynomial< Degree , Real > ) polys;

	PPolynomial( void );
	PPolynomial( const PPolynomial< Degree , Real >& p );
	~PPolynomial(void);

	PPolynomial& operator = (const PPolynomial& p);

	int size(void) const;

	void set(const size_t& size);
	// Note: this method will sort the elements in sps
	void set( StartingPolynomial< Degree , Real >* sps , const int& count );
	void reset( const size_t& newSize );

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

	Real integralSmoothSine		( const int& smoothness , const Real& tMin , const Real& tMax ) const;
	Real integralSmoothCosine	( const int& smoothness , const Real& tMin , const Real& tMax ) const;
	Real integralSmoothCosecant	( const Real& tMin , const Real& tMax , const int& samples ) const;
	Real integralSmoothSecant	( const Real& tMin , const Real& tMax , const int& samples ) const;
	Real integralSmoothSine		( const int& smoothness , const Real& a , const Real& tMin , const Real& tMax ) const;
	Real integralSmoothCosine	( const int& smoothness , const Real& a , const Real& tMin , const Real& tMax ) const;
	Real integralSmoothCosecant	( const Real& a , const Real& tMin , const Real& tMax , const int& samples ) const;
	Real integralSmoothSecant	( const Real& a , const Real& tMin , const Real& tMax , const int& samples ) const;

	Real Integral			( void ) const;
	Real IntegralSine		( void ) const;
	Real IntegralCosine		( void ) const;
	Real IntegralCosecant	( const int& samples ) const;
	Real IntegralSecant		( const int& samples ) const;

	template< int Degree2 >
	PPolynomial< Degree , Real >& operator = (const PPolynomial< Degree2 , Real >& p);

	template< int Degree2 , class Real2 >
	operator PPolynomial< Degree2 , Real2 > ();

	PPolynomial  operator + (const PPolynomial& p) const;
	PPolynomial  operator - (const PPolynomial& p) const;

	template< int Degree2 >
	PPolynomial< Degree+Degree2 , Real > operator * ( const Polynomial< Degree2 , Real >& p) const;

	template< int Degree2 >
	PPolynomial< Degree+Degree2 , Real > operator * ( const PPolynomial< Degree2 , Real >& p) const;


	PPolynomial& operator += (const Real& s);
	PPolynomial& operator -= (const Real& s);
	PPolynomial& operator *= (const Real& s);
	PPolynomial& operator /= (const Real& s);
	PPolynomial  operator +  (const Real& s) const;
	PPolynomial  operator -  (const Real& s) const;
	PPolynomial  operator *  (const Real& s) const;
	PPolynomial  operator /  (const Real& s) const;

	PPolynomial& addScaled(const PPolynomial& poly,const Real& scale);

	PPolynomial scale(const Real& s) const;
	PPolynomial shift(const Real& t) const;

	PPolynomial< Degree-1 , Real > derivative(void) const;
	PPolynomial< Degree+1 , Real > integral(void) const;

	void getSolutions(const Real& c,std::vector<Real>& roots,const Real& EPS,const Real& min=-DBL_MAX,const Real& max=DBL_MAX) const;

	void printnl(void) const;

	PPolynomial< Degree+1 , Real > MovingAverage(const Real& radius);

	static PPolynomial ConstantFunction(const Real& radius=0.5);
	static PPolynomial GaussianApproximation(const Real& radius=0.5);
	void write(FILE* fp,const int& samples,const Real& min,const Real& max) const;
};
#include "PPolynomial.inl"
#endif // P_POLYNOMIAL_INCLUDED
