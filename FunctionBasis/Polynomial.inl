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

#include <float.h>
#include <math.h>
#include <algorithm>

////////////////
// Polynomial //
////////////////
template< int Degree , class Real >
Polynomial< Degree , Real >::Polynomial( void ){memset(coefficients,0,sizeof(Real)*(Degree+1));}
template< int Degree , class Real >
template< int Degree2 , class Real2 >
Polynomial< Degree , Real >::Polynomial(const Polynomial< Degree2 , Real2 >& P)
{
	memset( coefficients , 0 , sizeof(Real)*(Degree+1));
	for(int i=0;i<=Degree && i<=Degree2;i++) coefficients[i]=Real( P.coefficients[i] );
}


template< int Degree , class Real >
template< int Degree2 >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::operator  = (const Polynomial< Degree2 , Real > &p){
	int d=Degree<Degree2?Degree:Degree2;
	memset(coefficients,0,sizeof(Real)*(Degree+1));
	memcpy(coefficients,p.coefficients,sizeof(Real)*(d+1));
	return *this;
}

template< int Degree , class Real >
Polynomial< Degree-1 , Real > Polynomial< Degree , Real >::derivative(void) const{
	Polynomial< Degree-1 , Real > p;
	for( int i=0 ; i<Degree ; i++ ) p.coefficients[i] = coefficients[i+1]*Real(i+1);
	return p;
}

template< int Degree , class Real >
Polynomial< Degree+1 , Real > Polynomial< Degree , Real >::integral( void ) const
{
	Polynomial< Degree+1 , Real > p;
	p.coefficients[0]=0;
	for( int i=0 ; i<=Degree ; i++ ) p.coefficients[i+1] = coefficients[i] / Real( i+1 );
	return p;
}

template<> inline float  Polynomial< 0 , float  >::integralSine  ( const float&  tMin , const float&  tMax ) const { return coefficients[0] * ( cos( tMin ) - cos( tMax ) ); }
template<> inline double Polynomial< 0 , double >::integralSine  ( const double& tMin , const double& tMax ) const { return coefficients[0] * ( cos( tMin ) - cos( tMax ) ); }
template<> inline float  Polynomial< 0 , float  >::integralCosine( const float&  tMin , const float&  tMax ) const { return coefficients[0] * ( sin( tMax ) - sin( tMin ) ); }
template<> inline double Polynomial< 0 , double >::integralCosine( const double& tMin , const double& tMax ) const { return coefficients[0] * ( sin( tMax ) - sin( tMin ) ); }
template<> inline float  Polynomial< 0 , float  >::integralSine  ( const float&  a , const float&  tMin , const float&  tMax ) const { return - coefficients[0] * ( cos( tMax*a ) - cos( tMin*a ) ) / a; }
template<> inline double Polynomial< 0 , double >::integralSine  ( const double& a , const double& tMin , const double& tMax ) const { return - coefficients[0] * ( cos( tMax*a ) - cos( tMin*a ) ) / a; }
template<> inline float  Polynomial< 0 , float  >::integralCosine( const float&  a , const float&  tMin , const float&  tMax ) const { return   coefficients[0] * ( sin( tMax*a ) - sin( tMin*a ) ) / a; }
template<> inline double Polynomial< 0 , double >::integralCosine( const double& a , const double& tMin , const double& tMax ) const { return   coefficients[0] * ( sin( tMax*a ) - sin( tMin*a ) ) / a; }
template< int Degree , class Real>
Real Polynomial< Degree , Real >::integralSine( const Real& tMin,const Real& tMax ) const
{
	return -( (*this)( tMax ) * cos( tMax ) - (*this)( tMin ) * cos( tMin ) ) + derivative().integralCosine( tMin , tMax );
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralCosine( const Real& tMin,const Real& tMax ) const
{
	return  ( (*this)( tMax ) * sin( tMax ) - (*this)( tMin ) * sin( tMin ) ) - derivative().integralSine( tMin , tMax );
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralSine( const Real& a , const Real& tMin,const Real& tMax ) const
{
	return -( (*this)( tMax ) * cos( tMax*a ) - (*this)( tMin ) * cos( tMin*a ) ) / a + derivative().integralCosine( a , tMin , tMax ) / a;
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralCosine( const Real& a , const Real& tMin,const Real& tMax ) const
{
	return  ( (*this)( tMax ) * sin( tMax*a ) - (*this)( tMin ) * sin( tMin*a ) ) / a - derivative().integralSine  ( a , tMin , tMax ) / a;
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralCosecant( const Real& tMin , const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for(int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real(samples) * Real( i+0.5 );
		sum += (*this)( t ) / sin( t );
	}
	return sum / Real( samples ) * (tMax-tMin);
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralSecant( const Real& tMin,const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for(int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real(samples) * Real( i+0.5 );
		sum += (*this)( t ) / cos( t );
	}
	return sum / samples * (tMax-tMin);
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralCosecant( const Real& a , const Real& tMin,const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for(int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real(samples) * Real( i+0.5 );
		sum += (*this)( t ) / sin( t*a );
	}
	return sum / Real( samples ) * (tMax-tMin);
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralCosecant( const Real& a , const Real& tMin , const Real& tMax , bool noZero ) const
{
	Real sum = 0;
	if( noZero ) for( int d=1 ; d<=Degree ; d++ ) sum += coefficients[d] * Polynomial< 10 , Real >::IntegratePolyCosecant( a , tMin , tMax , d );
	else         for( int d=0 ; d<=Degree ; d++ ) sum += coefficients[d] * Polynomial< 10 , Real >::IntegratePolyCosecant( a , tMin , tMax , d );
	return sum;
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integralSecant( const Real& a , const Real& tMin,const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for(int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / (samples-1) * i;
		sum += (*this)( t ) / cos( t*a );
	}
	return sum / samples * (tMax-tMin);
}

template< int Degree , class Real >
Real Polynomial< Degree , Real >::operator() (const Real& t) const{
	Real temp=1;
	Real v=0;
	for(int i=0;i<=Degree;i++){
		v+=temp*coefficients[i];
		temp*=t;
	}
	return v;
}
template< int Degree , class Real >
Real Polynomial< Degree , Real >::integral( const Real& tMin , const Real& tMax ) const
{
	Real v=0;
	Real t1,t2;
	t1=tMin;
	t2=tMax;
	for( int i=0 ; i<=Degree ; i++ )
	{
		v += coefficients[i]*(t2-t1)/Real(i+1);
		if( t1!=-Real(DBL_MAX) && t1!=Real(DBL_MAX) ) t1*=tMin;
		if( t2!=-Real(DBL_MAX) && t2!=Real(DBL_MAX) ) t2*=tMax;
	}
	return v;
}
template< int Degree , class Real >
int Polynomial< Degree , Real >::operator == (const Polynomial& p) const{
	for(int i=0;i<=Degree;i++){if(coefficients[i]!=p.coefficients[i]){return 0;}}
	return 1;
}
template< int Degree , class Real >
int Polynomial< Degree , Real >::operator != (const Polynomial& p) const{
	for(int i=0;i<=Degree;i++){if(coefficients[i]==p.coefficients[i]){return 0;}}
	return 1;
}
template< int Degree , class Real >
int Polynomial< Degree , Real >::isZero(void) const{
	for(int i=0;i<=Degree;i++){if(coefficients[i]!=0){return 0;}}
	return 1;
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::setZero(void){memset(coefficients,0,sizeof(Real)*(Degree+1));}

template< int Degree , class Real >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::addScaled(const Polynomial& p,const Real& s){
	for(int i=0;i<=Degree;i++){coefficients[i]+=p.coefficients[i]*s;}
	return *this;
}
template< int Degree , class Real >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::operator += (const Polynomial< Degree , Real >& p){
	for(int i=0;i<=Degree;i++){coefficients[i]+=p.coefficients[i];}
	return *this;
}
template< int Degree , class Real >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::operator -= (const Polynomial< Degree , Real >& p){
	for(int i=0;i<=Degree;i++){coefficients[i]-=p.coefficients[i];}
	return *this;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::operator + (const Polynomial< Degree , Real >& p) const{
	Polynomial q;
	for(int i=0;i<=Degree;i++){q.coefficients[i]=(coefficients[i]+p.coefficients[i]);}
	return q;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::operator - (const Polynomial< Degree , Real >& p) const{
	Polynomial q;
	for(int i=0;i<=Degree;i++)	{q.coefficients[i]=coefficients[i]-p.coefficients[i];}
	return q;
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::Scale(const Polynomial& p,const Real& w,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p.coefficients[i]*w;}
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::AddScaled(const Polynomial& p1,const Real& w1,const Polynomial& p2,const Real& w2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]*w1+p2.coefficients[i]*w2;}
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::AddScaled(const Polynomial& p1,const Real& w1,const Polynomial& p2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]*w1+p2.coefficients[i];}
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::AddScaled(const Polynomial& p1,const Polynomial& p2,const Real& w2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]+p2.coefficients[i]*w2;}
}

template< int Degree , class Real >
void Polynomial< Degree , Real >::Subtract(const Polynomial &p1,const Polynomial& p2,Polynomial& q){
	for(int i=0;i<=Degree;i++){q.coefficients[i]=p1.coefficients[i]-p2.coefficients[i];}
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::Negate(const Polynomial& in,Polynomial& out){
	out=in;
	for(int i=0;i<=Degree;i++){out.coefficients[i]=-out.coefficients[i];}
}

template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::operator - (void) const{
	Polynomial q=*this;
	for(int i=0;i<=Degree;i++){q.coefficients[i]=-q.coefficients[i];}
	return q;
}
template< int Degree , class Real >
template< int Degree2 >
Polynomial< Degree+Degree2 , Real > Polynomial< Degree , Real >::operator * (const Polynomial< Degree2 , Real >& p) const
{
	Polynomial< Degree+Degree2 , Real > q;
	for(int i=0;i<=Degree;i++){for(int j=0;j<=Degree2;j++){q.coefficients[i+j]+=coefficients[i]*p.coefficients[j];}}
	return q;
}

template< int Degree , class Real >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::operator += (const Real& s){
	coefficients[0]+=s;
	return *this;
}
template< int Degree , class Real >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::operator -= (const Real& s){
	coefficients[0]-=s;
	return *this;
}
template< int Degree , class Real >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::operator *= (const Real& s){
	for(int i=0;i<=Degree;i++){coefficients[i]*=s;}
	return *this;
}
template< int Degree , class Real >
Polynomial< Degree , Real >& Polynomial< Degree , Real >::operator /= (const Real& s){
	for(int i=0;i<=Degree;i++){coefficients[i]/=s;}
	return *this;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::operator + (const Real& s) const{
	Polynomial< Degree , Real > q=*this;
	q.coefficients[0]+=s;
	return q;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::operator - (const Real& s) const{
	Polynomial q=*this;
	q.coefficients[0]-=s;
	return q;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::operator * (const Real& s) const{
	Polynomial q;
	for(int i=0;i<=Degree;i++){q.coefficients[i]=coefficients[i]*s;}
	return q;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::operator / (const Real& s) const{
	Polynomial q;
	for(int i=0;i<=Degree;i++){q.coefficients[i]=coefficients[i]/s;}
	return q;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::scale(const Real& s) const{
	Polynomial q=*this;
	Real s2=1.0;
	for(int i=0;i<=Degree;i++){
		q.coefficients[i]*=s2;
		s2/=s;
	}
	return q;
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::shift(const Real& t) const{
	Polynomial< Degree , Real > q;
	for(int i=0;i<=Degree;i++){
		Real temp=1;
		for(int j=i;j>=0;j--){
			q.coefficients[j] += coefficients[i]*temp;
			temp *= -t*Real(j);
			temp /= Real(i-j+1);
		}
	}
	return q;
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::printnl(void) const{
	for(int j=0;j<=Degree;j++)
	{
		printf( "%6.4f x^%d " , double( coefficients[j] ) , j );
		if( j<Degree && coefficients[j+1]>=Real(0) ) printf( "+" );
	}
	printf("\n");
}
template< int Degree , class Real >
void Polynomial< Degree , Real >::getSolutions(const Real& c,std::vector<Real>& roots,const Real& EPS) const {
	Real r[4][2];
	int rCount=0;
	roots.clear();
	switch(Degree){
	case 1:
		rCount=Factor(coefficients[1],coefficients[0]-c,r,EPS);
		break;
	case 2:
		rCount=Factor(coefficients[2],coefficients[1],coefficients[0]-c,r,EPS);
		break;
	case 3:
		rCount=Factor(coefficients[3],coefficients[2],coefficients[1],coefficients[0]-c,r,EPS);
		break;
//	case 4:
//		rCount=Factor(coefficients[4],coefficients[3],coefficients[2],coefficients[1],coefficients[0]-c,r,EPS);
//		break;
	default:
		printf("Can't solve polynomial of degree: %d\n",Degree);
	}
	for(int i=0;i<rCount;i++){
		if(fabs(r[i][1])<=EPS){
			roots.push_back(r[i][0]);
//printf("%d] %f\t%f\n",i,r[i][0],(*this)(r[i][0])-c);
		}
	}
}
template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::Cosecant( Real center )
{
	Polynomial< 5 , Real > p;
	Real c = cos( center );
	Real s = Real( 1. ) / sin( center );
	Real s2 = s*s;
	Real cs = c*s;

	p.coefficients[0] =  s;
	p.coefficients[1] = -cs;
	p.coefficients[2] =  s2 / Real( 2. );
	p.coefficients[3] =  - Real(  2. ) * cs * s2 / Real( 3. * 2. );
	p.coefficients[4] = (  Real(  6. ) * s2 - Real( 4. ) ) * s2 / Real( 4. * 3. * 2. );
	p.coefficients[5] = ( -Real( 24. ) * s2 + Real( 8. ) ) * s2 * cs / Real( 5. * 4. * 3. * 2. );
	return Polynomial( p ).shift( center );
}
template<> inline Polynomial< 10 , float > Polynomial< 10 , float >::XCscX( void )
{
	Polynomial< 10 , float > p;
	p.coefficients[ 0] =   1.0f;
	p.coefficients[ 2] =   1.0f       /  6.f;
	p.coefficients[ 4] =   7.0f       / 30.f /                                      4.f / 3.f;
	p.coefficients[ 6] =  31.0f       / 42.f /                          6.f / 5.f / 4.f / 3.f;
	p.coefficients[ 8] = 127.0f       / 30.f /              8.f / 7.f / 6.f / 5.f / 4.f / 3.f;
	p.coefficients[10] = 511.0f * 5.f / 66.f / 10.f / 9.f / 8.f / 7.f / 6.f / 5.f / 4.f / 3.f;
/*
csc(x) = 1/x + \sum_{k=1}^{\infty} (-1)^{k-1}*2*(2^{2k-1}-1)B_{2k} x^{2k-1} / (2k)!
B_0 =  1
B_2 =  1/6
B_4 = -1/30
B_6 =  1/42
B_8 = -1/30
B_10 = 5/66
*/
	return p;
}
template<> inline Polynomial< 10 , double > Polynomial< 10 , double >::XCscX( void )
{
	Polynomial< 10 , double > p;
	p.coefficients[ 0] =   1.0;
	p.coefficients[ 2] =   1.0     / 6;
	p.coefficients[ 4] =   7.0     / 30 /                          4 / 3;
	p.coefficients[ 6] =  31.0     / 42 /                  6 / 5 / 4 / 3;
	p.coefficients[ 8] = 127.0     / 30 /          8 / 7 / 6 / 5 / 4 / 3;
	p.coefficients[10] = 511.0 * 5 / 66 / 10 / 9 / 8 / 7 / 6 / 5 / 4 / 3;
	return p;
}

template< int Degree , class Real >
Real Polynomial< Degree , Real >::IntegratePolyCosecant( const Real& a , const Real& tMin , const Real& tMax , const int& polyDegree )
{
	if( polyDegree==0 ) return - ( log( fabs( Real( 1. ) / sin( tMax*a ) + cos( tMax*a ) / sin( tMax*a ) ) ) - log( fabs( Real( 1. ) / sin( tMin*a ) + cos( tMin*a ) / sin( tMin*a ) ) ) ) / a;
	else
	{
		Real integral = 0;
		Real coefficients[ Degree+1 ];
		coefficients[0] = Real( 1. ) / Real( polyDegree );
		Real sign = 1.;
		for( int k=1 ; k<=Degree ; k++ )
		{
			coefficients[k] = sign * Real( 2. ) * Real( ( 1<<(2*k-1) ) - 1. ) * Bernoulli< Real >( 2*k ) / Real( 2*k + polyDegree );
			sign *= -1.;
		}
		Real _tMin = 1. , _tMax = 1.;
		_tMin /= a , _tMax /= a;
		for( int k=0 ; k<polyDegree ; k++ ) _tMin *= tMin , _tMax *= tMax;
		Real tMin2 = tMin * tMin * a * a , tMax2 = tMax * tMax * a * a;
#if 1
		Real integralMin = 0. , integralMax = 0.;
		for( int k=Degree ; k>0 ; k-- )
		{
			Real fact = 1.;
			if( k ) fact = Real( ( 2*k ) * ( 2*k-1 ) );
			integralMin = ( integralMin + coefficients[k] ) * tMin2 / fact;
			integralMax = ( integralMax + coefficients[k] ) * tMax2 / fact;
		}
		integralMin += coefficients[0];
		integralMax += coefficients[0];
		integralMin *= _tMin;
		integralMax *= _tMax;
		integral = integralMax - integralMin;
#else
		for( int k=0 ; k<=Degree ; k++ )
		{
			integral += coefficients[k] * ( _tMax - _tMin );
			_tMin *= tMin2 / Real( (2*k+2)*(2*k+1) );
			_tMax *= tMax2 / Real( (2*k+2)*(2*k+1) );
		}
#endif
		return integral;
//csc(x) = 1/x + \sum_{k=1}^{\infty} (-1)^{k-1}*2*(2^{2k-1}-1)B_{2k} x^{2k-1} / (2k)!
	}
}

template< int Degree , class Real >
Polynomial< Degree , Real > Polynomial< Degree , Real >::XCscX( void )
{
	return Polynomial< Degree , Real >( Polynomial< 10 , Real >::XCscX() );
}
template< class Real >
Real Bernoulli( int N )
{
	Real* table = new Real[N+1];
	for( int i=0 ; i<=N ; i++ )
	{
		table[i] = Real(1.) /  Real( i+1 );
		for( int j=i ; j>0 ; j-- ) table[j-1] = Real( j ) * ( table[j-1]-table[j] );
	}
	Real ret = table[0];
	delete[] table;
	return ret;
}