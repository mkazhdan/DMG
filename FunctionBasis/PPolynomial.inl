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


////////////////////////
// StartingPolynomial //
////////////////////////
template< int Degree , class Real >
template<int Degree2>
StartingPolynomial< Degree+Degree2 , Real > StartingPolynomial< Degree , Real >::operator * ( const StartingPolynomial< Degree2 , Real >& p) const
{
	StartingPolynomial< Degree+Degree2 , Real > sp;
	if(start>p.start) sp.start=start;
	else              sp.start=p.start;
	sp.p=this->p*p.p;
	return sp;
}
template< int Degree , class Real >
StartingPolynomial< Degree , Real > StartingPolynomial< Degree , Real >::scale(const Real& s) const{
	StartingPolynomial q;
	q.start=start*s;
	q.p=p.scale(s);
	return q;
}
template< int Degree , class Real >
StartingPolynomial< Degree , Real > StartingPolynomial< Degree , Real >::shift(const Real& s) const{
	StartingPolynomial q;
	q.start=start+s;
	q.p=p.shift(s);
	return q;
}

template< int Degree , class Real >
int StartingPolynomial< Degree , Real >::operator < (const StartingPolynomial< Degree , Real >& sp) const{
	if(start<sp.start){return 1;}
	else{return 0;}
}
template< int Degree , class Real >
int StartingPolynomial< Degree , Real >::Compare(const void* v1,const void* v2){
	Real d=((StartingPolynomial*)(v1))->start-((StartingPolynomial*)(v2))->start;
	if		( d<Real(0) )	return -1;
	else if	( d>Real(0) )	return  1;
	else					return  0;
}

/////////////////
// PPolynomial //
/////////////////
template< int Degree , class Real >
PPolynomial< Degree , Real >::PPolynomial(void)
{
	polyCount = 0;
	polys = NullPointer< StartingPolynomial< Degree , Real > >( );
}
template< int Degree , class Real >
PPolynomial< Degree , Real >::PPolynomial( const PPolynomial< Degree , Real >& p )
{
	polyCount=0;
	polys = NullPointer< StartingPolynomial< Degree , Real > >( );
	set( p.polyCount );
	memcpy( polys , p.polys , sizeof(StartingPolynomial< Degree , Real >)*p.polyCount);
}

template< int Degree , class Real >
PPolynomial< Degree , Real >::~PPolynomial( void )
{
	FreePointer( polys );
	polyCount = 0;
}
template< int Degree , class Real >
void PPolynomial< Degree , Real >::set( StartingPolynomial< Degree , Real >* sps , const int& count )
{
	int i,c=0;
	set( count );
	qsort( sps , count , sizeof( StartingPolynomial< Degree , Real > ) , StartingPolynomial< Degree , Real >::Compare );
	for( i=0 ; i<count ; i++ )
		if(!c || sps[i].start!=polys[c-1].start)	polys[c++]    = sps[i];
		else										polys[c-1].p += sps[i].p;
	reset( c );
}
template< int Degree , class Real >
int PPolynomial< Degree , Real >::size(void) const{return (int)(sizeof(StartingPolynomial< Degree , Real >)*polyCount);}

template< int Degree , class Real >
void PPolynomial< Degree , Real >::set(const size_t &size)
{
	FreePointer( polys );
	polyCount = 0;
	polyCount = size;
	if( size )
	{
		polys = AllocPointer< StartingPolynomial< Degree , Real > >( (int)(size) );
		memset( polys , 0 , sizeof( StartingPolynomial< Degree , Real > ) * size );
	}
}
template< int Degree , class Real >
void PPolynomial< Degree , Real >::reset( const size_t& newSize )
{
	Pointer( StartingPolynomial< Degree , Real > ) temp = AllocPointer< StartingPolynomial< Degree , Real > >( (int)( newSize ) );
	memcpy( temp , polys , sizeof( StartingPolynomial< Degree , Real > ) * newSize );
	FreePointer( polys );
	polys = temp;
	polyCount = newSize;
}

template< int Degree , class Real >
PPolynomial< Degree , Real >& PPolynomial< Degree , Real >::operator = ( const PPolynomial< Degree , Real >& p )
{
	set( p.polyCount );
	memcpy( polys , p.polys , sizeof( StartingPolynomial< Degree , Real > ) * p.polyCount );
	return *this;
}

template< int Degree , class Real >
template< int Degree2 >
PPolynomial< Degree , Real >& PPolynomial< Degree , Real >::operator  = (const PPolynomial< Degree2 , Real >& p){
	set( p.polyCount );
	for( int i=0 ; i<(int)(polyCount) ; i++ )
	{
		polys[i].start = p.polys[i].start;
		polys[i].p = p.polys[i].p;
	}
	return *this;
}
template< int Degree , class Real >
template< int Degree2 , class Real2 >
PPolynomial< Degree , Real >::operator PPolynomial< Degree2 , Real2 > ()
{
	PPolynomial< Degree2 , Real2 > poly;
	poly.set(polyCount);

	for(int i=0;i<(int)(polyCount);i++){
		poly.polys[i].start=polys[i].start;
		poly.polys[i].p=polys[i].p;
	}
	return poly;
}

template< int Degree , class Real >
Real PPolynomial< Degree , Real >::operator ()(const Real& t) const{
	Real v=0;
	for(int i=0;i<(int)(polyCount) && t>polys[i].start;i++){v+=polys[i].p(t);}
	return v;
}

template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integral( const Real& tMin , const Real& tMax ) const
{
	int m=1;
	Real start,end,s,v=0;
	start=tMin;
	end=tMax;
	if( tMin > tMax )
	{
		m = -1;
		start = tMax;
		end = tMin;
	}
	for( int i=0 ; i<(int)(polyCount) && polys[i].start<end ; i++ )
	{
		if( start < polys[i].start )	s=polys[i].start;
		else							s=start;
		v += polys[i].p.integral( s , end );
	}
	return v * Real( m );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSine( const Real& tMin , const Real& tMax ) const
{
	int m=1;
	Real start,end,s,v=0;
	start=tMin;
	end=tMax;
	if( tMin > tMax )
	{
		m=-1;
		start=tMax;
		end=tMin;
	}
	for(int i=0;i<(int)(polyCount) && polys[i].start<end;i++)
	{
		if( start < polys[i].start )	s=polys[i].start;
		else							s=start;
		v+=polys[i].p.integralSine( s , end );
	}
	return v * Real( m );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSine( const Real& a , const Real& tMin , const Real& tMax ) const
{
	int m = 1;
	Real start , end , s , v = 0;
	start = tMin;
	end   = tMax;
	if( tMin > tMax )
	{
		m = -1;
		start = tMax;
		end   = tMin;
	}
	for( int i=0 ; i<(int)(polyCount) && polys[i].start<end ; i++ )
	{
		if( start < polys[i].start )	s = polys[i].start;
		else							s = start;
		v += polys[i].p.integralSine( a , s , end );
	}
	return v * Real( m );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralCosine( const Real& tMin , const Real& tMax ) const
{
	int m=1;
	Real start,end,s,v=0;
	start=tMin;
	end=tMax;
	if( tMin > tMax )
	{
		m=-1;
		start=tMax;
		end=tMin;
	}
	for(int i=0;i<(int)(polyCount) && polys[i].start<end;i++)
	{
		if( start < polys[i].start )	s=polys[i].start;
		else							s=start;
		v+=polys[i].p.integralCosine( s , end );
	}
	return v * Real( m );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralCosine( const Real& a , const Real& tMin , const Real& tMax ) const
{
	int m=1;
	Real start,end,s,v=0;
	start=tMin;
	end=tMax;
	if( tMin > tMax )
	{
		m=-1;
		start=tMax;
		end=tMin;
	}
	for(int i=0;i<(int)(polyCount) && polys[i].start<end;i++)
	{
		if( start < polys[i].start )	s=polys[i].start;
		else							s=start;
		v+=polys[i].p.integralCosine( a , s , end );
	}
	return v * Real( m );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralCosecant( const Real& tMin , const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for( int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real( samples ) * Real(i+0.5);
		sum += (*this)( t ) / sin( t );
	}
	return sum / Real( samples ) * (tMax-tMin);
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralCosecant( const Real& a , const Real& tMin , const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for(int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real( samples ) * Real(i+0.5);
		sum += (*this)( t ) / sin( t*a );
	}
	return sum / Real( samples ) * (tMax-tMin);
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralCosecant( const Real& a , const Real& tMin , const Real& tMax , bool noZero ) const
{
	int m=1;
	Real start , end , s , v=0;
	start = tMin;
	end   = tMax;
	if( tMin > tMax )
	{
		m=-1;
		start = tMax;
		end   = tMin;
	}
	for( int i=0 ; i<(int)(polyCount) && polys[i].start<end ; i++ )
	{
		if( start < polys[i].start )	s = polys[i].start;
		else							s = start;
		v += polys[i].p.integralCosecant( a , s , end , noZero );
	}
	return v * Real( m );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSecant( const Real& tMin , const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for(int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real(samples-1) * Real(i+0.5);
		sum += (*this)( t ) / cos( t );
	}
	return sum / Real( samples ) * (tMax-tMin);
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSecant( const Real& a , const Real& tMin , const Real& tMax , const int& samples ) const
{
	Real sum = 0;
	for(int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real(samples-1) * Real(i+0.5);
		sum += (*this)( t ) / cos( t*a );
	}
	return sum / Real( samples ) * (tMax-tMin);
}

// Integrators for smooth functions
template<> inline float  PPolynomial< 0 , float  >::integralSmoothSine  ( const int& , const float&  tMin , const float&  tMax ) const { return integralSine( tMin , tMax ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothSine  ( const int& , const double& tMin , const double& tMax ) const { return integralSine( tMin , tMax ); }
template<> inline float  PPolynomial< 0 , float  >::integralSmoothCosine( const int& , const float&  tMin , const float&  tMax ) const { return integralCosine( tMin , tMax ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothCosine( const int& , const double& tMin , const double& tMax ) const { return integralCosine( tMin , tMax ); }
template<> inline float  PPolynomial< 0 , float  >::integralSmoothCosecant( const float&  tMin , const float&  tMax , const int& samples ) const { return integralCosecant( tMin , tMax , samples ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothCosecant( const double& tMin , const double& tMax , const int& samples ) const { return integralCosecant( tMin , tMax , samples ); }
template<> inline float  PPolynomial< 0 , float  >::integralSmoothSecant( const float&  tMin , const float&  tMax , const int& samples ) const { return integralSecant( tMin , tMax , samples ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothSecant( const double& tMin , const double& tMax , const int& samples ) const { return integralSecant( tMin , tMax , samples ); }
template<> inline float  PPolynomial< 0 , float  >::integralSmoothSine  ( const int& , const float&  a , const float&  tMin , const float&  tMax ) const { return integralSine( a , tMin , tMax ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothSine  ( const int& , const double& a , const double& tMin , const double& tMax ) const { return integralSine( a , tMin , tMax ); }
template<> inline float  PPolynomial< 0 , float  >::integralSmoothCosine( const int& , const float&  a , const float&  tMin , const float&  tMax ) const { return integralCosine( a , tMin , tMax ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothCosine( const int& , const double& a , const double& tMin , const double& tMax ) const { return integralCosine( a , tMin , tMax ); }
template<> inline float  PPolynomial< 0 , float  >::integralSmoothCosecant( const float&  a , const float&  tMin , const float&  tMax , const int& samples ) const { return integralCosecant( a , tMin , tMax , samples ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothCosecant( const double& a , const double& tMin , const double& tMax , const int& samples ) const { return integralCosecant( a , tMin , tMax , samples ); }
template<> inline float  PPolynomial< 0 , float  >::integralSmoothSecant( const float&  a , const float&  tMin , const float&  tMax , const int& samples ) const { return integralSecant( a , tMin , tMax , samples ); }
template<> inline double PPolynomial< 0 , double >::integralSmoothSecant( const double& a , const double& tMin , const double& tMax , const int& samples ) const { return integralSecant( a , tMin , tMax , samples ); }
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSmoothSine( const int& smoothness , const Real& tMin , const Real& tMax ) const
{
	if( smoothness>0 ) return -( (*this)( tMax ) * cos( tMax ) - (*this)( tMin ) * cos( tMin ) ) + derivative().integralSmoothCosine( smoothness-1 , tMin , tMax );
	else return integralSine( tMin , tMax );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSmoothSine( const int& smoothness , const Real& a , const Real& tMin , const Real& tMax ) const
{
	if( smoothness>0 ) return -( (*this)( tMax ) * cos( tMax*a ) - (*this)( tMin ) * cos( tMin*a ) ) / a + derivative().integralSmoothCosine( smoothness-1 , a , tMin , tMax ) / a;
	else return integralSine( a , tMin , tMax );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSmoothCosine( const int& smoothness , const Real& tMin , const Real& tMax ) const
{
	if( smoothness>0 ) 	return  ( (*this)( tMax ) * sin( tMax ) - (*this)( tMin ) * sin( tMin ) ) - derivative().integralSmoothSine( smoothness-1 , tMin , tMax );
	else return integralCosine( tMin , tMax );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSmoothCosine( const int& smoothness , const Real& a , const Real& tMin , const Real& tMax ) const
{
	if( smoothness>0 ) return ( (*this)( tMax ) * sin( tMax*a ) - (*this)( tMin ) * sin( tMin*a ) ) / a - derivative().integralSmoothSine( smoothness-1 , a , tMin , tMax ) / a;
	else return integralCosine( a , tMin , tMax );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSmoothCosecant( const Real& tMin , const Real& tMax , const int& samples ) const
{
	// \int 1/sin(x) = \log( \sin( x /2 ) ) - \log( \cos( x / 2 ) )
	return
		(*this)( tMax ) * ( log( sin( tMax/2 ) ) - log( cos( tMax/2 ) ) ) - (*this)( tMin ) * ( log( sin( tMin/2 ) ) - log( cos( tMin/2 ) ) )
		- derivative().integralCosecant( tMin , tMax , samples );
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::integralSmoothCosecant( const Real& a , const Real& tMin , const Real& tMax , const int& samples ) const
{
	// \int 1/sin(x) = \log( \sin( x /2 ) ) - \log( \cos( x / 2 ) )
	// Assuming that if tMin = 0 or tMax = 0 then the value of the polynomial and its derivative at zero are zero.
	PPolynomial< Degree-1 , Real > dPoly = derivative( );
	Real start , end;
	Real a2 = a / Real( 2 );
	if( tMin!=Real(0) ) start = (*this)( tMin ) * log( tan( tMin*a2 ) );
	else                start = dPoly.derivative()( 0 );
	if( tMax!=Real(0) ) end   = (*this)( tMax ) * log( tan( tMax*a2 ) );
	else                end   = dPoly.derivative()( 0 );
	Real integral = ( end - start ) / a;
	Real sum = 0;
	for( int i=0 ; i<samples ; i++ )
	{
		Real t = tMin + (tMax-tMin) / Real( samples ) * Real(i+0.5);
		sum += dPoly( t ) * log( tan( t*a2 ) );
	}
	integral -= sum * (tMax-tMin) / Real( samples ) / a;
	return integral;
}
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::Integral		( void ) const { return integral		( polys[0].start , polys[polyCount-1].start ); }
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::IntegralSine	( void ) const { return integralSine	( polys[0].start , polys[polyCount-1].start ); }
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::IntegralCosine	( void ) const { return integralCosine	( polys[0].start , polys[polyCount-1].start ); }
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::IntegralCosecant( const int& samples ) const { return integralCosecant	( polys[0].start , polys[polyCount-1].start , samples ); }
template< int Degree , class Real >
Real PPolynomial< Degree , Real >::IntegralSecant	( const int& samples ) const { return integralSecant	( polys[0].start , polys[polyCount-1].start , samples ); }


template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::operator + (const PPolynomial< Degree , Real >& p) const{
	PPolynomial q;
	int i,j;
	size_t idx=0;
	q.set(polyCount+p.polyCount);
	i=j=-1;

	while(idx<q.polyCount){
		if		(j>=(int)(p.polyCount)-1)				{q.polys[idx]=  polys[++i];}
		else if	(i>=(int)(  polyCount)-1)				{q.polys[idx]=p.polys[++j];}
		else if(polys[i+1].start<p.polys[j+1].start){q.polys[idx]=  polys[++i];}
		else										{q.polys[idx]=p.polys[++j];}
//		if(idx && polys[idx].start==polys[idx-1].start)	{polys[idx-1].p+=polys[idx].p;}
//		else{idx++;}
		idx++;
	}
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::operator - (const PPolynomial< Degree , Real >& p) const{
	PPolynomial q;
	int i,j;
	size_t idx=0;
	q.set(polyCount+p.polyCount);
	i=j=-1;

	while(idx<q.polyCount){
		if		(j>=(int)(p.polyCount)-1)				{q.polys[idx]=  polys[++i];}
		else if	(i>=(int)(  polyCount)-1)				{q.polys[idx].start=p.polys[++j].start;q.polys[idx].p=p.polys[j].p*(-1.0);}
		else if(polys[i+1].start<p.polys[j+1].start){q.polys[idx]=  polys[++i];}
		else										{q.polys[idx].start=p.polys[++j].start;q.polys[idx].p=p.polys[j].p*(-1.0);}
//		if(idx && polys[idx].start==polys[idx-1].start)	{polys[idx-1].p+=polys[idx].p;}
//		else{idx++;}
		idx++;
	}
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real >& PPolynomial< Degree , Real >::addScaled(const PPolynomial< Degree , Real >& p,const Real& scale){
	int i,j;
	StartingPolynomial< Degree , Real >* oldPolys=polys;
	size_t idx=0,cnt=0,oldPolyCount=polyCount;
	polyCount=0;
	polys=NULL;
	set(oldPolyCount+p.polyCount);
	i=j=-1;
	while(cnt<polyCount){
		if		(j>=(int)( p.polyCount)-1)				{polys[idx]=oldPolys[++i];}
		else if	(i>=(int)(oldPolyCount)-1)				{polys[idx].start= p.polys[++j].start;polys[idx].p=p.polys[j].p*scale;}
		else if	(oldPolys[i+1].start<p.polys[j+1].start){polys[idx]=oldPolys[++i];}
		else											{polys[idx].start= p.polys[++j].start;polys[idx].p=p.polys[j].p*scale;}
		if(idx && polys[idx].start==polys[idx-1].start)	{polys[idx-1].p+=polys[idx].p;}
		else{idx++;}
		cnt++;
	}
	free(oldPolys);
	reset(idx);
	return *this;
}
template< int Degree , class Real >
template< int Degree2 >
PPolynomial<Degree+Degree2 , Real > PPolynomial< Degree , Real >::operator * (const PPolynomial< Degree2 , Real >& p) const
{
	PPolynomial< Degree+Degree2 , Real > q;
	StartingPolynomial< Degree+Degree2 , Real > *sp;
	int i,j,spCount=(int)(polyCount*p.polyCount);

	sp = (StartingPolynomial< Degree+Degree2 , Real >* )malloc( sizeof( StartingPolynomial< Degree+Degree2 , Real > ) * spCount );
	for(i=0;i<(int)(polyCount);i++)
		for(j=0;j<(int)(p.polyCount);j++)
			sp[i*p.polyCount+j]=polys[i]*p.polys[j];
	q.set( sp , spCount );
	free( sp );
	return q;
}
template< int Degree , class Real >
template<int Degree2>
PPolynomial< Degree+Degree2 , Real > PPolynomial< Degree , Real >::operator * (const Polynomial< Degree2 , Real >& p) const{
	PPolynomial< Degree+Degree2 , Real > q;
	q.set(polyCount);
	for(int i=0;i<(int)(polyCount);i++)
	{
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p*p;
	}
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::scale(const Real& s) const{
	PPolynomial q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){q.polys[i]=polys[i].scale(s);}
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::shift(const Real& s) const{
	PPolynomial q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){q.polys[i]=polys[i].shift(s);}
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree-1 , Real > PPolynomial< Degree , Real >::derivative(void) const{
	PPolynomial<Degree-1 , Real > q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p.derivative();
	}
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree+1 , Real > PPolynomial< Degree , Real >::integral(void) const{
	int i;
	PPolynomial< Degree+1 , Real > q;
	q.set(polyCount);
	for(i=0;i<(int)(polyCount);i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p.integral();
		q.polys[i].p-=q.polys[i].p(q.polys[i].start);
	}
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real >& PPolynomial< Degree , Real >::operator  += (const Real &s){polys[0].p+=s;}
template< int Degree , class Real >
PPolynomial< Degree , Real >& PPolynomial< Degree , Real >::operator  -= (const Real &s){polys[0].p-=s;}
template< int Degree , class Real >
PPolynomial< Degree , Real >& PPolynomial< Degree , Real >::operator  *= (const Real &s){
	for(int i=0;i<(int)(polyCount);i++){polys[i].p*=s;}
	return *this;
}
template< int Degree , class Real >
PPolynomial< Degree , Real >& PPolynomial< Degree , Real >::operator  /= (const Real &s){
	for(size_t i=0;i<polyCount;i++){polys[i].p/=s;}
	return *this;
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::operator + (const Real& s) const{
	PPolynomial q=*this;
	q+=s;
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::operator - (const Real& s) const{
	PPolynomial q=*this;
	q-=s;
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::operator * (const Real& s) const{
	PPolynomial q=*this;
	q*=s;
	return q;
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::operator / (const Real& s) const{
	PPolynomial q=*this;
	q/=s;
	return q;
}

template< int Degree , class Real >
void PPolynomial< Degree , Real >::printnl(void) const{
	Polynomial< Degree , Real > p;

	if(!polyCount)
	{
		Polynomial< Degree , Real > p;
		printf("[-Infinity,Infinity]\n");
	}
	else{
		for(size_t i=0;i<polyCount;i++)
		{
			if		( polys[i  ].start==Real( DBL_MAX) ) printf( "[Infinity," );
			else if	( polys[i  ].start==Real(-DBL_MAX) ) printf( "[-Infinity," );
			else                                         printf( "[%f," , double( polys[i].start ) );
			if( i+1==polyCount )                         printf( "Infinity]\t" );
			else if ( polys[i+1].start==Real( DBL_MAX) ) printf( "Infinity]\t" );
			else if	( polys[i+1].start==Real(-DBL_MAX) ) printf( "-Infinity]\t" );
			else                                         printf( "%f]\t" , double(polys[i+1].start) );
			p=p+polys[i].p;
			p.printnl();
		}
	}
	printf("\n");
}
template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::ConstantFunction(const Real& radius)
{
	PPolynomial q;
	q.set(2);

	q.polys[0].start=-radius;
	q.polys[1].start= radius;

	q.polys[0].p.coefficients[0]= 1.0;
	q.polys[1].p.coefficients[0]=-1.0;
	return q;
}

template<> inline PPolynomial< 0 , float  > PPolynomial< 0 , float  >::GaussianApproximation( const float&  width ) { return ConstantFunction(width); }
template<> inline PPolynomial< 0 , double > PPolynomial< 0 , double >::GaussianApproximation( const double& width ) { return ConstantFunction(width); }

template< int Degree , class Real >
PPolynomial< Degree , Real > PPolynomial< Degree , Real >::GaussianApproximation(const Real& width)
{
	return PPolynomial< Degree-1 , Real >::GaussianApproximation().MovingAverage(width);
}
template< int Degree , class Real >
PPolynomial< Degree+1 , Real > PPolynomial< Degree , Real >::MovingAverage(const Real& radius)
{
	PPolynomial< Degree+1 , Real > A;
	Polynomial< Degree+1 , Real > p;
	StartingPolynomial< Degree+1 , Real >* sps;

	sps=(StartingPolynomial< Degree+1 , Real >*)malloc(sizeof(StartingPolynomial< Degree+1 , Real >)*polyCount*2);

	for(int i=0;i<(int)(polyCount);i++){
		sps[2*i  ].start=polys[i].start-radius;
		sps[2*i+1].start=polys[i].start+radius;
		p=polys[i].p.integral()-polys[i].p.integral()(polys[i].start);
		sps[2*i  ].p=p.shift(-radius);
		sps[2*i+1].p=p.shift( radius)*-1;
	}
	A.set(sps,(int)(polyCount*2));
	free(sps);
	return A*Real( 1.0 ) / ( Real( 2. ) * radius );
}

template< int Degree , class Real >
void PPolynomial< Degree , Real >::getSolutions(const Real& c,std::vector<Real>& roots,const Real& EPS,const Real& min,const Real& max) const{
	Polynomial< Degree , Real > p;
	std::vector<Real> tempRoots;

	p.setZero();
	for(size_t i=0;i<polyCount;i++){
		p+=polys[i].p;
		if(polys[i].start>max){break;}
		if(i<polyCount-1 && polys[i+1].start<min){continue;}
		p.getSolutions(c,tempRoots,EPS);
		for(size_t j=0;j<tempRoots.size();j++){
			if(tempRoots[j]>polys[i].start && (i+1==polyCount || tempRoots[j]<=polys[i+1].start)){
				if(tempRoots[j]>min && tempRoots[j]<max){roots.push_back(tempRoots[j]);}
			}
		}
	}
}

template< int Degree , class Real >
void PPolynomial< Degree , Real >::write(FILE* fp,const int& samples,const Real& min,const Real& max) const{
	fwrite(&samples,sizeof(int),1,fp);
	for(int i=0;i<samples;i++){
		Real x=min+i*(max-min)/(samples-1);
		float v=(*this)(x);
		fwrite(&v,sizeof(float),1,fp);
	}
}
