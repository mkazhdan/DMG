#include <winsock2.h>
#include <ws2tcpip.h>
#include "MultigridSolver.h"
#include "Util/MultiStreamIO.h"
#include "Util/Socket.h"
#include "Util/Half/half.h"

void SetPaddedSize( int width , int height , int minRes , int& paddedWidth , int& paddedHeight )
{
	int blockSize;
	long long _paddedWidth  = minRes;
	long long _paddedHeight = minRes;
	int domainW = FiniteElements1D<float,Type,Degree>::DomainSize( width );
	int domainH = FiniteElements1D<float,Type,Degree>::DomainSize( height );
	while( _paddedWidth <domainW )	_paddedWidth  *= 2;
	while( _paddedHeight<domainH )	_paddedHeight *= 2;

	blockSize = _paddedWidth/minRes;
	blockSize *= 2;
	while( _paddedWidth-blockSize>domainW ) _paddedWidth -= blockSize;

	blockSize = _paddedHeight/minRes;
	blockSize *= 2;
	while( _paddedHeight-blockSize>domainH ) _paddedHeight -= blockSize;

	_paddedWidth  = FiniteElements1D<float,Type,Degree>::Dimension(_paddedWidth);
	_paddedHeight = FiniteElements1D<float,Type,Degree>::Dimension(_paddedHeight);

	paddedWidth  = _paddedWidth;
	paddedHeight = _paddedHeight;
}
bool IsValidBandSize( int width , int height , int iters )
{
	if( height< 2*Degree ) return false;
	if( width < 16*(iters+2)*WordPerDegree ) return false;
	return true;
}
bool IsDownSamplable( int width , int height , int iters , int& newWidth , int& newHeight )
{
	if( (width&3) || (height&3) || (width<20) )	return false;
	newWidth  = width>>1;
	newHeight = height>>1;
	return true;
}
bool IsDownSamplable(int width,int height,int& newWidth,int& newHeight)
{
	if( (width&3) || (height&3) || (width<20) )	return false;
	newWidth  = width>>1;
	newHeight = height>>1;
	return true;
}
bool IsDownSamplable(int dCount,ClientSocket* cSockets,int sockCount)
{
	for(int i=0;i<sockCount;i++)
		if ( (cSockets[i].clientData.start>>dCount)&3 || ((cSockets[i].clientData.end-cSockets[i].clientData.start)>>dCount)<20)	return false;
	return true;
}
template<class Real>
void SetRectangularLaplacianMatrix( SparseMatrix< Real > & lap , int width , int height , double iWeight=0 )
{
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotStencil1,d2DotStencil1,dotStencil2,d2DotStencil2;
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,dotStencil1,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,d2DotStencil1,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,dotStencil2,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,d2DotStencil2,1,1,false);
	SetRectangularLaplacianMatrix( dotStencil1 , d2DotStencil1 , dotStencil2 , d2DotStencil2 , lap , width , height , iWeight );
}
template<class Real>
void SetRectangularLaplacianMatrix( DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
								    SparseMatrix< Real > & lap , int width , int height , double iWeight=0 )
{
	lap.Resize( width*height );

	int ii , jj;
	int xStart , xEnd , yStart , yEnd;
	for(int y=0;y<height;y++)
	{
		if		(y<Degree)			jj = y							, yStart = -y		, yEnd = Degree;
		else if	(y>height-1-Degree)	jj = 2*Degree+(y-(height-1))	, yStart = -Degree	, yEnd = height-1-y;
		else						jj = Degree						, yStart = -Degree	, yEnd = Degree;

		for(int x=0;x<width;x++)
		{
			if		(x<Degree)			ii = x						, xStart = -x		, xEnd = Degree;
			else if	(x>width-1-Degree)	ii = 2*Degree+(x-(width-1))	, xStart = -Degree	, xEnd = width-1-x;
			else						ii = Degree					, xStart = -Degree	, xEnd = Degree;

			int idx = x+y*width;
			lap.SetGroupSize( idx , (yEnd-yStart+1)*(xEnd-xStart+1) );
			int _i = 0;
			for(int yy=yStart;yy<=yEnd;yy++)
				for(int xx=xStart;xx<=xEnd;xx++)
				{
					lap.m_ppElements[idx][_i  ].N=(x+xx)+(y+yy)*width;
					lap.m_ppElements[idx][_i++].Value=
						dotMajor.caseTable[ii].values[xx+Degree]*d2DotMinor.caseTable[jj].values[yy+Degree]+
						d2DotMajor.caseTable[ii].values[xx+Degree]*dotMinor.caseTable[jj].values[yy+Degree]+
						(dotMajor.caseTable[ii].values[xx+Degree]*dotMinor.caseTable[jj].values[yy+Degree])*iWeight;
				}
		}
	}
}

template<class Real>
void SetSphericalLaplacianMatrix(SparseMatrix<Real>& lap,int width,int height)
{
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotStencil1,d2DotStencil1,dotStencil2,d2DotStencil2;
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,dotStencil1,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,d2DotStencil1,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,dotStencil2,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,d2DotStencil2,1,1,false);
	SetSphericalLaplacianMatrix( dotStencil1 , d2DotStencil1 , dotStencil2 , d2DotStencil2 , lap , width , height );
}
template<class Real>
void SetSphericalLaplacianMatrix(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
								 SparseMatrix<Real>& lap,int width,int height)
{
	lap.Resize(width*height);

	for(int x=0;x<width;x++)
		for(int y=0;y<height;y++)
		{
			int idx1 = x + y*width;
			lap.SetGroupSize( idx1 , (2*Degree+1)*(2*Degree+1) );
			for(int i=-Degree;i<=Degree;i++)
				for(int j=-Degree;j<=Degree;j++)
				{
					int idx2 = (i+Degree)+(j+Degree)*(2*Degree+1);
					int xx=(x+i+width)%width;
					int yy= y+j;
					if(yy<0)
					{
						yy = -yy-1;
						if(xx<width/2)	xx = width/2-1-xx;
						else			xx = width/2+width-1-xx;
					}
					else if(yy>=height)
					{
						yy = height-1-(yy-height);
						xx = width-1-xx;
					}
					lap.m_ppElements[idx1][idx2].N = xx+yy*width;
					lap.m_ppElements[idx1][idx2].Value = 
						dotMajor.caseTable[Degree].values[i+Degree]*d2DotMinor.caseTable[Degree].values[j+Degree]+
						d2DotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree];
				}
		}
}
template<class Real>
void SetSphericalLaplacianMatrix(SparseMatrix<Real>& lap,int width,int height,double iWeight)
{
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotStencil1,d2DotStencil1,dotStencil2,d2DotStencil2;
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,dotStencil1,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,d2DotStencil1,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,dotStencil2,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,d2DotStencil2,1,1,false);
	SetSphericalLaplacianMatrix( dotStencil1 , d2DotStencil1 , dotStencil2 , d2DotStencil2 , lap , width , height , iWeight );
}
template<class Real>
void SetSphericalLaplacianMatrix(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
								 SparseMatrix<Real>& lap,int width,int height,double iWeight)
{

	lap.Resize(width*height);

	for(int x=0;x<width;x++)
		for(int y=0;y<height;y++)
		{
			int idx1 = x + y*width;
			lap.SetGroupSize( idx1 , (2*Degree+1)*(2*Degree+1) );
			for(int i=-Degree;i<=Degree;i++)
				for(int j=-Degree;j<=Degree;j++)
				{
					int idx2 = (i+Degree)+(j+Degree)*(2*Degree+1);
					int xx=(x+i+width)%width;
					int yy= y+j;
					if(yy<0)
					{
						yy = -yy-1;
						if(xx<width/2)	xx = width/2-1-xx;
						else			xx = width/2+width-1-xx;
					}
					else if(yy>=height)
					{
						yy = height-1-(yy-height);
						xx = width-1-xx;
					}
					lap.m_ppElements[idx1][idx2].N = xx+yy*width;
					lap.m_ppElements[idx1][idx2].Value = 
						dotMajor.caseTable[Degree].values[i+Degree]*d2DotMinor.caseTable[Degree].values[j+Degree]+
						d2DotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree]+
						(dotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree])*iWeight;
				}
		}
}

DWORD WINAPI AcceptThread(LPVOID lpParam)
{
	AcceptThreadData* params = (AcceptThreadData*)lpParam;
	*params->connectSocket=accept(*params->listenSocket,NULL,NULL);
	if (*params->connectSocket == INVALID_SOCKET)
	{
		fprintf( stderr , "accept failed: %s\n" , LastSocketError() );
		params->success = false;
	}
	else
	{
		int val=1;
		setsockopt(*params->connectSocket,IPPROTO_TCP,TCP_NODELAY,(char*)&val,sizeof(val));
		closesocket(*params->listenSocket);
		*params->listenSocket = INVALID_SOCKET;
		params->success = true;
	}
	return 0;
}
///////////////////////////
// SphericalSynchronizer //
///////////////////////////
template<int Channels>
SphericalSynchronizer<Channels>::SphericalSynchronizer(void)
{
	_cCount = 0;
	_clientSockets = NULL;
}
template<int Channels>
void SphericalSynchronizer<Channels>::init(int width,const int* widths,int cCount,int rPasses,int pPasses,int vCycles)
{
	_width=width;
	_cCount=cCount;
	_rPasses=rPasses;
	_pPasses=pPasses;
	_vCycles=vCycles;

	if(_clientSockets)	delete[] _clientSockets;

	_clientSockets = new ClientSocket[cCount];
	if(!_clientSockets)
	{
		fprintf(stderr,"Failed to allocate SphericalSynchronizer::_clientSockets\n");
		exit(0);
	}
	for(int i=0;i<cCount;i++)
	{
		_clientSockets[i].client = GetListenSocket(_clientSockets[i].port);
		if (_clientSockets[i].client == INVALID_SOCKET)	fprintf(stderr,"Failed to set server in SphericalSynchronizer\n")	,	exit(0);
		if(i)	_clientSockets[i].start = _clientSockets[i-1].end;
		else	_clientSockets[i].start = 0;
		_clientSockets[i].end = _clientSockets[i].start+widths[i];
	}
}
template<int Channels>
SphericalSynchronizer<Channels>::~SphericalSynchronizer(void)
{
	delete[] _clientSockets;
	_clientSockets=NULL;
}
template<int Channels>	int SphericalSynchronizer<Channels>::clients(void)		const	{ return _cCount; }
template<int Channels>	int SphericalSynchronizer<Channels>::port(int cIndex)	const	{ return _clientSockets[cIndex].port; }
template<int Channels>
template<class DataType>
void SphericalSynchronizer<Channels>::Run(void)
{
	for(int i=0;i<_cCount;i++)
	{
		SOCKET temp = accept(_clientSockets[i].client,NULL,NULL);
		if (temp == INVALID_SOCKET)	fprintf( stderr , "accept failed in SphericalSynchronizer::Run(): %s\n" , LastSocketError() )	,	exit(0);
		int val=1;
		setsockopt(temp,IPPROTO_TCP,TCP_NODELAY,(char*)&val,sizeof(val));
		closesocket(_clientSockets[i].client);
		_clientSockets[i].client=temp;
	}

	DataType *buffer1,*buffer2;
	buffer1 = new DataType[_width*Channels*Degree];
	buffer2 = new DataType[_width*Channels*Degree];
	for ( int v = 0 ; v < _vCycles ; v++ )
	{
		//////////////////////////////////////////
		// First synchronize on the restriction //
		//////////////////////////////////////////

		// Synchronize the head of the stream
		for(int k=0;k<_rPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on head restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width/2;x++)
				{
					int xx = _width/2-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
				for(int x=_width/2;x<_width;x++)
				{
					int xx = _width/2+_width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on head restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}
		// Synchronize the tail of the stream
		for(int k=0;k<_rPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on tail restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width;x++)
				{
					int xx = _width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on tail restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}

		/////////////////////////////////////////
		// Now synchronize on the prolongation //
		/////////////////////////////////////////

		// Synchronize the head of the stream
		for(int k=0;k<_pPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on head prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width/2;x++)
				{
					int xx = _width/2-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
				for(int x=_width/2;x<_width;x++)
				{
					int xx = _width/2+_width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on head prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}
		// Synchronize the tail of the stream
		for(int k=0;k<_pPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on tail prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width;x++)
				{
					int xx = _width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on tail prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}
	}
	delete[] buffer1;
	delete[] buffer2;
}
template<int Channels>
template<class DataType>
DWORD WINAPI SphericalSynchronizer<Channels>::RunThread(LPVOID lpParam)
{
	SphericalSynchronizer<Channels>* synchronizer = (SphericalSynchronizer*)lpParam;
	synchronizer->Run<DataType>();
	return 0;
}

/////////////////////////////
// SocketedMultigridServer //
/////////////////////////////
template<int Channels>
SocketedMultigridServer<Channels>::SocketedMultigridServer(void)
{
	// Initialize Winsock
	MyWinSock::Load();
	clientSockets=NULL;
	_icSynchronizers = _oocSynchronizers = NULL;
	_synchronizerHandles = NULL;
	_syncSockets = NULL;
}
template<int Channels>
SocketedMultigridServer<Channels>::~SocketedMultigridServer(void)
{
	if(clientSockets)	delete[] clientSockets;
	clientSockets=NULL;
	MyWinSock::UnLoad();
	if(_icSynchronizers)	delete[] _icSynchronizers;
	if(_oocSynchronizers)	delete[] _oocSynchronizers;
	_icSynchronizers = _oocSynchronizers = NULL;
	if(_globalData.spherical)
	{
		WaitForMultipleObjects(_oocLevels+1+_icDCount,_synchronizerHandles,true,INFINITE);
		for(int i=0;i<_oocLevels+1+_icDCount;i++)	CloseHandle(_synchronizerHandles[i]);
		delete[] _synchronizerHandles;
		_synchronizerHandles = NULL;
		if(_syncSockets)
		{
			for(int i=0;i<_icDCount;i++)	if(_syncSockets[i] != INVALID_SOCKET)	closesocket(_syncSockets[i]);
			delete[] _syncSockets;
		}
		_syncSockets = NULL;
	}
}
template<int Channels>
bool SocketedMultigridServer<Channels>::SetUp( int port , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles,
											   int quality , int lanes , bool verbose , bool spherical , bool sharpen , double iWeight , double gWeight ,
											   bool showProgress , bool noCG )
{
	SOCKET listenSocket = INVALID_SOCKET;
	clientSockets = new ClientSocket[clientCount];

	_noCG		= noCG;
	_clientCount= clientCount;
	_inCoreRes	= inCoreRes;
	_minMGRes	= minMGRes;

	_globalData.spherical	= spherical;
	_globalData.verbose		= verbose;
	_globalData.iters		= iters;
	_globalData.vCycles		= vCycles;
	_globalData.sharpen		= sharpen;
	_globalData.iWeight		= iWeight;
	_globalData.gWeight		= gWeight;
	_globalData.quality		= quality;
	_globalData.lanes		= lanes;
	_globalData.progress	= showProgress;

	// Create a listening SOCKET for connecting to server
	listenSocket = GetListenSocket( port );
	if ( listenSocket == INVALID_SOCKET ) return false;

	char hostAddress[512];
	GetHostAddress( hostAddress );
	printf( "Server Address: %s:%d\n", hostAddress , port );

	// Establish a connection to the clients
	for( int i=0 ; i<clientCount ; i++ ) clientSockets[i].client  = AcceptSocket( listenSocket );

	// Get the client info
	int width = 0 , height;
	for( int i=0 ; i<clientCount ; i++ )
	{
		ReceiveOnSocket( clientSockets[i].client , &clientSockets[i].clientData , sizeof(ClientData) ) ;
		if(!i) height = clientSockets[i].clientData.height;
		else if( height!=clientSockets[i].clientData.height ) fprintf( stderr , "Band heights differ: %d != %d\n" , height , clientSockets[i].clientData.height ) , exit(0);
		width += clientSockets[i].clientData.width;
	}

	// Sort the clients by index order
	qsort( clientSockets , clientCount , sizeof(ClientSocket) , ClientSocket::Sort );

	// Store the clipping width and then computed the padded dimensions
	_globalData.cWidth  = width;
	_globalData.cHeight = height;

	if(!spherical)	_globalData.width = width+1, _globalData.height = height+1;
	else			_globalData.width = width  , _globalData.height = height;

	SetPaddedSize( _globalData.width , _globalData.height , _minMGRes, _globalData.width , _globalData.height );

	// Send the general solver information
	for( int i=0 ; i<clientCount ; i++ ) SendOnSocket( clientSockets[i].client , &_globalData , sizeof(_globalData) );

	// Set the local information for the bands
	int start=0;
	for( int i=0 ; i<clientCount ; i++ )
	{
		int w = clientSockets[i].clientData.width;
		clientSockets[i].clientData.start = start;
		clientSockets[i].clientData.end = clientSockets[i].clientData.cEnd = start+w;
		if( i==clientCount-1 ) clientSockets[i].clientData.end  = _globalData.width;
		start = clientSockets[i].clientData.end;
		SendOnSocket( clientSockets[i].client , &clientSockets[i].clientData , sizeof(ClientData) );
	}

	// Figure out how many different out-of-core levels are needed
	_oocLevels = 1;
	while( long long( width>>_oocLevels ) * long long( height>>_oocLevels ) >= _inCoreRes*_inCoreRes ) _oocLevels++;

	// WARNING!!! It's possible that _oocLevels = 0 if the image fits in-core
	std::vector< std::vector< std::pair<int,int> > > bandSizes;
	bandSizes.resize( _oocLevels );

	// Load in the initial band sizes and check that they are not too small.
	bandSizes[0].resize( _clientCount );
	for( int j=0 ; j<_clientCount ; j++ )
	{
		int width = clientSockets[j].clientData.end - clientSockets[j].clientData.start;
		if( !IsValidBandSize( width , _globalData.height , _globalData.iters ) )
		{
			fprintf( stderr , "Initial stream[%d] is too small: %d x %d\n" , j , width , _globalData.height );
			return false;
		}
		bandSizes[0][j].first = width;
		bandSizes[0][j].second = 0;
	}

	// Now check if we need to merge bands as we go to coarser resolutions
	for( int i=1 ; i<_oocLevels ; i++ )
	{
		bool downSample = false;
		for( int j=0 ; j<bandSizes[i-1].size() ; j++ ) if( !IsValidBandSize( bandSizes[i-1][j].first/2 , _globalData.height>>i , _globalData.iters ) ) downSample = true;
		if( downSample )
		{
			bandSizes[i].resize( bandSizes[i-1].size()>>1 );
			for( int j=0 ; j<bandSizes[i  ].size() ; j++ ) bandSizes[i][j].first = bandSizes[i][j].second = 0;
			for( int j=0 ; j<bandSizes[i-1].size() ; j++ ) bandSizes[i][j>>1].first += bandSizes[i-1][j].first , bandSizes[i][j>>1].second++;
			for( int j=0 ; j<bandSizes[i  ].size() ; j++ )
			{
				bandSizes[i][j].first >>= 1;
				if( !IsValidBandSize( bandSizes[i][j].first , _globalData.height>>i , _globalData.iters ) )
				{
					fprintf( stderr , "Stream[%d][%d] is too small: %d x %d\n" , i , j , bandSizes[i][j].first , _globalData.height>>i );
					return false;
				}
			}
		}
		else
		{
			bandSizes[i].resize( bandSizes[i-1].size() );
			for( int j=0 ; j<bandSizes[i-1].size() ; j++ ) bandSizes[i][j].first = bandSizes[i-1][j].first>>1 , bandSizes[i][j].second = 1;
		}
	}


	// Find the levels at which we change the number of threads
	std::vector< int > xLevels;
	xLevels.push_back( 0 );
	for( int i=1 ; i<_oocLevels ; i++ ) if( bandSizes[i].size() != bandSizes[i-1].size() ) xLevels.push_back( i );
for( int i=0 ; i<bandSizes.size() ; i++ )
{
	for( int j=0 ; j<bandSizes[i].size() ; j++ ) printf( "%d\t" , bandSizes[i][j].first);
	printf("\n");
}
for( int i=0 ; i<xLevels.size() ; i++ ) printf( "%d] [%d x %d] %d\n" , i , _globalData.width>>xLevels[i] , _globalData.height>>xLevels[i] , bandSizes[xLevels[i]].size() );

	// Construct the thread information
	processSockets.resize( xLevels.size() );
	for( int i=0 ; i<xLevels.size() ; i++ ) processSockets[i].resize( bandSizes[ xLevels[i] ].size() );

	for( int i=0 ; i<_clientCount ; i++ )
	{
		int idx = i;
		int tCount = 1;
		for( int j=1 ; j<xLevels.size() && !(idx&1) ; j++ ) tCount++ , idx>>=1;
		SendOnSocket( clientSockets[i].client , &tCount , sizeof(tCount) );

		ProcessData *pd = new ProcessData[tCount];
		pd[0].index = 0 ;
		pd[0].offset = i;
		pd[0].startDepth = 0;
		pd[0].start = 0;
		for( int k=0 ;  k<i ; k++ ) pd[0].start += bandSizes[0][k].first;
		pd[0].stop = pd[0].start + bandSizes[0][i].first;
		pd[0].height = _globalData.height;
		pd[0].children = bandSizes[0][i].second;

		idx = i;
		for( int j=1 ; j<tCount ; j++ )
		{
			idx >>= 1;
			pd[j].index  = j;
			pd[j].offset = idx;
			pd[j-1].endDepth = pd[j].startDepth = xLevels[j]-1;
			pd[j].start = 0;
			for( int k=0 ; k<idx ; k++ ) pd[j].start += bandSizes[xLevels[j]][k].first;
			pd[j].stop  = pd[j].start + bandSizes[xLevels[j]][idx].first;
			pd[j].height = _globalData.height >> xLevels[j];
			pd[j].children = bandSizes[xLevels[j]][idx].second;
		}
		if( tCount<xLevels.size() ) pd[tCount-1].endDepth = xLevels[tCount]-1;
		else						pd[tCount-1].endDepth = _oocLevels-1;
		SendOnSocket( clientSockets[i].client , pd , sizeof(ProcessData) * tCount );
		delete[] pd;
	}
	/*
	for( int i=0 ; i<processSockets.size() ; i++ )
		for( int j=0 ; j<processSockets[i].size() ; j++ )
		{
			// Establish a connection to the process
			int index , offset;
			SOCKET sock = AcceptSocket( listenSocket );
			ReceiveOnSocket( sock , &index , sizeof(index) );
			ReceiveOnSocket( sock , &offset , sizeof(offset) );
			processSockets[index][offset].sock = sock;
			ProcessData& pData = processSockets[index][offset].pData;
			if( index )
			{
				ReceiveOnSocket ( sock , &pData.childAddress , sizeof(&pData.childAddress) );
				ReceiveOnSocket ( sock , &pData.childPort , sizeof(&pData.childPort) );
			}
			if( spherical || offset )
			{
				ReceiveOnSocket ( sock , &pData.leftAddress , sizeof(&pData.leftAddress) );
				ReceiveOnSocket ( sock , &pData.leftPort , sizeof(&pData.leftPort) );
			}
		}
	for( int i=0 ; i<processSockets.size() ; i++ )
		for( int j=0 ; j<processSockets[i].size() ; j++ )
		{
			ProcessData& pData = processSockets[i][j].pData;
			{
				// ERROR!!! You have to special-case the lowest-resolution ooc processors
				SOCKET& sock = processSockets[i+1][j>>1].sock;
				SendOnSocket ( sock , &pData.childAddress , sizeof(&pData.childAddress) );
				SendOnSocket ( sock , &pData.childPort , sizeof(&pData.childPort) );
			}
			if( spherical || j<processSockets[i].size()-1 )
			{
				SOCKET& sock = processSockets[i][(j+1)%processSockets[i].size()].sock;
				SendOnSocket ( sock , &pData.leftAddress , sizeof(&pData.leftAddress) );
				SendOnSocket ( sock , &pData.leftPort , sizeof(&pData.leftPort) );
			}
		}
		*/
/*
	for( int i=0 ; i<processSockets.size() ; i++ )
		if(spherical)
		{
			// Get the addresses and ports of the clients and pass them on to the neighbors
			for( int j=0 ; j<processSockets[i].size() ; j++ )
			{
				ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
				ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
			}
			for(int i=0;i<clientCount;i++)
			{
				SendOnSocket ( clientSockets[(i-1+clientCount)%clientCount].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
				SendOnSocket ( clientSockets[(i-1+clientCount)%clientCount].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
			}
	}
	else
	{
		// Get the addresses and ports of the clients and pass them on to the neighbors
		for(int i=1;i<clientCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
		for(int i=1;i<clientCount;i++)
		{
			SendOnSocket ( clientSockets[i-1].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			SendOnSocket ( clientSockets[i-1].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
	}
*/
	////////////////
	// UP TO HERE //
	////////////////
#if 0
	// Send the image dimensions to the client
	for( int i=0 ; i<clientCount ; i++ )
	{
		ServerData sd;
		sd.width  = _width;
		sd.height = _height;
		sd.start  = clientSockets[i].start;
		sd.end    = clientSockets[i].end;
		if(_cWidth>sd.end)						sd.cEnd=sd.end;
		else if(_cWidth>clientSockets[i].start)	sd.cEnd=_cWidth;
		else									sd.cEnd=clientSockets[i].start;
		sd.iters=_iters;
		sd.vCycles=_vCycles;
		sd.quality=quality;
		sd.lanes=lanes;
		sd.verbose=verbose;
		sd.spherical=spherical;
		sd.progress=showProgress;
		sd.sharpen=sharpen;
		sd.iWeight=iWeight;
		sd.gWeight=gWeight;
		SendOnSocket ( clientSockets[i].client,&sd,sizeof(sd) );
	}
	if(spherical)
	{
		// Get the addresses and ports of the clients and pass them on to the neighbors
		for(int i=0;i<clientCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
		for(int i=0;i<clientCount;i++)
		{
			SendOnSocket ( clientSockets[(i-1+clientCount)%clientCount].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			SendOnSocket ( clientSockets[(i-1+clientCount)%clientCount].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
	}
	else
	{
		// Get the addresses and ports of the clients and pass them on to the neighbors
		for(int i=1;i<clientCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
		for(int i=1;i<clientCount;i++)
		{
			SendOnSocket ( clientSockets[i-1].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			SendOnSocket ( clientSockets[i-1].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
	}
	// Now everybody is talking to each other...

	// Figure out how many different out-of-core resolutions are needed
	_oocDCount=0;
	while(IsDownSamplable(_oocDCount,clientSockets,clientCount) && long long(width>>_oocDCount)*long long(height>>_oocDCount)>=inCoreRes*inCoreRes)	_oocDCount++;
	for(int i=0;i<clientCount;i++)	SendOnSocket ( clientSockets[i].client,&_oocDCount,sizeof(_oocDCount) );

	// Figure out how many different in-core resolutions are needed
	_icDCount = 1;
	int w = _width>>(_oocDCount-1) , h = _height>>(_oocDCount-1);
	int ww = w , hh = h , ww2 , hh2;
	while(IsDownSamplable(ww,hh,ww2,hh2) && ww2*hh2>=_minMGRes*_minMGRes && ww2>=16 && hh2>=16)
//	while(IsDownSamplable(ww,hh,ww2,hh2) && ww2*hh2>=_minMGRes*_minMGRes)
	{
		_icDCount++;
		ww=ww2;
		hh=hh2;
	}

	if(_spherical)
	{
		_oocSynchronizers = new SphericalSynchronizer<Channels>[2*_oocDCount+1];
		_icSynchronizers = new SphericalSynchronizer<Channels>[2*_icDCount];
		if(!_oocSynchronizers || !_icSynchronizers)
		{
			fprintf(stderr,"Failed to allocate synchronizers\n");
			return false;
		}

		int *oocWidths,icWidths;
		oocWidths = new int[_cCount];

		for(int i=0;i<_cCount;i++)	oocWidths[i] = (clientSockets[i].end-clientSockets[i].start);
		_oocSynchronizers[0].init(_width,oocWidths,_cCount,1,0,1);
		for(int d=0;d<_oocDCount;d++)
		{
			for(int i=0;i<_cCount;i++)	oocWidths[i] = (clientSockets[i].end-clientSockets[i].start)>>d;
			_oocSynchronizers[2*d+1].init(_width>>d,oocWidths,_cCount,_iters*Degree+1,_iters*Degree+1,_vCycles);
			if(_verbose)	_oocSynchronizers[2*d+2].init(_width>>d,oocWidths,_cCount,Degree,Degree,_vCycles);
			else			_oocSynchronizers[2*d+2].init(_width>>d,oocWidths,_cCount,Degree,     0,_vCycles);
		}
		delete[] oocWidths;

		for(int d=0;d<_icDCount;d++)
		{
			icWidths = _width>>(d+_oocDCount-1);
			_icSynchronizers[2*d].init(icWidths,&icWidths,1,_iters*Degree+1,_iters*Degree+1,_vCycles);
			if(_verbose)	_icSynchronizers[2*d+1].init(icWidths,&icWidths,1,Degree,Degree,_vCycles);
			else			_icSynchronizers[2*d+1].init(icWidths,&icWidths,1,Degree,     0,_vCycles);
		}
		_synchronizerHandles = new HANDLE[2*(_oocDCount+_icDCount)+1];
		{
			DWORD synchThreadID;
			if( sharpen )
			{
				_synchronizerHandles[0]=CreateThread(
					NULL,												// default security attributes
					0,													// use default stack size  
					SphericalSynchronizer<Channels>::RunThread<float>,	// thread function 
					&_oocSynchronizers[0],								// argument to thread function 
					0,													// use default creation flags 
					&synchThreadID);									// returns the thread identifier
			}
			else
			{
				_synchronizerHandles[0]=CreateThread(
					NULL,																			// default security attributes
					0,																				// use default stack size  
					SphericalSynchronizer<Channels>::RunThread<ImageData<float,unsigned __int16> >,	// thread function 
					&_oocSynchronizers[0],															// argument to thread function 
					0,																				// use default creation flags 
					&synchThreadID);																// returns the thread identifier
			}
			if(!_synchronizerHandles[0])	fprintf(stderr,"Failed to create sync thread\n")	,	exit(0);
		}
		for(int i=0;i<2*_oocDCount;i++)
		{
			DWORD synchThreadID;
			_synchronizerHandles[i+1]=CreateThread( 
				NULL,												// default security attributes
				0,													// use default stack size  
				SphericalSynchronizer<Channels>::RunThread<float>,	// thread function 
				&_oocSynchronizers[i+1],							// argument to thread function 
				0,													// use default creation flags 
				&synchThreadID);									// returns the thread identifier
			if(!_synchronizerHandles[i+1])	fprintf(stderr,"Failed to create sync thread\n")	,	exit(0);
		}
		for(int i=0;i<2*_icDCount;i++)
		{
			DWORD synchThreadID;
			_synchronizerHandles[i+2*_oocDCount+1]=CreateThread( 
				NULL,												// default security attributes
				0,													// use default stack size  
				SphericalSynchronizer<Channels>::RunThread<float>,	// thread function 
				&_icSynchronizers[i],								// argument to thread function 
				0,													// use default creation flags 
				&synchThreadID);									// returns the thread identifier
			if(!_synchronizerHandles[i+2*_oocDCount+1])	fprintf(stderr,"Failed to create sync thread\n")	,	exit(0);
		}
		int* oocPorts = new int[2*_oocDCount+1];
		int* icPorts = new int[2*_icDCount];
		for(int d=0;d<2*_icDCount;d++)	icPorts[d] = _icSynchronizers[d].port(0);
		for(int i=0;i<_cCount;i++)
		{
			for(int d=0;d<2*_oocDCount+1;d++)	oocPorts[d] = _oocSynchronizers[d].port(i);
			SendOnSocket( clientSockets[i].client , oocPorts , sizeof(int)*2*(_oocDCount+1) );
		}
		delete[] oocPorts;

		_syncSockets = new SOCKET[2*_icDCount];
		for(int i=0;i<2*_icDCount;i++)
		{
			_syncSockets[i] = GetConnectSocket(hostAddress,icPorts[i]);
			if(_syncSockets[i] == INVALID_SOCKET)	return false;
		}
		delete[] icPorts;
	}

	closesocket(listenSocket);
	listenSocket = INVALID_SOCKET;
#endif
	return true;
}
template<int Channels>
void SocketedMultigridServer<Channels>::Run(void)
{
#if 0
	printf("Running server\n");
	DotProductStencil stencils[4];
	SolverInfo<Channels>*  solverInfo = new SolverInfo<Channels>[_oocDCount];
	SolverInfo<Channels>* tSolverInfo = new SolverInfo<Channels>[_oocDCount];
	double zeroAverage[Channels];
	for(int c=0;c<Channels;c++)	zeroAverage[c]=0;
	Vector<float> globalB,globalX;
	int gWidth  = _width  >> (_oocDCount-1);
	int gHeight = _height >> (_oocDCount-1);

	globalB.Resize( gWidth * gHeight * Channels );
	globalX.Resize( gWidth * gHeight * Channels );

	if(_verbose)	MultigridSolver<float,ZERO_DERIVATIVE,2,Channels>::verbose = MultigridSolver<float,ZERO_DERIVATIVE,2,Channels>::FULL_VERBOSE;

	memset(_average, 0 , sizeof(_average) );

	for(int ii=0;ii<_vCycles;ii++)
	{
		double t;

		t=Time();
		memset(solverInfo,0,sizeof(SolverInfo<Channels>)*_oocDCount);

#if SOCKET_BACKED_GRID
		for( int j=0 ; j<gHeight ; j++ )
			for(int i=0;i<_cCount;i++)
			{
				Vector<float> localB,localX;
				int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
				int lOffset = clientSockets[i].start>>(_oocDCount-1);
				localX.Resize( lWidth * Channels );
				localB.Resize( lWidth * Channels );
				ReceiveOnSocket( clientSockets[i].clientB , &localB[0] , sizeof(float) * lWidth * Channels , true );
				ReceiveOnSocket( clientSockets[i].clientX , &localX[0] , sizeof(float) * lWidth * Channels , true );
				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*j*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*j*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*c];
					float *pLocalX  = &localX [lWidth*c];
					memcpy(pGlobalB,pLocalB,lWidth*sizeof(float));
					memcpy(pGlobalX,pLocalX,lWidth*sizeof(float));
				}
		}
#endif // SOCKET_BACKED_GRID
		for(int i=0;i<_cCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client , tSolverInfo , sizeof(SolverInfo<Channels>) * _oocDCount );
			for(int j=0;j<_oocDCount;j++)
			{
				solverInfo[j].bSquareNorm += tSolverInfo[j].bSquareNorm;
				solverInfo[j].rSquareNorm += tSolverInfo[j].rSquareNorm;
				solverInfo[j].xSquareNorm += tSolverInfo[j].xSquareNorm;
				for(int c=0;c<Channels;c++)	solverInfo[j].solutionSum[c] += tSolverInfo[j].solutionSum[c];
			}
			if(ii==0)
			{
				double avg[Channels];
				ReceiveOnSocket ( clientSockets[i].client , avg ,sizeof(avg) );
				for(int c=0 ; c<Channels ; c++)	_average[c]+=avg[c];
			}
		}
		if(ii==0)	for(int c=0 ; c<Channels ; c++)	_average[c] /= _cWidth	,	_average[c] /= _cHeight;

		// Get the low-res data off of the sockets
#if !SOCKET_BACKED_GRID
		for(int i=0;i<_cCount;i++)
		{
			Vector<float> localB,localX;
			int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
			int lHeight = _height>>(_oocDCount-1);
			int lOffset = clientSockets[i].start>>(_oocDCount-1);
			localX.Resize( lWidth*lHeight*Channels );
			localB.Resize( lWidth*lHeight*Channels );
			ReceiveOnSocket (clientSockets[i].client,&localB[0],sizeof(float)*localB.Dimensions() ) ;
			ReceiveOnSocket (clientSockets[i].client,&localX[0],sizeof(float)*localX.Dimensions() ) ;
			for(int y=0;y<gHeight;y++)
				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*y*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*y*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*y*Channels+lWidth*c];
					float *pLocalX  = &localX [lWidth*y*Channels+lWidth*c];
					memcpy(pGlobalB,pLocalB,lWidth*sizeof(float));
					memcpy(pGlobalX,pLocalX,lWidth*sizeof(float));
				}
		}
#endif // !SOCKET_BACKED_GRID
		for(int i=0;i<_cCount;i++) ReceiveOnSocket ( clientSockets[i].client , stencils , sizeof( DotProductStencil ) * 4 );

		printf("Out-Of-Core Restriction:    %f\n",Time()-t);
		if(_verbose)
			for(int i=_oocDCount-1;i>=0;i--)
			{
				printf("\tError[%d x %d] %g -> %g\t",_width>>(_oocDCount-1-i),_height>>(_oocDCount-1-i),sqrt(solverInfo[i].bSquareNorm),sqrt(solverInfo[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}


		// Get the base solution
		if(ii==_vCycles-1)	SolveInCore(stencils[0],stencils[1],stencils[2],stencils[3],globalB,globalX,_average);
		else				SolveInCore(stencils[0],stencils[1],stencils[2],stencils[3],globalB,globalX,zeroAverage);
		// Send the low-res data to the sockets
#if SOCKET_BACKED_GRID
		for( int j=0 ; j<gHeight ; j++ )
			for(int i=0;i<_cCount;i++)
			{
				Vector<float> localB,localX;
				int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
				int lOffset = clientSockets[i].start>>(_oocDCount-1);
				localX.Resize( lWidth * Channels );
				localB.Resize( lWidth * Channels );

				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*j*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*j*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*c];
					float *pLocalX  = &localX [lWidth*c];
					memcpy(pLocalB,pGlobalB,lWidth*sizeof(float));
					memcpy(pLocalX,pGlobalX,lWidth*sizeof(float));
				}
				SendOnSocket( clientSockets[i].clientX , &localX[0] , sizeof(float) * lWidth * Channels , true );
			}
#else // !SOCKET_BACKED_GRID
		for(int i=0;i<_cCount;i++)
		{
			Vector<float> localB,localX;
			int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
			int lHeight = _height>>(_oocDCount-1);
			int lOffset = clientSockets[i].start>>(_oocDCount-1);
			localX.Resize( lWidth*lHeight*Channels );
			localB.Resize( lWidth*lHeight*Channels );

			for(int y=0;y<gHeight;y++)
				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*y*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*y*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*y*Channels+lWidth*c];
					float *pLocalX  = &localX [lWidth*y*Channels+lWidth*c];
					memcpy(pLocalB,pGlobalB,lWidth*sizeof(float));
					memcpy(pLocalX,pGlobalX,lWidth*sizeof(float));
				}
			SendOnSocket (clientSockets[i].client,&localB[0],sizeof(float)*localB.Dimensions() ) ;
			SendOnSocket (clientSockets[i].client,&localX[0],sizeof(float)*localX.Dimensions() ) ;
		}
#endif // SOCKET_BACKED_GRID
		t=Time();
		memset(solverInfo,0,sizeof(SolverInfo<Channels>)*_oocDCount);
		for(int i=0;i<_cCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client , tSolverInfo , sizeof(SolverInfo<Channels>)*_oocDCount);
			for(int j=0;j<_oocDCount;j++)
			{
				solverInfo[j].bSquareNorm += tSolverInfo[j].bSquareNorm;
				solverInfo[j].rSquareNorm += tSolverInfo[j].rSquareNorm;
				solverInfo[j].xSquareNorm += tSolverInfo[j].xSquareNorm;
				for(int c=0;c<Channels;c++)	solverInfo[j].solutionSum[c] += tSolverInfo[j].solutionSum[c];
			}
		}
		printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
		if(_verbose)
			for(int i=0;i<_oocDCount;i++)
			{
				printf("\tError[%d x %d] %g -> %g\t",_width>>(_oocDCount-1-i),_height>>(_oocDCount-1-i),sqrt(solverInfo[i].bSquareNorm),sqrt(solverInfo[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
	}
	delete[] solverInfo;
	delete[] tSolverInfo;
#endif
}
template<int Channels>
void SocketedMultigridServer<Channels>::SolveInCore(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,
													DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
													Vector<float>& in,Vector<float>& out,double average[Channels])
{
#if 0
	double t;
	int w = _width>>(_oocDCount-1) , h = _height>>(_oocDCount-1);
	SocketedMultiGridStreamingSolver<Channels>* solvers;

	Vector<float> lowB,lowX;
	StreamingGrid *B,*X;

	solvers=new SocketedMultiGridStreamingSolver<Channels>[_icDCount];

	for(int i=1;i<_icDCount;i++)	solvers[i].parent=&solvers[i-1];
	for(int i=0;i<_icDCount-1;i++)	solvers[i].rChild =solvers[i].pChild =&solvers[i+1];
	for(int i=0;i<_icDCount;i++)
	{
		solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
		for(int c=0;c<Channels;c++)	solvers[i].solutionSum[c]=0;
		solvers[i].setResidual=true;
	}
	SOCKET leftSocket = INVALID_SOCKET , rightSocket = INVALID_SOCKET;
	if(_spherical)
	{
		int port=0;
	    struct sockaddr_in right;
		SOCKET listenSocket=GetListenSocket(port);
		if (listenSocket == INVALID_SOCKET)	fprintf(stderr,"Failed to generate listen socket\n")	,	exit(0);

		char hostName[512];
		gethostname(hostName,512);
		hostent* host=gethostbyname(hostName);

		rightSocket = socket(AF_INET, SOCK_STREAM, 0);
		if ( rightSocket == INVALID_SOCKET )	fprintf(stderr,"Error at socket(): %s\n", LastSocketError() )	,	exit(0);

		memset(&right, 0, sizeof(right));
		right.sin_family = AF_INET;
		memcpy(&right.sin_addr,host->h_addr,sizeof(struct in_addr));
		right.sin_port= htons ( port );

		// Accept on the left socket
		AcceptThreadData params;
		HANDLE acceptThread = NULL;
		DWORD accepthThreadID;
		params.listenSocket = &listenSocket;
		params.connectSocket = &leftSocket;
		acceptThread=CreateThread( 
			NULL,				// default security attributes
			0,					// use default stack size  
			AcceptThread,		// thread function 
			&params,			// argument to thread function 
			0,					// use default creation flags 
			&accepthThreadID);	// returns the thread identifier
		if(!acceptThread)	fprintf(stderr,"Failed to create accept thread\n")	,	exit(0);

		while (connect( rightSocket, (const sockaddr*)&right, sizeof(right) ) == SOCKET_ERROR)
		{
			fprintf(stderr,"Connecting...\n");
			Sleep(1);
		}
		WaitForSingleObject(acceptThread,INFINITE);
		CloseHandle(acceptThread);
		if(!params.success)	fprintf(stderr,"Failed to accept self\n")	,	exit(0);

	}
	// BADNESS!!! For spherical domains, the server socket has not been set here.
	if( _sharpen )	solvers[_icDCount-1].Initialize(dotMajor,d2DotMajor,dotMinor,d2DotMinor,_iWeight,0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);
	else			solvers[_icDCount-1].Initialize(dotMajor,d2DotMajor,dotMinor,d2DotMinor,         0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);

	lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].major*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].major*Channels*sizeof(float),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].major*Channels*sizeof(float),solvers[0].minor);
	solvers[_icDCount-1].inB=new MemoryBackedGrid(&in[0],w*Channels*sizeof(float),h);
	solvers[_icDCount-1].inX=new MemoryBackedGrid(&out[0],w*Channels*sizeof(float),h);
	solvers[0].outB=B;

	// Data for the interleaved streaming multigrid
	t=Time();
	// Initialize
	solvers[_icDCount-1].InitRestriction();
	solvers[_icDCount-1].SetRestriction();
	solvers[_icDCount-1].SolveRestriction();

	// Clean up
	delete solvers[_icDCount-1].inB;
	solvers[_icDCount-1].inB=NULL;
	delete solvers[_icDCount-1].inX;
	solvers[_icDCount-1].inX=NULL;

	if(_verbose)
	{
		printf("In-Core Restriction: %f\n",Time()-t);
		for(int i=_icDCount-1;i>=0;i--)
		{
			printf("\tError[%d x %d] %g -> %g\t",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
			for(int c=0;c<Channels;c++)
			{
				if(c==0)	printf("(");
				else		printf(" ");
				printf("%g",((solvers[i].solutionSum[c])/solvers[i].major)/solvers[i].minor);
				if(c==Channels-1)	printf(")");
				else				printf("");
			}
			printf("\n");
		}
	}
	solvers[_icDCount-1].UnSetRestriction();
	// Get the base solution
	SparseMatrix<double> lap;
	if(_spherical)
		if(!_sharpen || _iWeight==0)	SetSphericalLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor);
		else							SetSphericalLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor,_iWeight);
	else
		if(!_sharpen || _iWeight==0)	SetRectangularLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor);
		else							SetRectangularLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor,_iWeight);
	{
		Vector<double> myLowX,myLowB;
		myLowB.Resize(solvers[0].major*solvers[0].minor);
		myLowX.Resize(solvers[0].major*solvers[0].minor);
		lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
		for(int c=0;c<Channels;c++)
		{
			for(int i=0;i<solvers[0].major;i++)
				for(int j=0;j<solvers[0].minor;j++)
					myLowB[i+j*solvers[0].major]=lowB[i+j*solvers[0].major*Channels+c*solvers[0].major];
		if(!_noCG)
			if(!_sharpen || _iWeight==0) MySolveConjugateGradient<Type>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
//			else						 MySolveConjugateGradient<ZERO_VALUE>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
			else						 MySolveConjugateGradient<Type>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);

			for(int i=0;i<solvers[0].major;i++)
				for(int j=0;j<solvers[0].minor;j++)
				{
					lowX[i+j*solvers[0].major*Channels+c*solvers[0].major]=myLowX[i+j*solvers[0].major];
					myLowX[i+j*solvers[0].major]=0;
				}
		}
	}

	for(int i=0;i<_icDCount;i++)
	{
		solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
		for(int c=0;c<Channels;c++)	solvers[i].solutionSum[c]=0;
		solvers[i].setResidual=_verbose;
	}

	// Solve the prolongation
	t=Time();
	solvers[0].inX=X;
	// Set the child dependencies
	solvers[_icDCount-1].outX=new MemoryBackedGrid(&out[0],w*Channels*sizeof(float),h);
	solvers[_icDCount-1].InitProlongation();
	solvers[0].SetProlongation();
	// Solve
	solvers[0].SolveProlongation();
	if(_verbose)
	{
		printf("In-Core Prolongation:    %f\n",Time()-t);
		for(int i=0;i<_icDCount;i++)
		{
			printf("\tError[%d x %d] %g -> %g\t",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
			for(int c=0;c<Channels;c++)
			{
				if(c==0)	printf("(");
				else		printf(" ");
				printf("%g",((solvers[i].solutionSum[c])/solvers[i].major)/solvers[i].minor);
				if(c==Channels-1)	printf(")");
				else				printf("");
			}
			printf("\n");
		}
	}
	solvers[0].UnSetProlongation();
	delete solvers[_icDCount-1].outX;
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;

	int ww=(long long(_cWidth) *w+_width -1)/_width;
	int hh=(long long(_cHeight)*h+_height-1)/_height;
	for(int c=0;c<Channels;c++)
	{
		double newAverage=0;
		int pCount=0;

		for(int j=0;j<hh;j++)
		{
			float* o=&out[j*w*Channels+c*w];
			for(int i=0;i<ww;i++)
			{
				newAverage+=o[i];
				pCount++;
			}
		}

		newAverage/=pCount;
		newAverage=average[c]-newAverage;

		if(!_sharpen || _iWeight==0)
			for(int j=0;j<h;j++)
			{
				float* o=&out[j*w*Channels+c*w];
				for(int i=0;i<w;i++)	o[i]+=newAverage;
			}
	}

	if(leftSocket!=INVALID_SOCKET)	closesocket(leftSocket);
	if(rightSocket!=INVALID_SOCKET)	closesocket(rightSocket);
#endif
}

/////////////////////////////
// SocketedMultigridClient //
/////////////////////////////
template< int Channels, class StorageType , class PixelType , class LabelType >
SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::SocketedMultigridClient(void)
{
	MyWinSock::Load();

	_serverSocket = INVALID_SOCKET;
	_processData = NULL;
	_processes = NULL;
	_processCount = 0;
}
template< int Channels, class StorageType , class PixelType , class LabelType >
SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::~SocketedMultigridClient(void)
{
	if( _serverSocket != INVALID_SOCKET )	closesocket( _serverSocket ) , _serverSocket = INVALID_SOCKET;
	if( _processData )	delete[] _processData	, _processData = NULL;
	if( _processes )	delete[] _processes		, _processes = NULL;
	_processCount = 0;

	MyWinSock::UnLoad();
}
template< int Channels, class StorageType , class PixelType , class LabelType >
bool SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::sharpen( void ) const { return _globalData.sharpen; }
template< int Channels, class StorageType , class PixelType , class LabelType >
int SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::quality( void ) const { return _globalData.quality; }
template< int Channels, class StorageType , class PixelType , class LabelType >
int SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::cSize( void ) const { return _clientData.cEnd - _clientData.start; }
template< int Channels, class StorageType , class PixelType , class LabelType >
int SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::size( void ) const { return _clientData.end - _clientData.start; }
template< int Channels, class StorageType , class PixelType , class LabelType >
int SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::cHeight( void ) const { return _globalData.cHeight; }
template< int Channels, class StorageType , class PixelType , class LabelType >
int SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::height( void ) const { return _globalData.height; }

template< int Channels, class StorageType , class PixelType , class LabelType >
bool SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::SetUp( char* address , int port , int index , int w , int h )
{
	SOCKET listenSocket = INVALID_SOCKET;
    struct sockaddr_in local;

	printf( "Connecting to server" );
	_serverSocket  = GetConnectSocket( address , port , 500 , true );

	// Send the band information
	_clientData.index	= index;
	_clientData.width	= w;
	_clientData.height	= h;
	SendOnSocket( _serverSocket , &_clientData , sizeof(_clientData) );

	// Get the global information
	ReceiveOnSocket( _serverSocket , &_globalData , sizeof(_globalData) );

	// Get the local information
	ReceiveOnSocket( _serverSocket , &_clientData , sizeof(_clientData) );
printf("(%d , %d) x (%d , %d)\n",_globalData.cWidth,_globalData.width,_globalData.cHeight,_globalData.height);
printf("[ %d , %d ] x [ %d , %d ]\n" , _clientData.start,_clientData.cEnd,0,_globalData.cHeight);
printf("[ %d , %d ] x [ %d , %d ]\n" , _clientData.start,_clientData.end,0,_globalData.height);


	// Get back the total number of threads this client is responsible for and the information about these threads
	ReceiveOnSocket( _serverSocket , &_processCount , sizeof(_processCount) );
	_processData = new ProcessData[ _processCount ];
	ReceiveOnSocket( _serverSocket , _processData , sizeof(ProcessData) * _processCount );


for( int i=0 ; i<_processCount ; i++ )
{
	printf( "%d] [ (%d--%d) x %d] [%d -- %d] (%d,%d)\t%d\n" , i , _processData[i].start , _processData[i].stop , _processData[i].height , _processData[i].startDepth , _processData[i].endDepth , _processData[i].index , _processData[i].offset , _processData[i].children );
}
	////////////////
	// UP TO HERE //
	////////////////

#if 0
	// Get back the other dimension data
	{
		ServerData sd;
		ReceiveOnSocket ( _serverSocket,&sd,sizeof(sd) );
		_width		= sd.width;
		_height		= sd.height;
		_end		= sd.end;
		_cEnd		= sd.cEnd;
		_start		= sd.start;
		_cHeight	= h;
s	}
	// Send information about the left socket
	if(_start!=0 || _spherical)
	{
		int port=0;
		listenSocket=GetListenSocket(port);
		if (listenSocket == INVALID_SOCKET)	return false;

		char hostName[512];
		gethostname(hostName,512);
		hostent* host=gethostbyname(hostName);

		SendOnSocket ( _serverSocket,host->h_addr,sizeof(struct in_addr) );
		SendOnSocket ( _serverSocket,&port,sizeof(port) );
//		printf("Sending left: %s : %d\n",inet_ntoa(*(struct in_addr*)host->h_addr),port);
	}
	// Get the information about the right socket
	if(_end!=_width || _spherical)
	{
		int port=0;
		_rightSocket = socket(AF_INET, SOCK_STREAM, 0);
		if ( _rightSocket == INVALID_SOCKET )
		{
			fprintf(stderr,"Error at socket(): %s\n", LastSocketError() );
			return false;
		}
		memset(&right, 0, sizeof(right));
		right.sin_family = AF_INET;
		ReceiveOnSocket ( _serverSocket,&right.sin_addr,sizeof(right.sin_addr) );
		ReceiveOnSocket ( _serverSocket,&port,sizeof(port) );
		right.sin_port= htons ( port );
//		printf("Reading Right: %s : %d\n",inet_ntoa(right.sin_addr),port);
	}

	// Accept on the left socket
	AcceptThreadData params;
	HANDLE acceptThread = NULL;
	DWORD accepthThreadID;
	params.listenSocket = &listenSocket;
	params.connectSocket = &_leftSocket;
	if(_start!=0 || _spherical)
	{
		acceptThread=CreateThread( 
			NULL,				// default security attributes
			0,					// use default stack size  
			AcceptThread,		// thread function 
			&params,			// argument to thread function 
			0,					// use default creation flags 
			&accepthThreadID);	// returns the thread identifier
		if(!acceptThread)
		{
			fprintf(stderr,"Failed to create accept thread\n");
			exit(0);
		}
	}

	// Connect on the right socket
	if(_end!=_width || _spherical)
	{
		while (connect( _rightSocket, (const sockaddr*)&right, sizeof(right) ) == SOCKET_ERROR)
		{
			fprintf(stderr,"Error at connect(): %s\n", LastSocketError() );
			fprintf(stderr,"retrying...\n");
			Sleep(1000);
		}
		int val=1;
		setsockopt(_rightSocket,IPPROTO_TCP,TCP_NODELAY,(char*)&val,sizeof(val));
	}

	if(acceptThread)
	{
		WaitForSingleObject(acceptThread,INFINITE);
		CloseHandle(acceptThread);
		if(!params.success)	return false;
	}

	// Get the depth of the multigrid solver
	ReceiveOnSocket ( _serverSocket,&_dCount,sizeof(_dCount) );

	if(_spherical)
	{
		int* ports = new int[2*_dCount+1];
		_syncSockets = new SOCKET[2*_dCount+1];
		ReceiveOnSocket ( _serverSocket , ports , sizeof(int)*2*(_dCount+1) );
		for(int i=0;i<2*_dCount+1;i++)
		{
			_syncSockets[i] = GetConnectSocket(address,ports[i]);
			if(_syncSockets[i] == INVALID_SOCKET)	return false;
		}

		delete[] ports;
	}
#endif
	return true;
}
template< int Channels, class StorageType , class PixelType , class LabelType >
void SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::Run(StreamingGrid* pixels,StreamingGrid* labels,StreamingGrid* out)
{
#if 0
	double t;
	DotProductStencil stencils[4];
	SolverInfo<Channels>* solverInfo = new SolverInfo<Channels>[_dCount];
	SocketedStreamingDivergence<Channels,PixelType,LabelType,StorageType>* sLaplacian = new SocketedStreamingDivergence<Channels,PixelType,LabelType,StorageType>();
	SocketedMultiGridStreamingSolver<Channels,StorageType>* solvers;
	SocketedStreamingSolver<Channels>::server.Reset();

#if !SOCKET_BACKED_GRID
	Vector<float> lowB,lowX;
#endif // !SOCKET_BACKED_GRID
	StreamingGrid *B,*X;
	// Also, the template should be StorageType, not Real
	MultiStreamIOClient* in;

	if(_vCycles>1)	in=new MultiStreamIOClient( size()*sizeof(float)*Channels , _height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true);
	else			in=NULL;
 
	solvers = new SocketedMultiGridStreamingSolver<Channels,StorageType>[_dCount];
	solvers[_dCount-1].showProgress=_progress;
	for(int i=1;i<_dCount;i++)		solvers[i].parent = &solvers[i-1];
	for(int i=0;i<_dCount-1;i++)	solvers[i].rChild =  solvers[i].pChild = &solvers[i+1];
	for(int i=0;i<_dCount;i++)		solvers[i].laneNum = _lanes;


	sLaplacian->parent = & solvers[_dCount-1];
	solvers[_dCount-1].rChild = sLaplacian;
	sLaplacian->Initialize(pixels,labels,_start,_end,_width,_height,_iters,_leftSocket,_syncSockets,_rightSocket,true,_spherical);

#if SOCKET_BACKED_GRID
	B = new SocketBackedGrid( _serverBSocket , solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor , true );
	X = new SocketBackedGrid( _serverXSocket , solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor , true );
#else // !SOCKET_BACKED_GRID
	lowX.Resize(solvers[0].size()*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].size()*solvers[0].minor*Channels);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
#endif // SOCKET_BACKED_GRID
	for(int ii=0;ii<_vCycles;ii++)
	{
		if(!ii)
		{
			pixels->SetServer(&SocketedStreamingSolver<Channels>::server);
			labels->SetServer(&SocketedStreamingSolver<Channels>::server);
		}
		/////////////////
		// RESTRICTION //
		/////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=true;
		}
		if(ii)	solvers[_dCount-1].inX=in;
		else	solvers[_dCount-1].inX=NULL;
		solvers[_dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		solvers[0].outX=X;
		solvers[_dCount-1].rChild=sLaplacian;
		// Data for the interleaved streaming multigrid
		t=Time();
		// Initialize
		if(sLaplacian)	sLaplacian->InitRestriction()	,	sLaplacian->SetRestriction();
		else	solvers[_dCount-1].InitRestriction()	,	solvers[_dCount-1].SetRestriction();
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii)	in->SetServer(&SocketedStreamingSolver<Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
		if(sLaplacian)	sLaplacian->SolveRestriction();
		else			solvers[_dCount-1].SolveRestriction();
		SocketedStreamingSolver<Channels>::server.Reset();

		printf("Out-Of-Core Restriction:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}
		if(_verbose)
			for(int i=_dCount-1;i>=0;i--)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		if(ii==_vCycles-1 && in)	delete in,	in=NULL;
		if(sLaplacian)
		{
			SendOnSocket ( _serverSocket , sLaplacian->average , sizeof(sLaplacian->average) );
			sLaplacian->UnSetRestriction();
			delete sLaplacian;
			sLaplacian=NULL;
		}
		else	solvers[_dCount-1].UnSetRestriction();

		t=Time();
#if SOCKET_BACKED_GRID
#else // !SOCKET_BACKED_GRID
		// Send the low res stuff over to the server so that it can perform the in-core solve.
		SendOnSocket (_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		SendOnSocket (_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;
#endif // !SOCKET_BACKED_GRID

		stencils[0] = solvers[0].dotMajor;
		stencils[1] = solvers[0].d2DotMajor;
		stencils[2] = solvers[0].dotMinor;
		stencils[3] = solvers[0].d2DotMinor;
		SendOnSocket ( _serverSocket , stencils , sizeof( DotProductStencil ) * 4 );

		// Get back the in-core solution.
#if SOCKET_BACKED_GRID
#else // !SOCKET_BACKED_GRID
		ReceiveOnSocket	(_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		ReceiveOnSocket	(_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;
#endif // !SOCKET_BACKED_GRID
		printf("In-Core:    %f\n",Time()-t);

		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=_verbose;
		}
		solvers[_dCount-1].inX=NULL;
		solvers[0].outB=NULL;
		solvers[0].outX=NULL;
		if(ii==_vCycles-1)	solvers[_dCount-1].outX=out;
		else				solvers[_dCount-1].outX=in;
		solvers[0].inX=X;

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		solvers[_dCount-1].InitProlongation();
		solvers[0].SetProlongation();
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii<_vCycles-1)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
		solvers[0].SolveProlongation();
		SocketedStreamingSolver<Channels>::server.Reset();

		printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}
		if(_verbose)
			for(int i=0;i<_dCount;i++)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		solvers[0].UnSetProlongation();
	}
	delete[] solvers;
	delete B;
	delete X;
	B = X = NULL;
	delete[] solverInfo;
#endif
}
template< int Channels, class StorageType , class PixelType , class LabelType >
void SocketedMultigridClient< Channels , StorageType , PixelType , LabelType >::Run(StreamingGrid* pixels,StreamingGrid* out)
{
#if 0
	double t;
	DotProductStencil stencils[4];
	SolverInfo<Channels>* solverInfo = new SolverInfo<Channels>[_dCount];
	SocketedStreamingLaplacian<Channels,PixelType,StorageType>* sLaplacian = new SocketedStreamingLaplacian<Channels,PixelType,StorageType>();
	SocketedMultiGridStreamingSolver<Channels,StorageType>* solvers;
	SocketedStreamingSolver<Channels>::server.Reset();

	Vector<float> lowB,lowX;
	StreamingGrid *B,*X;
	// Also, the template should be StorageType, not Real
	MultiStreamIOClient* in;

	if(_vCycles>1)	in=new MultiStreamIOClient( size()*sizeof(float)*Channels , _height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true);
	else			in=NULL;
 
	solvers = new SocketedMultiGridStreamingSolver<Channels,StorageType>[_dCount];
	solvers[_dCount-1].showProgress=_progress;
	for(int i=1;i<_dCount;i++)		solvers[i].parent = &solvers[i-1];
	for(int i=0;i<_dCount-1;i++)	solvers[i].rChild =  solvers[i].pChild = &solvers[i+1];
	for(int i=0;i<_dCount;i++)		solvers[i].laneNum = _lanes;


	sLaplacian->parent = & solvers[_dCount-1];
	solvers[_dCount-1].rChild = sLaplacian;
	sLaplacian->Initialize(pixels,_iWeight,_gWeight,_start,_end,_width,_height,_iters,_leftSocket,_syncSockets,_rightSocket,true,_spherical);

	lowX.Resize(solvers[0].size()*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].size()*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
	for(int ii=0;ii<_vCycles;ii++)
	{
		if(!ii)
			pixels->SetServer(&SocketedStreamingSolver<Channels>::server);
		/////////////////
		// RESTRICTION //
		/////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=true;
		}
		if(ii)	solvers[_dCount-1].inX=in;
		else	solvers[_dCount-1].inX=NULL;
		solvers[_dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		solvers[0].outX=X;
		solvers[_dCount-1].rChild=sLaplacian;
		// Data for the interleaved streaming multigrid
		t=Time();
		// Initialize
		if(sLaplacian)	sLaplacian->InitRestriction()	,	sLaplacian->SetRestriction();
		else	solvers[_dCount-1].InitRestriction()	,	solvers[_dCount-1].SetRestriction();
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii)	in->SetServer(&SocketedStreamingSolver<Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
		if(sLaplacian)	sLaplacian->SolveRestriction();
		else			solvers[_dCount-1].SolveRestriction();
		SocketedStreamingSolver<Channels>::server.Reset();


		printf("Out-Of-Core Restriction:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}

		if(_verbose)
			for(int i=_dCount-1;i>=0;i--)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		if(ii==_vCycles-1 && in)	delete in,	in=NULL;
		if(sLaplacian)
		{
			SendOnSocket ( _serverSocket , sLaplacian->average , sizeof(sLaplacian->average) );
			sLaplacian->UnSetRestriction();
			delete sLaplacian;
			sLaplacian=NULL;
		}
		else	solvers[_dCount-1].UnSetRestriction();

		t=Time();
		// Send the low res stuff over to the server so that it can perform the in-core solve.
		SendOnSocket (_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		SendOnSocket (_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;

		stencils[0] = solvers[0].dotMajor;
		stencils[1] = solvers[0].d2DotMajor;
		stencils[2] = solvers[0].dotMinor;
		stencils[3] = solvers[0].d2DotMinor;
		SendOnSocket ( _serverSocket , stencils , sizeof( DotProductStencil ) * 4 );

		// Get back the in-core solution.
		ReceiveOnSocket	(_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		ReceiveOnSocket	(_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;
		printf("In-Core:    %f\n",Time()-t);



		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=_verbose;
		}
		solvers[_dCount-1].inX=NULL;
		solvers[0].outB=NULL;
		solvers[0].outX=NULL;
		if(ii==_vCycles-1)	solvers[_dCount-1].outX=out;
		else				solvers[_dCount-1].outX=in;
		solvers[0].inX=X;

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		solvers[_dCount-1].InitProlongation();
		solvers[0].SetProlongation();
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii<_vCycles-1)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
		solvers[0].SolveProlongation();
		SocketedStreamingSolver<Channels>::server.Reset();

		printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}
		if(_verbose)
			for(int i=0;i<_dCount;i++)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		solvers[0].UnSetProlongation();
	}
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;
	delete[] solverInfo;
#endif
}
