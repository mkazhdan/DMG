#ifndef SOCKET_DATA_INCLUDED
#define SOCKET_DATA_INCLUDED

#include "Util/Array.h"


class GlobalData
{
public:
	int periodicType;
	bool verbose , showProgress , shortSync , gammaCorrection , hdr , lump;
	int iters , unknownType , vCycles , quality , lanes , width , height , cWidth , cHeight , tileWidth , tileHeight;
	double iWeight , gWeight , gScale;
	char tileExt[16];
};
template<int Channels>
class SolverInfo
{
public:
	double bSquareNorm;
	double rSquareNorm;
	double xSquareNorm;
	double solutionSum[Channels];
};
class ProcessData
{
public:
	int index , offset , maxIndex , maxOffset;
	int start , stop , width , height;
	int startDepth , endDepth;
	int children;
	int depths( void ) { return endDepth - startDepth + 1; }
};
class ProcessSocketData
{
public:
	ProcessSocketData( void ) { left = childX = childB = _INVALID_SOCKET_; }
	union
	{
		Socket left;
		Socket right;
	};
	union
	{
		Socket childX;
		Socket parentX;
	};
	union
	{
		Socket childB;
		Socket parentB;
	};
};

class ClientData
{
public:
	int index;
	union
	{
		int start;
		int width;
	};
	union
	{
		int end;
		int height;
	};
	union
	{
		int cEnd;
		int subClients;
	};
};
class IPData
{
public:
	IPData( void ) { port = 0; }
	EndpointAddress address;
	int port;
};
template< class Data >
class IndexedSocket
{
public:
	Socket sock;
	int index;
	Data data;

	static int Sort(const void* v1,const void* v2)
	{
		IndexedSocket *s1 = (IndexedSocket*)v1;
		IndexedSocket *s2 = (IndexedSocket*)v2;
		return s1->index-s2->index;
	}
};


class ProcessConnectionInfo
{
public:
	AcceptorSocket listenLeft , listenChildX , listenChildB;
	DataStream *left, *right;
	Socket *childX , *childB , parentX , parentB;
	int children;
	ProcessConnectionInfo( void )
	{
		listenLeft = listenChildX = listenChildB = _INVALID_ACCEPTOR_SOCKET_;
		left = right = NULL;
		parentX = parentB = _INVALID_SOCKET_;
		childX = childB = NULL;
		children = 0;
	}
	~ProcessConnectionInfo( void )
	{
		if( childX ) delete[] childX , childX = NULL;
		if( childB ) delete[] childB , childB = NULL;
	}
	void init( int children )
	{
		if( children )
		{
			this->children = children;
			childX = new Socket[children];
			childB = new Socket[children];
			for( int i=0 ; i<children ; i++ ) childX[i] = childB[i] = _INVALID_SOCKET_;
		}
	}
};

class AcceptThreadData
{
public:
	int connectCount;
	AcceptorSocket listenSocket;
	Socket *connectSocket;
	bool success;
};
int AcceptThread( void* vparams )
{
	AcceptThreadData* params = (AcceptThreadData*)vparams;
	for( int i=0 ; i<params->connectCount ; i++ )
	{
		params->connectSocket[i] = AcceptSocket( params->listenSocket );
		if ( params->connectSocket[i] == _INVALID_SOCKET_ )
		{
			params->success = false;
			delete params;
			return 0;
		}
	}
	CloseAcceptorSocket( params->listenSocket );
	params->success = true;
	delete params;
	return 0;
}
ThreadHandle SpawnListenerSocket( AcceptorSocket listenSocket , Socket* connectSocket , int connectCount )
{
 	AcceptThreadData* params = new AcceptThreadData();	// Have to allocate otherwise memory may be de-allocated (while the thread is running)
	params->listenSocket = listenSocket;
	params->connectSocket = connectSocket;
	params->connectCount = connectCount;
	ThreadHandle acceptThread = RunThread( AcceptThread , params );
	if( !TestThreadHandle( acceptThread ) ) fprintf( stderr , "[ERRPR] CreateThread failed" ) , exit( 0 );
	return acceptThread;

}

template< class Real , int Channels >
SocketBackedGrid* GetSocketBackedGrid( Socket sock , int index , int width , int height )
{
	SendOnSocket( sock , ( ConstPointer( int ) )( GetPointer( index ) ) , sizeof(index) );
	SendOnSocket( sock , ( ConstPointer( int ) )( GetPointer( width ) ) , sizeof(width) );
	return new SocketBackedGrid( sock , width * Channels * sizeof(Real) , height );
}

template< class Real , int Channels >
MultiSocketBackedGrid< Real , Channels >* GetMultiSocketBackedGrid( Socket* sockets , int socketCount , int rows )
{
	IndexedSocket< int >* _sockets = new IndexedSocket< int >[socketCount];
	int *widths = new int[socketCount];

	for( int k=0 ; k<socketCount ; k++ )
	{
		_sockets[k].sock = sockets[k];
		ReceiveOnSocket( _sockets[k].sock , GetPointer( _sockets[k].index ) , sizeof(_sockets[k].index) );
		ReceiveOnSocket( _sockets[k].sock , GetPointer( _sockets[k].data  ) , sizeof(_sockets[k].data ) );
	}
	qsort( _sockets , socketCount , sizeof( IndexedSocket<int> ) , IndexedSocket<int>::Sort );
	for( int k=0 ; k<socketCount ; k++ )
	{
		sockets[k] = _sockets[k].sock;
		widths[k] = _sockets[k].data;
	}
	MultiSocketBackedGrid< Real , Channels >* grid = new MultiSocketBackedGrid< Real , Channels >( sockets , widths , socketCount , rows );
	delete[] _sockets;
	delete[] widths;
	return grid;
}
#endif // SOCKET_DATA_INCLUDED