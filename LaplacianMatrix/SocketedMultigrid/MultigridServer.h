#ifndef MULTIGRID_SERVER_INCLUDED
#define MULTIGRID_SERVER_INCLUDED

#define MISHA_CODE_CLEAN_UP 0
enum
{
	NO_PERIODIC,
	CYLINDRICAL_PERIODIC,
	SPHERICAL_PERIODIC
};
#include "Util/XPlatform.h"
#include "LaplacianMatrix/SocketedStreamingSolver.h"
#include "MultigridThread.h"
#include "SocketData.h"

class ClientSocket
{
public:
	struct in_addr address;
	int port;
	Socket socket;
	AcceptorSocket listener;
	ClientData clientData;

	ClientSocket(void)
	{
		socket = _INVALID_SOCKET_;
		listener = _INVALID_ACCEPTOR_SOCKET_;
		port = 0;
	}
	~ClientSocket(void)
	{
		if( socket!=_INVALID_SOCKET_ ) CloseSocket( socket );
		if( listener!=_INVALID_ACCEPTOR_SOCKET_ ) CloseAcceptorSocket( listener );
	}
	static int Sort(const void* v1,const void* v2)
	{
		ClientSocket *s1=(ClientSocket*)v1;
		ClientSocket *s2=(ClientSocket*)v2;
		return s1->clientData.index-s2->clientData.index;
	}
};
class SphericalSynchronizer
{
	int _width,_cCount;
	ClientSocket* _clientSockets;
	int _rPasses,_pPasses,_vCycles;
public:

	SphericalSynchronizer( void );
	~SphericalSynchronizer( void );
	void init( int width , const int* widths , int cCount , int rPasses , int pPasses , int vCycles );
	int clients(void)		const;
	int port(int cIndex)	const;
	template<class DataType>	void Run(void);

	template< class DataType > static int RunThread( void* );
};

class ProcessSocket
{
public:
	IPData childAddressX , childAddressB , leftAddress;
	Socket sock;
	int index , offset;

	ProcessSocket( void ){ sock = _INVALID_SOCKET_; }
	~ProcessSocket( void ){ CloseSocket( sock ); }
	static int Sort( const void* v1 , const void* v2 ){ return ( (ProcessSocket*)v1 )->index - ( (ProcessSocket*)v2 )->index; }
};

#if MISHA_CODE_CLEAN_UP
template< int PixelChannels , int LabelChannels , class SyncType , class LabelType >
class SocketedMultigridServer
{
	static void _SolveInCore( DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
		std::vector< float >& in , std::vector< float >& out , double average[PixelChannels] );
public:
	static void Run
	(
		char* prefix , int port , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles , int minBandSize ,
		int tileWidth , int tileHeight , const char* tileExt ,
		bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType ,
		double iWeight , bool lump , double gWeight , double gScale ,
		bool removeAverageGradient , int unknownType , bool showProgress ,
		bool noCG , bool shortSync
	);
}
#else // !MISHA_CODE_CLEAN_UP
template< int PixelChannels , int LabelChannels , class SyncType >
class SocketedMultigridServer
{
	bool _noCG;
	int _oocLevels , _icLevels , _inCoreRes , _minMGRes;
	GlobalData _globalData;
	ProcessPartitionHierarchy< ProcessSocket > processPartition;
	ProcessConnectionInfo _inCoreConnectInfo;
	double _average[PixelChannels];
	std::vector< SphericalSynchronizer* > _oocSynchronizers;
	SphericalSynchronizer *_icSynchronizers , imageSynchronizer;
	Socket* _icSyncSockets;
	ThreadHandle* _synchronizerHandles;


	void SolveInCore(
		DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
		std::vector< float >& in , std::vector< float >& out , double average[PixelChannels] );

public:
	SocketedMultigridServer( void );
	~SocketedMultigridServer( void );
	template< class LabelType >
	bool SetUp( char* address , int port , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles , int minBandSize ,
		int tileWidth , int tileHeight , const char* tileExt , 
		bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType , double iWeight , bool lump , double gWeight , double gScale ,
		bool removeAverageGradient , int unknownType , bool showProgress ,
		bool noCG , bool shortSync , int* clipDimensions );

	template< class LabelType >
	bool SetUp( char* address , AcceptorSocket listenSocket , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles , int minBandSize ,
		int tileWidth , int tileHeight , const char* tileExt , 
		bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType , double iWeight , bool lump , double gWeight , double gScale ,
		bool removeAverageGradient , int unknownType , bool showProgress ,
		bool noCG , bool shortSync , int* clipDimensions );

	void Run( void );
};

class SocketedSuperMultigridServer
{
	GlobalData _globalData;
	bool _removeAverageGradient;
	bool _noCG;
	int _port , _clientCount , _inCoreRes , _minMGRes , _minBandSize;
	char _address[512];

public:
	// Set the parameters of the sockected super server.
	bool SetUp( char* prefix , int port , int clientCount , int iters , int inCoreRes , int minMGRes , int vCycles , int minBandSize ,
		int tileWidth , int tileHeight , const char* tileExt ,
		bool gammaCorrection , int quality , int lanes , bool verbose , int periodicType ,
		double iWeight , bool lump , double gWeight , double gScale ,
		bool removeAverageGradient ,
		int unknownType , bool showProgress ,
		bool noCG , bool shortSync );

	// Connect to the clients, partition the workload, and then set-up and run the socketed server.
	template< int PixelChannels , int LabelChannels , class SyncType , class LabelType > bool Run( void );
};
#endif // MISHA_CODE_CLEAN_UP
#include "MultigridServer.inl"
#endif // MULTIGRID_SERVER_INCLUDED