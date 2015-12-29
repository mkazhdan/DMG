#ifndef MULTIGRID_CLIENT_INCLUDED
#define MULTIGRID_CLIENT_INCLUDED
#include "LaplacianMatrix/SocketedStreamingSolver.h"
#include "MultigridThread.h"
#include "SocketData.h"

template< int PixelChannels , int LabelChannels , class StorageType , class SyncType , class PixelType , class LabelType >
class SocketedMultigridClient
{
	EndpointAddress _myAddress;
	int _processCount;
	Pointer( ProcessingBlockData ) _threadBlockData;
	Pointer( ProcessConnectionInfo ) _pConnectInfo;
	ClientData _clientData;
	int _start,_end,_width,_height,_cEnd,_cHeight,_dCount;
	bool _inCore;
	GlobalData _globalData;
	std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > > _averageMap;
public:
	SocketedMultigridClient(void);
	~SocketedMultigridClient(void);
	bool SetUp( char* serverAddress , char* prefix , int port , int index , int w , int h , const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& averageMap , bool inCore=false ); 
	void Run( StreamingGrid* lowPixels , StreamingGrid* pixels , StreamingGrid* labels , StreamingGrid* outGuess , StreamingGrid* out , MultiStreamIOServer* ioServer , bool showProgress );
	int width	(void)	const;	// the width of the image
	int cWidth	(void)	const;	// the clipped width of the image
	int size	(void)	const;	// the width of the band
	int cSize	(void)	const;	// the clipped width of the band
	int height	(void)	const;	// the height of the image/band
	int cHeight	(void)	const;	// the clipped height of the image/band
	int start	(void)	const;
	int end		(void)	const;
	int cEnd	(void)	const;

	bool gammaCorrection( void ) const;
	int quality	(void)	const;
	int iters	(void)	const;
	int vCycles	(void)	const;
	bool shortSync( void ) const;
	int         tileWidth    ( void ) const;
	int         tileHeight   ( void ) const;
	const char* tileExtension( void ) const;
	int processCount( void ) const;
};

template< class StorageType , class PixelType , class LabelType >
class SocketedSuperMultigridClient
{
	char* _address;
	char* _prefix;
	int _port;
	int _clientNum;
	bool _inCore;
	int _blockSize;
	GlobalData _globalData;
	ClientData _clientData;

	template< int PixelChannels , int LabelChannels , class SyncType >
	class _SocketedMultigridClientData
	{
	public:
		_SocketedMultigridClientData( void ) { lowPixels = pixels = outGuess = labels = out = NULL ; averageMapPtr = NULL; }
		~_SocketedMultigridClientData( void )
		{
			if( lowPixels ) delete lowPixels;
			if( pixels ) delete pixels;
			if( labels ) delete labels;
			if( outGuess ) delete outGuess;
			if( out ) delete out;
			outGuess = lowPixels = pixels = labels = NULL;
		}
		SharedMemoryGrid *lowPixels , *pixels , *outGuess , *labels , *out;
		const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >* averageMapPtr;
		SocketedMultigridClient< PixelChannels , LabelChannels , StorageType , SyncType , PixelType , LabelType > client;
		MultiStreamIOServer* ioServer;
		bool showProgress;

		char *address , *prefix;
		int port , index , width , cWidth , height , cHeight;
		bool inCore;

		static int ProcessThread( void* vparams )
		{
			_SocketedMultigridClientData* data = (_SocketedMultigridClientData*)vparams;
			data->client.Run( data->lowPixels , data->pixels , data->labels , data->outGuess , data->out , data->ioServer , data->showProgress );
			return 0;
		}
		ThreadHandle spawnProcess( void ){ return RunThread( ProcessThread , this ); }
		static int SetUpThread( void* vparams )
		{
			_SocketedMultigridClientData* data = (_SocketedMultigridClientData*)vparams;
			data->client.SetUp( data->address , data->prefix , data->port , data->index , data->width , data->height , *(data->averageMapPtr) , data->inCore );
			return 0;
		}
		ThreadHandle spawnSetUp( void ){ return RunThread( SetUpThread , this ); }
	};

	template< int PixelChannels , int LabelChannels , class SyncType >
	void _Run( StreamingGrid* lowPixels , StreamingGrid* pixels , StreamingGrid* labels , StreamingGrid* outGuess , StreamingGrid* out , const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& averageMap );
	MultiStreamIOServer* sharedIOServer;
public:
	bool SetUp( char* address , char* prefix , int port , int index , int w , int h , int subClients , bool inCore , MultiStreamIOServer* ioServer , bool& getGradientAverage , bool& useGrayImage );
	template< int PixelChannels , int LabelChannels >
	void Run( StreamingGrid* lowPixels , StreamingGrid* pixels , StreamingGrid* labels , StreamingGrid* outGuess , StreamingGrid* out , const std::vector< std::pair< LabelData< LabelType , LabelChannels > , GradientAverage< PixelChannels > > >& averageMap );
	int size	(void)	const;	// the width of the band
	int cSize	(void)	const;	// the clipped width of the band
	int height	(void)	const;	// the height of the image/band
	int cHeight	(void)	const;	// the clipped height of the image/band
	bool gammaCorrection( void ) const;
	int quality	(void)	const;
	bool showProgress( void ) const;
	int         tileWidth    ( void ) const;
	int         tileHeight   ( void ) const;
	const char* tileExtension( void ) const;
	bool shortSync( void ) const;
	int periodicType( void ) const;
	bool verbose( void ) const;
};
#include "MultigridProcess.inl"
#endif // MULTIGRID_CLIENT_INCLUDED