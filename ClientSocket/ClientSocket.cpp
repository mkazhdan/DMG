#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <Util/XPlatform.h>
#include <Util/Time.h>
#include <Util/CmdLineParser.h>
#include <Util/Half/half.h>
#include <LaplacianMatrix/SocketedMultigrid/MultigridServer.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif // !_WIN32
#include "LaplacianMatrix/SocketedMultigrid/MultigridProcess.h"
#include "LaplacianMatrix/SocketedMultigrid/MultigridServer.h"
#include "LaplacianMatrix/SocketedMultigrid/SocketData.h"
#include "Util/ImageStream.h"




#define DEFAULT_BUFLEN 1024
#define USE_SERVER 0	// [Question] Why does the server no longer work?

typedef uint16_t LType;

cmdLineInt Port( "port" , 0 );
cmdLineString Address( "address" , "127.0.0.1" ) , Prefix( "prefix" ) , TempDir( "temp" , "." );
cmdLineString Storage( "storage" );
cmdLineString LowPixels( "lowPixels" ) , Pixels( "pixels" ) , Labels( "labels" ) , Out( "out" ) , OutGuess( "outGuess" );
cmdLineInt Threads( "threads" , 1 ) , IOThreads( "ioThreads" , 1 );
cmdLineReadable InCore( "inCore" );
cmdLineInt Index( "index" , 0 );
cmdLineReadable* params[] =
{
	&LowPixels , &Pixels , &OutGuess , &Labels , &Out , &Port , &Address , &Prefix , &Index , &TempDir , &Storage , 
	&Threads , &IOThreads , &InCore ,
};
void ShowUsage( char* ex )
{
	printf( "Usage %s:\n",ex) , fflush( stdout );
	printf( "\t --%s <low frequency image pixels>\n" , LowPixels.name ) , fflush( stdout );
	printf( "\t --%s <composited image pixels>\n" , Pixels.name ) , fflush( stdout );
	printf( "\t --%s <image labels>\n" , Labels.name ) , fflush( stdout );
	printf( "\t --%s <output image>\n" , Out.name ) , fflush( stdout );
	printf( "\t --%s <server port>\n" , Port.name ) , fflush( stdout );
	printf( "\t[--%s <server connection address>=%s]\n" , Address.name , Address.value ) , fflush( stdout );
	printf( "\t[--%s <preferred address prefix>]\n" , Prefix.name ) , fflush( stdout );
	printf( "\t[--%s <client index>=%d]\n" , Index.name , Index.value ) , fflush( stdout );
	printf( "\t[--%s <scratch directory>=%s]\n" , TempDir.name , TempDir.value ) , fflush( stdout );
	printf( "\t[--%s <temporary storage type>]\n" , Storage.name ) , fflush( stdout );
	printf( "\t[--%s <number of computation threads>=%d]\n" , Threads.name , Threads.value ) , fflush( stdout );
	printf( "\t[--%s <number of I/O threads>=%d]\n" , IOThreads.name , IOThreads.value ) , fflush( stdout );
	printf( "\t[--%s]\n" , InCore.name ) , fflush( stdout );
}

template< class PixelType , int PixelChannels , class StorageType >
int _Execute( SocketedSuperMultigridClient< StorageType , PixelType , LType >& client , bool getGradientAverage , int w1 , int h1 , MultiStreamIOServer* ioServerPtr )
{
	int w2 , h2;
	std::vector< std::pair< LabelData< LType , 3 > , GradientAverage< PixelChannels > > > averageMap;
	if( getGradientAverage )
	{
		printf( "Computing gradient average(s)\n" );
		StreamingGrid *pixels = NULL , *labels = NULL;
		pixels = GetReadStream< PixelType , PixelChannels >( Pixels.value , w1 , h1 , client.gammaCorrection() , ioServerPtr , IOThreads.value );
		if( Labels.set ) labels = GetReadStream< LType , 3 >( Labels.value , w2 , h2 , false , ioServerPtr , IOThreads.value );
		SetGradientAverageMap< PixelChannels , 3 , PixelType , LType >( pixels , labels , averageMap , client.periodicType()==SPHERICAL_PERIODIC , client.showProgress() );

		if( pixels ) delete  pixels;
		if( labels ) delete labels;
	}
	StreamingGrid *lowPixels = NULL , *outGuess = NULL , *pixels = NULL , *labels = NULL , *out = NULL;
	pixels = GetReadStream< PixelType , PixelChannels >( Pixels.value , w1 , h1 , client.gammaCorrection() , ioServerPtr , IOThreads.value );
	if( LowPixels.set )
	{
		lowPixels = GetReadStream< PixelType , PixelChannels >( LowPixels.value , w2 , h2 , client.gammaCorrection() , ioServerPtr , IOThreads.value );
		if( w1!=w2 || h1!=h2 )
		{
			fprintf( stderr , "Image and low-image dimensions differ: %d x %d != %d x %d\n",w1,h1,w2,h2);
			return EXIT_FAILURE;
		}
	}
	if( OutGuess.set )
	{
		outGuess = GetReadStream< PixelType , PixelChannels >( OutGuess.value , w2 , h2 , client.gammaCorrection() , ioServerPtr , IOThreads.value );
		if( w1!=w2 || h1!=h2 )
		{
			fprintf( stderr , "Image and initial guess image dimensions differ: %d x %d != %d x %d\n" , w1 , h1 , w2 , h2 );
			return EXIT_FAILURE;
		}
	}
	if( Labels.set )
	{
		labels = GetReadStream< LType , 3 >( Labels.value , w2 , h2 , false , ioServerPtr , IOThreads.value );
		if( w1!=w2 || h1!=h2 )
		{
			fprintf(stderr,"Image and label dimensions differ: %d x %d != %d x %d\n",w1,h1,w2,h2);
			return EXIT_FAILURE;
		}
	}
	DefaultOutputTileWidth = client.tileWidth( );
	DefaultOutputTileHeight = client.tileHeight( );
	if( client.tileExtension() ) DefaultOutputTileExtension = client.tileExtension( );

	out = GetWriteStream< PixelType , PixelChannels >( Out.value , w1 , h1 , client.gammaCorrection() , client.quality() , ioServerPtr , client.hdr() , IOThreads.value );

	char tmpDir[1024];
	{
		sprintf( tmpDir , "TMP=%s" , TempDir.value );
		size_t l = strlen(tmpDir);
		if( tmpDir[l-1]!=FileSeparator ) tmpDir[l] = FileSeparator , tmpDir[l+1] = 0;
	}
#ifdef _WIN32
	_putenv( tmpDir );
#else // !_WIN32
	putenv( tmpDir );
#endif // _WIN32
	double t=Time();
	client.template Run< PixelChannels , 3 >( lowPixels , pixels , labels , outGuess , out , averageMap );
	size_t current , peak;
	WorkingSetInfo( current , peak );
	if( client.verbose() )
	{
		printf( "    Running Time: %f\n" , Time()-t ) , fflush( stdout );
		printf( "Peak working set: %llu MB\n" , (unsigned long long)(peak>>20) ) , fflush( stdout );
	}

	if( lowPixels ) delete lowPixels;
	if(    pixels ) delete    pixels;
	if(    labels ) delete    labels;
	if(  outGuess ) delete  outGuess;
	if(       out ) delete       out;
	IOServer::UnLoad();
	return EXIT_SUCCESS;
}
template< class PixelType , class StorageType >
int Execute( void )
{
#if USE_SERVER
	MultiStreamIOServer ioServer;
	MultiStreamIOServer* ioServerPtr = &ioServer;
#else // !USE_SERVER
	MultiStreamIOServer* ioServerPtr = NULL;
#endif // !USE_SERVER

	int w1 , h1;
	GetReadSize( Pixels.value , w1 , h1 );
	SocketedSuperMultigridClient< StorageType , PixelType , LType > client;
	bool getGradientAverage;
	bool useGrayImage;
	if( !client.SetUp( Address.value , Prefix.value , Port.value , Index.value , w1 , h1 , Threads.value , InCore.set , ioServerPtr , getGradientAverage , useGrayImage ) )
	{
		fprintf( stderr , "SocketedSuperMultigridClient::SetUp failed\n" );
		return EXIT_FAILURE;
	}
	if( useGrayImage ) return _Execute< PixelType , 1 , StorageType >( client , getGradientAverage , w1 , h1 , ioServerPtr );
	else               return _Execute< PixelType , 3 , StorageType >( client , getGradientAverage , w1 , h1 , ioServerPtr );
}
int main( int argc , char* argv[])
{
#ifdef _WIN32
	_setmaxstdio( 1024 );
#else // !_WIN32
	struct rlimit rl;
	getrlimit( RLIMIT_NOFILE , &rl ); 
	rl.rlim_cur = 1024+3; 
	setrlimit( RLIMIT_NOFILE , &rl );
#endif // _WIN32
	IOServer::Load();
	int paramNum = sizeof( params ) / sizeof( cmdLineReadable* );
	if( !cmdLineParse( argc-1 , &argv[1] , paramNum , params , true ) ) return EXIT_FAILURE;

	if( !Pixels.set || !Out.set || ( !Address.set && !Port.set ) )
	{
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
	if( !Port.set )
	{
		int a1 , a2 , a3 , a4;
		if( sscanf( Address.value , " %d.%d.%d.%d:%d " , &a1 , &a2 , &a3 , &a4 , &Port.value ) != 5 )
		{
			ShowUsage(argv[0]);
			return EXIT_FAILURE;
		}
		sprintf( Address.value, "%d.%d.%d.%d" , a1 , a2 , a3 , a4 );
	}

	if( Storage.set )
		if		( !strcasecmp( Storage.value , "float" ) )	return Execute< float , float >	( );
		else if	( !strcasecmp( Storage.value , "double" ) )	return Execute< float , double >( );
		else if	( !strcasecmp( Storage.value , "half" ) )	return Execute< float , half >	( );
		else
		{
			fprintf( stderr , "Unrecognized storage type: %s\n" ,  Storage.value );
			return EXIT_FAILURE;
		}
	else return Execute< float , half >( );
}
