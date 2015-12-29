#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <Util/XPlatform.h>
#include <Util/Time.h>
#include <Util/CmdLineParser.h>
#include <Util/Half/half.h>
#include <LaplacianMatrix/SocketedMultigrid/MultigridServer.h>
#include <LaplacianMatrix/SocketedMultigrid/SocketData.h>


// [WARNING] This definition needs to be synchronized with ClientSocket
typedef uint16_t LType;

#define DEFAULT_BUFLEN 1024

cmdLineInt Count("count") , Port( "port" , 0 ) , Iters( "iters" , 5 ) , InCoreRes( "inCoreRes" , 1024 ) , MinMGRes( "minMGRes" , 64 ) , VCycles( "vCycles" , 1 ) , Quality( "quality" , 100 ) , Lanes( "lanes" , 1 );
cmdLineInt MinBandSize( "minBandSize" , 0 ) , TileWidth( "tileWidth" , 8192 ) , TileHeight( "tileHeight" , 8192 ) , UnknownType( "unknownType" , UNKNOWN_BLACK );
cmdLineReadable Verbose("verbose") , Spherical( "spherical" ) , Cylindrical( "cylindrical" ) , Progress( "progress" ) , GammaCorrection( "gCorrection" );
cmdLineReadable NoConjugateGradient( "noCG" ) , Lump( "lump" );
cmdLineReadable NoDisk( "noDisk" ) , NoNetwork( "noNetwork" ) , ShortSync( "shortSync" );
cmdLineString Prefix( "prefix" ) , TileExtension( "tileExt" );
cmdLineFloat IWeight( "iWeight" , 0 ) , GWeight( "gWeight" , 1 ) , GScale( "gScale" , 1 );
cmdLineReadable Deramp( "deramp" );
cmdLineReadable GrayImage( "gray" );
cmdLineReadable* params[]=
{
	&Port , &Count , &Iters , &InCoreRes , &MinMGRes , &VCycles , &Verbose , &Spherical , &Cylindrical , &Progress , &Quality , &Lanes , &NoConjugateGradient , &MinBandSize , &Prefix , &UnknownType , &GammaCorrection ,
	&IWeight , &GScale , &GWeight , &Lump ,
	&NoDisk , &NoNetwork , &ShortSync , &TileWidth , &TileHeight , &TileExtension , &Deramp ,
	&GrayImage ,
};

void ShowUsage( char* ex )
{
	printf( "Usage %s:\n" , ex ) , fflush( stdout );
	printf( "\t--%s <client count>\n" , Count.name ) , fflush( stdout );
	printf( "\t[--%s <preferred address prefix>]\n" , Prefix.name ) , fflush( stdout );
	printf( "\t[--%s <listen port>=%d]\n" , Port.name , Port.value ) , fflush( stdout );
	printf( "\t[--%s <minimum multigrid resolution>=%d]\n" , MinMGRes.name , MinMGRes.value ) , fflush( stdout );
	printf( "\t[--%s <Gauss-Seidel iterations>=%d]\n" , Iters.name , Iters.value ) , fflush( stdout );
	printf( "\t[--%s <in core resolution>=%d]\n" , InCoreRes.name , InCoreRes.value ) , fflush( stdout );
	printf( "\t[--%s <v-cycles>=%d]\n" , VCycles.name , VCycles.value ) , fflush( stdout );
	printf( "\t[--%s <image quality>=%d]\n" , Quality.name , Quality.value ) , fflush( stdout );
	printf( "\t[--%s <minimum band size>=%d]\n" , MinBandSize.name , MinBandSize.value ) , fflush( stdout );

	printf( "\t[--%s <default output tile width>=%d]\n" , TileWidth.name , TileWidth.value ) , fflush( stdout );
	printf( "\t[--%s <default output tile height>=%d]\n" , TileHeight.name , TileHeight.value ) , fflush( stdout );
	printf( "\t[--%s <default output file extension>]\n" , TileExtension.name ) , fflush( stdout );

	printf("\t[--%s <solver lanes>=%d]\n",Lanes.name,Lanes.value) , fflush( stdout );

	printf( "\t[--%s <unknown label behavior>=%d]\n" , UnknownType.name , UnknownType.value );
	printf( "\t\t%d] Do nothing\n" , UNKNOWN_NONE );
	printf( "\t\t%d] Fill with black\n" , UNKNOWN_BLACK );
	printf( "\t\t%d] Fill with harmonic\n" , UNKNOWN_HARMONIC );

	printf( "\t[--%s]\n" , Spherical.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , Cylindrical.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , Verbose.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , Progress.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , NoDisk.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , NoNetwork.name ) , fflush( stdout );
	printf( "\t[--%s <value interpolation weight>=%f]\n" , IWeight.name , IWeight.value ) , fflush( stdout );
	printf( "\t[--%s <gradient interpolation weight>=%f]\n" , GWeight.name , GWeight.value ) , fflush( stdout );
	printf( "\t[--%s <gradient scale>=%f]\n" , GScale.name , GScale.value ) , fflush( stdout );
	printf( "\t[--%s]\n" , GammaCorrection.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , Deramp.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , GrayImage.name ) , fflush( stdout );
	printf( "\t[--%s]\n" , Lump.name ) , fflush( stdout );
}
template< int PixelChannels >
int Execute( void )
{
	IOServer::Load();
	if( Verbose.set ) PrintHostAddresses();

#if MISHA_CODE_CLEAN_UP
	double t=Time();
	if( ShortSync.set )
		SocketedMultigridServer< PixelChannels , 3 , half  , LType > server
		(
			Prefix.value , Port.value , Count.value , Iters.value , InCoreRes.value , MinMGRes.value , VCycles.value , MinBandSize.value , TileWidth.value , TileHeight.value , TileExtension.set ? TileExtension.value : NULL , GammaCorrection.set , Quality.value , Lanes.value , Verbose.set , Spherical.set ? SPHERICAL_PERIODIC : (Cylindrical.set ? CYLINDRICAL_PERIODIC : NO_PERIODIC ) , IWeight.value , Lump.set , GWeight.value , GScale.value ,
			Deramp.set , UnknownType.value , Progress.set ,
			NoConjugateGradient.set , ShortSync.set
		);
	else
		SocketedMultigridServer< PixelChannels , 3 , float , LType > server
		(
			Prefix.value , Port.value , Count.value , Iters.value , InCoreRes.value , MinMGRes.value , VCycles.value , MinBandSize.value , TileWidth.value , TileHeight.value , TileExtension.set ? TileExtension.value : NULL , GammaCorrection.set , Quality.value , Lanes.value , Verbose.set , Spherical.set ? SPHERICAL_PERIODIC : (Cylindrical.set ? CYLINDRICAL_PERIODIC : NO_PERIODIC ) , IWeight.value , Lump.set , GWeight.value , GScale.value ,
			Deramp.set , UnknownType.value , Progress.set ,
			NoConjugateGradient.set , ShortSync.set
		);
#else // !MISHA_CODE_CLEAN_UP
	SocketedSuperMultigridServer server;
	if( !server.SetUp( Prefix.value , Port.value , Count.value , Iters.value , InCoreRes.value , MinMGRes.value , VCycles.value , MinBandSize.value , TileWidth.value , TileHeight.value , TileExtension.set ? TileExtension.value : NULL , GammaCorrection.set , Quality.value , Lanes.value , Verbose.set , Spherical.set ? SPHERICAL_PERIODIC : (Cylindrical.set ? CYLINDRICAL_PERIODIC : NO_PERIODIC ) , IWeight.value , Lump.set , GWeight.value , GScale.value ,
		Deramp.set , UnknownType.value , Progress.set ,
		NoConjugateGradient.set , ShortSync.set
		) )
		fprintf( stderr , "SocketedSuperMultigridServer::SetUp failed\n" );

	double t=Time();
	if( ShortSync.set ) server.Run< PixelChannels , 3 , half  , LType >();
	else                server.Run< PixelChannels , 3 , float , LType >();
#endif // MISHA_CODE_CLEAN_UP
	size_t current , peak;
	WorkingSetInfo( current , peak );
	if( Verbose.set )
	{
		printf( "    Running Time: %f\n",Time()-t) , fflush( stdout );
		printf( "Peak working set: %llu MB\n" , (unsigned long long)(peak>>20) ) , fflush( stdout );
	}

	IOServer::UnLoad();

	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{
	int paramNum = sizeof(params)/sizeof(cmdLineReadable*);
	if( !cmdLineParse( argc-1 , &argv[1] , paramNum , params , true ) ) return EXIT_FAILURE;
	if( !Count.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	if( GrayImage.set ) return Execute< 1 >( );
	else                return Execute< 3 >( );
}
