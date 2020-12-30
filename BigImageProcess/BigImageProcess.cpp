#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "Util/Half/half.h"
#include "Util/MultiStreamIO.h"
#include "Util/CmdLineParser.h"
#include "Util/ImageStream.h"
#include "LinearAlgebra/Vector.h"

cmdLineString In( "in" ) , Out( "out" ) , TileExtension( "tileExt" );
cmdLineIntArray<4> Crop( "crop" );
cmdLineInt Quality( "quality" , 100 ) , Down( "down" ) , TileWidth( "tileWidth" , 8192 ) , TileHeight( "tileHeight" , 8192 ) , Threads( "threads" , 1 ) , IOThreads( "ioThreads" , 1 );
cmdLineFloatArray<3> DCTermOffset( "dcOffset" );
cmdLineReadable Progress( "progress" ) , HDR( "hdr" ) , GammaCorrection( "gCorrect" ) , Clamp( "clamp" ) , WhiteBalance( "whiteBalance" ) , WhiteBalanceOffset( "whiteBalanceOffset" ) , Gray( "gray" );
cmdLineFloat Gamma( "gamma" , 1.0 ) , Brightness( "brighten" , 1.0 ) , WhiteBalanceValue( "whiteBalanceValue" , 0.05f );
cmdLineReadable* params[]=
{
	&In , &Out , &Down , &Quality , &Crop , &DCTermOffset , &Gamma , &Brightness , &TileWidth , &TileHeight , &Progress , &TileExtension , &HDR , &Clamp , &Threads , &IOThreads , &WhiteBalance , &WhiteBalanceOffset , &WhiteBalanceValue , &Gray , 
};

void ShowUsage(char* ex)
{
	printf("Usage %s:\n",ex);
	printf( "\t--%s <input image>\n",In.name);
	printf( "\t[--%s <output image>]\n",Out.name);
	printf( "\t[--%s <2X down sampling passes>]\n",Down.name);
	printf( "\t[--%s <startX , startY , width , height >]\n" , Crop.name);
	printf( "\t[--%s <JPEG compression quality>=%d]\n",Quality.name,Quality.value);
	printf( "\t[--%s <default output tile width>=%d]\n" , TileWidth.name , TileWidth.value );
	printf( "\t[--%s <default output tile height>=%d]\n" , TileHeight.name , TileHeight.value );
	printf( "\t[--%s <default output file extension>]\n" , TileExtension.name );
	printf( "\t[--%s <offset for the DC term>]\n",DCTermOffset.name);
	printf( "\t[--%s <gamma correction term>=%f]\n",Gamma.name,Gamma.value);
	printf( "\t[--%s <white balance percent>=%.3f]\n" , WhiteBalanceValue.name , WhiteBalanceValue.value );
	printf( "\t[--%s <brightness multiplier>=%f]\n" , Brightness.name , Brightness.value );
	printf( "\t\tPixel <- [ (Pixel + Offset) * Brightness ]^Gamma\n");
	printf( "\t[--%s <number of threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <number of I/O threads>=%d]\n" , IOThreads.name , IOThreads.value );
	printf( "\t[--%s]\n" , WhiteBalance.name );
	printf( "\t[--%s]\n" , WhiteBalanceOffset.name );
	printf( "\t[--%s]\n" , Gray.name );
	printf( "\t[--%s]\n" , HDR.name );
	printf( "\t[--%s]\n" , Clamp.name );
	printf( "\t[--%s]\n" , Progress.name );
}
template< class PixelType >
class MyImage
{
public:
	StreamingGrid* data;
	int width , height;

	MyImage(void)
	{
		data=NULL;
		width=height=0;
	}
	~MyImage(void)
	{
		if(data)	delete data;
		data=NULL;
		width=height=0;
	}
};

double DownStencil2[2][2]=
{
	{1,1},
	{1,1}
};
double DownStencil4[4][4]=
{
	{1,3,3,1},
	{3,9,9,3},
	{3,9,9,3},
	{1,3,3,1}
};
double DownStencil6[6][6]=
{
	{  1,  5, 10, 10,  5,  1},
	{  5, 25, 50, 50, 25,  5},
	{ 10, 50,100,100, 50, 10},
	{ 10, 50,100,100, 50, 10},
	{  5, 25, 50, 50, 25,  5},
	{  1,  5, 10, 10,  5,  1}
};
template< class Real , int Channels >
void SetRange( MyImage< Color< Real , Channels > >& img , Color< Real , Channels >& min , Color< Real , Channels >& max , double fraction , int bins=(1<<20) )
{
	img.data->reset( true , 1 );
	double* histograms[Channels];
	for( int c=0 ; c<Channels ; c++ ) histograms[c] = new double[bins] , memset( histograms[c] , 0 , sizeof(double) * bins );
	for( int j=0 ; j<img.height ; j++ )
	{
		if( Progress.set )
		{
			double r = double(j) / img.height * 100;
			printf("[%.1f%%]\r" , r );
		}
		ConstPointer( Color< Real , Channels > ) row = ( ConstPointer( Color< Real , Channels > ) )(*img.data)[j];
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<img.width ; i++ )
		{
			Color< Real , Channels > tmp = row[i];
			for( int c=0 ; c<Channels ; c++ )
			{
				double x = (double)tmp[c] * (bins-1);
				int ix1 = (int)floor(x) , ix2 = ix1+1;
				double dx = x - ix1;
				ix1 = std::max< int >( 0 , std::min< int >( ix1 , bins-1 ) );
				ix2 = std::max< int >( 0 , std::min< int >( ix2 , bins-1 ) );
				histograms[c][ix1] += (1.-dx);
				histograms[c][ix2] +=     dx ;
			}
		}
		img.data->advance();
	}
	double cutoff = fraction;
	cutoff *= img.width , cutoff *= img.height;

	for( int c=0 ; c<Channels ; c++ )
	{
		double sum;
		
		sum = 0;
		for( int i=0 ; i<bins ; i++ )
		{
			sum += histograms[c][i];
			if( sum>cutoff )
			{
				min[c] = (Real)( (double)i / bins );
				break;
			}
		}
		sum = 0;
		for( int i=bins-1 ; i>=0 ; i-- )
		{
			sum += histograms[c][i];
			if( sum>cutoff )
			{
				max[c] = (Real)( (double)i / bins );
				break;
			}
		}
	}
	for( int c=0 ; c<Channels ; c++ ) delete[] histograms[c];
}
template< class PixelType , class Real >
void DownSample2( const MyImage< PixelType >& high , MyImage< PixelType >& low , PixelType dcOffset , PixelType brightness , double gamma=1.0 )
{
	high.data->reset(true,1);
	low.data->reset(false,1);

	PixelType *_highRows[2] , *highRows[2];

	for( int i=0 ; i<2 ; i++ ) _highRows[i] = new PixelType[ high.width ];

	for( int j=0 ; j<low.height ; j++ )
	{
		if( Progress.set )
		{
			double r = double(j) / low.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		Pointer( PixelType ) lowRow = ( Pointer( PixelType ) )(*low.data)[j];
		// Copy in the next two high rows.
		for(int l=0;l<2;l++)
		{
			int jj=2*j+l;
			if(jj>=0 && jj<high.height)
			{
				Pointer( PixelType ) highRow = ( Pointer( PixelType ) )(*high.data)[jj];
#pragma omp parallel for num_threads( Threads.value )
				for( int k=0 ; k<high.width ; k++ )
				{
					highRow[k] = ( highRow[k] + dcOffset ) * brightness;
					highRow[k] = pow( highRow[k] , gamma );
				}
				memcpy( _highRows[jj%2] , highRow , sizeof(PixelType) * high.width );
			}
			high.data->advance();
		}
		for(int l=0;l<2;l++)
		{
			int jj=2*j+l;
			if(jj<0 || jj>=high.height)	highRows[l]=NULL;
			else						highRows[l]=_highRows[jj%2];
		}

#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<low.width ; i++ )
		{
			PixelType value;
			value *= 0;
			double sum = 0;
			for(int l=0;l<2;l++)
			{
				int jj=j*2+l;
				if(jj<0 || jj>=high.height)	continue;
				for(int k=0;k<2;k++)
				{
					int ii=2*i+k;
					if(ii<0 || ii>=high.width)	continue;
					value += highRows[l][ii]*(Real)DownStencil2[k][l];
					sum += DownStencil2[k][l];
				}
			}
			lowRow[i] = value/(Real)sum;
		}
		low.data->advance();
	}
	for(int i=0;i<2;i++)	delete[] _highRows[i];
}
template< class PixelType , class Real >
void DownSample2( const MyImage< PixelType >& high , MyImage< PixelType >& low )
{
	PixelType offset( (Real)0 ) , brightness( (Real)1 );
	DownSample2< PixelType , Real >( high , low , offset , brightness );
}

template< class PixelType , class Real >
void DownSample4( const MyImage< PixelType >& high , MyImage< PixelType >& low , PixelType dcOffset , PixelType brightness , double gamma=1.0 )
{
	high.data->reset(true,1);
	low.data->reset(false,1);

	PixelType *_highRows[4],*highRows[4];
	for(int i=0;i<4;i++)	_highRows[i]=new PixelType[high.width];

	// Buffer the first row
	memcpy(_highRows[0],(*high.data)[0],sizeof(PixelType)*high.width);
#pragma omp parallel for num_threads( Threads.value )
	for(int k=0;k<high.width;k++)
	{
		_highRows[0][k] = ( _highRows[0][k] + dcOffset ) * brightness;
		_highRows[0][k] = pow( _highRows[0][k] , gamma );
	}
	high.data->advance();

	for( int j=0 ; j<low.height ; j++ )
	{
		if( Progress.set )
		{
			double r = double(j) / low.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		PixelType* lowRow=(PixelType*)(*low.data)[j];
		// Copy in the next two high rows.
		for(int l=1;l<3;l++)
		{
			int jj=2*j+l;
			if(jj>=0 && jj<high.height)
			{
				memcpy(_highRows[jj%4],(PixelType*)(*high.data)[jj],sizeof(PixelType)*high.width);
#pragma omp parallel for num_threads( Threads.value )
				for(int k=0;k<high.width;k++)
				{
					_highRows[jj%4][k] = ( _highRows[jj%4][k] + dcOffset ) * brightness;
					_highRows[jj%4][k] = pow( _highRows[jj%4][k] , gamma );
				}
			}
			high.data->advance();
		}
		for(int l=-1;l<3;l++)
		{
			int jj=2*j+l;
			if(jj<0 || jj>=high.height)	highRows[l+1]=NULL;
			else						highRows[l+1]=_highRows[jj%4];
		}

#pragma omp parallel for num_threads( Threads.value )
		for(int i=0;i<low.width;i++)
		{
			PixelType value(0);
			double sum=0;
			for(int l=-1;l<3;l++)
			{
				int jj=j*2+l;
				if(jj<0 || jj>=high.height)	continue;
				for(int k=-1;k<3;k++)
				{
					int ii=2*i+k;
					if(ii<0 || ii>=high.width)	continue;
					 value += highRows[l+1][ii]*(Real)DownStencil4[k+1][l+1];
					sum += DownStencil4[k+1][l+1];
				}
			}
			lowRow[i] = value/(Real)sum;
		}
		low.data->advance();
	}
	for(int i=0;i<4;i++)	delete[] _highRows[i];
}
template< class PixelType , class Real >
void DownSample4( const MyImage< PixelType >& high , MyImage< PixelType >& low )
{
	PixelType offset(0) , brightness(1);
	DownSample4< PixelType , Real >( high , low , offset , brightness );
}

template< class PixelType , class Real >
void DownSample6( const MyImage< PixelType >& high , MyImage< PixelType >& low , PixelType dcOffset , PixelType brightness , double gamma=1.0 )
{
	high.data->reset(true,1);
	low.data->reset(false,1);

	PixelType *_highRows[6],*highRows[6];
	for(int i=0;i<6;i++)	_highRows[i] = new PixelType[high.width];

	// Buffer the first two row
	for(int j=0;j<2;j++)
	{
		PixelType* highRow=(PixelType*)(*high.data)[j];
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<high.width ; i++ )
		{
			_highRows[j][i] = ( highRow[i] + dcOffset ) * brightness;
			_highRows[j][i] = pow( _highRows[j][i] , gamma );
		}
		high.data->advance();
	}

	for(int j=0;j<low.height;j++)
	{
		if( Progress.set )
		{
			double r = double(j) / low.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		PixelType* lowRow=(PixelType*)(*low.data)[j];
		// Copy in the next two high rows.
		for(int l=2;l<4;l++)
		{
			int jj=2*j+l;
			if(jj>=0 && jj<high.height)
			{
				PixelType* highRow = (PixelType*)(*high.data)[jj];
				memcpy(_highRows[jj%6],highRow,sizeof(PixelType)*high.width);
#pragma omp parallel for num_threads( Threads.value )
				for( int k=0 ; k<high.width ; k++ )
				{
					_highRows[jj%6][k] = ( _highRows[jj%6][k] + dcOffset ) * brightness;
					_highRows[jj%6][k] = pow( _highRows[jj%6][k] , gamma );
				}
			}
			high.data->advance();
		}
		for(int l=-2;l<4;l++)
		{
			int jj=2*j+l;
			if(jj<0 || jj>=high.height)	highRows[l+2]=NULL;
			else						highRows[l+2]=_highRows[jj%6];
		}

#pragma omp parallel for num_threads( Threads.value )
		for(int i=0;i<low.width;i++)
		{
			PixelType value(0);
			double sum = 0;
			for(int l=-2;l<4;l++)
			{
				int jj=j*2+l;
				if(jj<0 || jj>=high.height)	continue;
				for(int k=-2;k<4;k++)
				{
					int ii=2*i+k;
					if(ii<0 || ii>=high.width)	continue;
					value += highRows[l+2][ii]*(Real)DownStencil6[k+2][l+2];
					sum += DownStencil6[k+2][l+2];
				}
			}
			lowRow[i] = value/(Real)sum;
		}
		low.data->advance();
	}
	for(int i=0;i<6;i++)	delete[] _highRows[i];
}
template< class PixelType , class Real>
void DownSample6( const MyImage< PixelType >& high , MyImage< PixelType >& low )
{
	PixelType offset(0) , brightness(1);
	DownSample6< PixelType , Real >( high , low , offset , brightness );
}

template< class PixelType , bool Remap >
void CropWindow( const MyImage< PixelType >& in , MyImage< PixelType >& out , int startX , int startY , PixelType dcOffset , PixelType brightness , double gamma=1.0 )
{
	for( int j=0 ; j<in.height ; j++ )
	{
		if( Progress.set )
		{
			double r = double(j) / in.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		Pointer( PixelType ) inRow = ( Pointer( PixelType ) )(*in.data)[j];
		if( j>=startY && j<startY+out.height )
		{
			Pointer( PixelType ) outRow = ( Pointer( PixelType ) )(*out.data)[j-startY];
			if( Remap )
#pragma omp parallel for num_threads( Threads.value )
				for( int i=0 ; i<out.width ; i++ )
				{
					outRow[i] = ( inRow[i+startX]+dcOffset ) * brightness;
					outRow[i] = pow( outRow[i] , gamma );
				}
			else memcpy( outRow , inRow + startX , sizeof( PixelType ) * out.width );
			out.data->advance();
		}
		in.data->advance();
	}
}
template< class PixelType , class Real >
void CropWindow( const MyImage< PixelType >& in , MyImage< PixelType >& out , int startX , int startY )
{
	PixelType offset( (Real)0) , brightness( (Real)1 );
	CropWindow< PixelType , false >( in , out , startX , startY , offset , brightness );
}
template< class Real , int Channels >
void ImageInfo( const MyImage< Color< Real , Channels > >& in )
{
	double lum , brightness=0 , variance=0;
	double average[Channels];
	for( int c=0 ; c<Channels ; c++ ) average[c] = 0;
	Color< double , Channels > min,max;

	for(int j=0;j<in.height;j++)
	{
		if( Progress.set )
		{
			double r = double(j) / in.height * 100;
			printf("[%.1f%%]\r" , r );
		}
		Pointer( Real ) inRow=( Pointer( Real ) )(*in.data)[j];
		for(int c=0;c<Channels;c++)
			for(int i=0;i<in.width;i++)
			{
				average[c]+=inRow[i+c*in.width];
				variance+=inRow[i+c*in.width]*inRow[i+c*in.width];
				brightness+=log(double(inRow[i+c*in.width]));

				if((!i && !j) || inRow[i+c*in.width]<min[c])	min[c]=inRow[i+c*in.width];
				if((!i && !j) || inRow[i+c*in.width]>max[c])	max[c]=inRow[i+c*in.width];
			}
		in.data->advance();
	}
	variance/=Channels;
	variance/=in.width;
	variance/=in.height;
	brightness/=Channels;
	brightness/=in.width;
	brightness/=in.height;
	for( int c=0 ; c<Channels ; c++ )
	{
		average[c]/=in.width;
		average[c]/=in.height;
	}
	lum = 0;
	for( int c=0 ; c<Channels ; c++ ) lum += average[c];
	lum /= Channels;
	variance-=lum*lum;
	printf( "Dimensions: %d x %d\n" , in.width , in.height );
	printf( "     Range:" ) ; for( int c=0 ; c<Channels ; c++ ) printf( " [%.5f,%.5f]" , min[c] , max[c] ) ; printf( "\n" );
	printf( "   Average:" ) ; for( int c=0 ; c<Channels ; c++ ) printf( " %f" , average[c] ) ; printf( "\n" );
	printf( " Luminance: %.5f\t%.5f\n",lum,sqrt(variance));
	printf( "Brightness: %.5f\n" , exp(brightness) );
}
template< class Real , int Channels >
int _run( void )
{
	Color< Real , Channels > brightness;
	for( int c=0 ; c<Channels ; c++ ) brightness[c] = (Real)1.;
	MultiStreamIOServer ioServer;
	char tmp[]  = "TMP=";
#ifdef _WIN32
	_putenv( tmp );
#else // !_WIN32
	putenv( tmp );
#endif // _WIN32

	if( Brightness.set ) for( int c=0 ; c<Channels ; c++ ) brightness[c] = (Real)Brightness.value;
	Color< Real , Channels > dcOffset;
	for( int c=0 ; c<Channels ; c++ ) dcOffset[c] = (Real)0.;
	if( DCTermOffset.set ) for( int c=0 ; c<Channels ; c++ ) dcOffset[c] = (Real)DCTermOffset.values[c];

	MyImage< Color< Real , Channels > > inImage , outImage;
	inImage.data = GetReadStream< Real , Channels >( In.value , inImage.width , inImage.height , GammaCorrection.set , NULL , IOThreads.value );
	if( !inImage.data )
	{
		fprintf(stderr,"Failed to read in image from: %s\n",In.value);
		return EXIT_FAILURE;
	}

	if( !Out.set ){ ImageInfo< Real , Channels >( inImage ) ; return EXIT_SUCCESS; }

	if( WhiteBalanceValue.set )
	{
		MyImage< Color< Real , Channels > > scratch;
		scratch.data = GetReadStream< Real , Channels >( In.value , scratch.width , scratch.height , GammaCorrection.set , NULL , IOThreads.value );
		Color< Real , Channels > min , max;
		SetRange< Real , Channels >( scratch , min , max , WhiteBalanceValue.value/100. );
		if( !WhiteBalanceOffset.set ) for( int c=0 ; c<Channels ; c++ ) min[c] = (Real)0.;
		// p_new = (p_old + off) * bright
		for( int c=0 ; c<Channels ; c++ )
		{
			dcOffset[c] = -min[c] , brightness[c] = (Real)( 1./(max[c]-min[c]) );
			printf( "%d] [ %g , %g ]\n" , c , (double)min[c] , (double)max[c] );
		}
	}
	bool cropSet = Crop.set;
	if( !Crop.set && (!Down.set || !Down.value) )
	{
		Crop.set=true;
		Down.set=false;
		Crop.values[0]=Crop.values[1]=0;
		Crop.values[2]=inImage.width;
		Crop.values[3]=inImage.height;
	}
	if( Down.set )
	{
		outImage.width=inImage.width;
		outImage.height=inImage.height;
		for(int i=0;i<Down.value;i++)
		{
			outImage.width=(outImage.width+1)/2;
			outImage.height=(outImage.height+1)/2;
		}
	}
	else
	{
		int startX,startY,cWidth,cHeight;
		startX  = Crop.values[0];
		startY  = Crop.values[1];
		cWidth  = Crop.values[2];
		cHeight = Crop.values[3];
		if( startX+cWidth >inImage.width  )	cWidth  = inImage.width  - startX;
		if( startY+cHeight>inImage.height ) cHeight = inImage.height - startY;
		if( cWidth<=0 || cHeight<=0 )
		{
			fprintf( stderr , "Window width and height must be positive: (%d , %d) %d x %d\n" , startX , startY , cWidth , cHeight );
			return EXIT_FAILURE;
		}
		outImage.width=cWidth;
		outImage.height=cHeight;
	}

	DefaultOutputTileWidth = TileWidth.value;
	DefaultOutputTileHeight = TileHeight.value;
	if( TileExtension.set ) DefaultOutputTileExtension = TileExtension.value;

	outImage.data = GetWriteStream< Real , Channels >( Out.value , outImage.width , outImage.height , false , Quality.value , NULL , HDR.set , IOThreads.value );
	if( !outImage.data )
	{
		fprintf( stderr , "Failed to write image to: %s\n" , Out.value );
		return EXIT_FAILURE;
	}

	if( Down.set && Down.value )
	{
		MyImage< Color< Real , Channels > > tempIn,tempOut;
		tempIn.data=inImage.data;
		tempIn.width=inImage.width;
		tempIn.height=inImage.height;
		inImage.data=NULL;
		for(int i=0;i<Down.value-1;i++)
		{
			tempOut.width=(tempIn.width+1)/2;
			tempOut.height=(tempIn.height+1)/2;
			tempOut.data=new MultiStreamIOClient(sizeof(Real)*Channels*tempOut.width,tempOut.height,2,NULL,true);
			if( Progress.set ) printf("Down-Sampling %d x %d -> %d x %d\n",tempIn.width,tempIn.height,tempOut.width,tempOut.height);
			if(!i)	DownSample2< Color< Real , Channels > , Real >( tempIn , tempOut , dcOffset , brightness , Gamma.value );
			else	DownSample2< Color< Real , Channels > , Real >( tempIn , tempOut );
			delete tempIn.data;
			tempIn.data=tempOut.data;
			tempIn.width=tempOut.width;
			tempIn.height=tempOut.height;
			tempOut.data=NULL;
		}
		if( Progress.set ) printf("Down-Sampling %d x %d -> %d x %d\n",tempIn.width,tempIn.height,outImage.width,outImage.height);
		if(Down.value==1)	DownSample2< Color< Real , Channels > , Real >( tempIn , outImage , dcOffset , brightness , Gamma.value );
		else				DownSample2< Color< Real , Channels > , Real >( tempIn , outImage );
	}
	else if( Crop.set )
	{
		if( cropSet ) printf( "Cropping: (%d , %d) x (%d , %d)\n",Crop.values[0],Crop.values[1],Crop.values[0]+outImage.width,Crop.values[1]+outImage.height );
		if( !DCTermOffset.set && !Gamma.set && !Brightness.set && !WhiteBalanceValue.set ) CropWindow< Color< Real , Channels > , Real >( inImage , outImage , Crop.values[0] , Crop.values[1] );
		else                                                                               CropWindow< Color< Real , Channels > , true >( inImage , outImage , Crop.values[0] , Crop.values[1] , dcOffset , brightness , Gamma.value );
	}
	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{
	int paramNum = sizeof(params)/sizeof(cmdLineReadable*);
	cmdLineParse( argc-1 , &argv[1] , paramNum , params , 0 );

	if( WhiteBalance.set ) WhiteBalanceValue.set = true;
	bool remapColors = Gray.set || DCTermOffset.set || Gamma.set || Brightness.set || Down.set || GammaCorrection.set || WhiteBalanceValue.set;
remapColors = true;
	if( !In.set )
	{
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
	if( Down.set && Crop.set )
	{
		fprintf( stderr , "Only one operation can be implemented at a time\n" );
		return EXIT_FAILURE;
	}
	if( WhiteBalanceValue.set && ( Brightness.set || DCTermOffset.set ) )
	{
		fprintf( stderr , "Only one operation can be implemented at a time\n" );
		return EXIT_FAILURE;
	}

	double t = Time();
	if( Gray.set )
		if( remapColors ) _run< float        , 1 >();
		else              _run< unsigned int , 1 >();
	else
		if( remapColors ) _run< float        , 3 >();
		else              _run< unsigned int , 3 >();
	if( Progress.set ) printf( "Image Process Time: %f\n" , Time() - t );
	return EXIT_SUCCESS;
}