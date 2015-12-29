/*
Copyright (c) 2008, Michael Kazhdan
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
#include "ChannelConverter.h"


//////////////////////
// WriteImageStream //
//////////////////////
template< class ChannelType , int Channels >
WriteImageStream< ChannelType , Channels >::WriteImageStream( ImageWriteInfo< ChannelType > iwInfo , char* fileName , int width , int height , bool gammaEncode , int quality , bool clamp , MultiStreamIOServer* ioServer ) 
{
	_InitWrite     = iwInfo.InitWrite;
	_WriteRow      = iwInfo.WriteRow;
	_FinalizeWrite = iwInfo.FinalizeWrite;
	current=0;
	_q=quality;
	_w=width;
	_h=height;
	_clamp=clamp;
	_gammaEncode = gammaEncode;
	info=NULL;
	_pixels = NullPointer< ChannelType >();
	_Init( fileName , ioServer );
}
template< class ChannelType , int Channels >
WriteImageStream< ChannelType , Channels >::WriteImageStream( void ) 
{
	_InitWrite		= NULL;
	_WriteRow		= NULL;
	_FinalizeWrite	= NULL;
	current = 0;
	_q = 100;
	_w = 0;
	_h = 0;
	_clamp = false;
	_gammaEncode = false;
	info = NULL;
	_pixels = NullPointer< ChannelType >();
}
template< class ChannelType , int Channels >
void WriteImageStream< ChannelType , Channels >::Init( ImageWriteInfo< ChannelType > iwInfo , char* fileName , int width , int height , bool gammaEncode , int quality , bool clamp , MultiStreamIOServer* ioServer ) 
{
	_InitWrite		= iwInfo.InitWrite;
	_WriteRow		= iwInfo.WriteRow;
	_FinalizeWrite	= iwInfo.FinalizeWrite;
	current=0;
	_q=quality;
	_w=width;
	_h=height;
	_clamp=clamp;
	_gammaEncode = gammaEncode;
	info=NULL;
	_pixels = NullPointer< ChannelType >();
	_Init( fileName , ioServer );
}

template< class ChannelType , int Channels >
void WriteImageStream< ChannelType , Channels >::_Init( char* fileName , MultiStreamIOServer* ioServer )
{
	FreePointer( _pixels );
	_pixels = AllocPointer< ChannelType >( _w * Channels );
	if( !_pixels ) fprintf( stderr , "[ERROR] WriteImageStream::_Init: Failed to allocate pixels: %d * %d\n" , _w , Channels ) , exit(0);
	info = _InitWrite( fileName , _w , _h , _q , ioServer );
}

template< class ChannelType , int Channels >
WriteImageStream< ChannelType , Channels >::~WriteImageStream(void)
{
	FreePointer( _pixels );
	_FinalizeWrite( info );
}
template< class ChannelType , int Channels > int  WriteImageStream< ChannelType , Channels >::rows( void ) const { return _h; }
template< class ChannelType , int Channels > int  WriteImageStream< ChannelType , Channels >::rowSize( void ) const { return _w * Channels * sizeof(ChannelType); }
template< class ChannelType , int Channels > Pointer( byte ) WriteImageStream< ChannelType , Channels >::operator[] ( int idx ){ return ( Pointer( byte ) )_pixels; }
template< class ChannelType , int Channels >
void WriteImageStream< ChannelType , Channels >::advance( void )
{
	if( current<_h )
	{
		if( _gammaEncode ) for( int x=0 ; x<_w*Channels ; x++ ) _pixels[ x ] = ConvertChannel< double , ChannelType >( pow( ConvertChannel< ChannelType , double >( _pixels[x] ) , 1.0/GAMMA ) );
		if( _clamp ) for( int x=0 ; x<_w*Channels ; x++ ) _pixels[x] = ConvertChannel< double , ChannelType >( std::max< double >( 0. , std::min< double >( 1. , ConvertChannel< ChannelType , double >( _pixels[x] ) ) ) );
		_WriteRow( _pixels , info , current );
	}
	current++;
}
/////////////////////
// ReadImageStream //
/////////////////////
template< class ChannelType , int Channels >
ReadImageStream< ChannelType , Channels >::ReadImageStream( ImageReadInfo< ChannelType > irInfo , char* fileName , int& width , int& height , bool gammaDecode , MultiStreamIOServer* ioServer )
{
	_InitRead		= irInfo.InitRead;
	_ReadRow		= irInfo.ReadRow;
	_FinalizeRead	= irInfo.FinalizeRead;

	_gammaDecode	= gammaDecode;
	current=-1;

	info = _InitRead( fileName , width , height , ioServer );
	_w=width;
	_h=height;

	_pixels = AllocPointer< ChannelType >( _w * Channels );
	if( !_pixels ) fprintf( stderr , "[ERROR] ReadImageStream::ReadImageStream: Failed to allocate pixels\n" ) , exit(0);
	advance();
}
template< class ChannelType , int Channels >
ReadImageStream< ChannelType , Channels >::ReadImageStream( void )
{
	_InitRead		= NULL;
	_ReadRow		= NULL;
	_FinalizeRead	= NULL;

	_gammaDecode	= false;
	current=-1;

	info = NULL;
	_w = 0;
	_h = 0;

	_pixels = NullPointer< ChannelType >();
}
template< class ChannelType , int Channels >
void ReadImageStream< ChannelType , Channels >::Init( ImageReadInfo< ChannelType > irInfo , char* fileName , int& width , int& height , bool gammaDecode , MultiStreamIOServer* ioServer )
{
	_InitRead		= irInfo.InitRead;
	_ReadRow		= irInfo.ReadRow;
	_FinalizeRead	= irInfo.FinalizeRead;

	_gammaDecode	= gammaDecode;
	current=-1;

	info = _InitRead(fileName,width,height , ioServer );
	_w=width;
	_h=height;

	_pixels = AllocPointer< ChannelType >( _w * Channels );
	if( !_pixels ) fprintf( stderr , "[ERROR] ReadImageStream::Init: Failed to allocate pixels\n" ) , exit(0);

	advance();
}
template< class ChannelType , int Channels >
ReadImageStream< ChannelType , Channels >::~ReadImageStream(void)
{
	FreePointer( _pixels );
	_FinalizeRead( info );
}
template< class ChannelType , int Channels > int  ReadImageStream< ChannelType , Channels >::rows( void ) const { return _h; }
template< class ChannelType , int Channels > int  ReadImageStream< ChannelType , Channels >::rowSize( void ) const { return _w * Channels * sizeof( ChannelType ); }
template< class ChannelType , int Channels > Pointer( byte ) ReadImageStream< ChannelType , Channels >::operator[] ( int idx ){ return ( Pointer( byte ) ) _pixels; }
template< class ChannelType , int Channels >
void ReadImageStream< ChannelType , Channels >::advance( void )
{
	current++;
	if( current<_h )
	{
		_ReadRow( _pixels , info , current );
		if( _gammaDecode ) for( int x=0 ; x<_w*Channels ; x++ ) _pixels[ x ] = ConvertChannel< double , ChannelType >( pow( ConvertChannel< ChannelType , double >( _pixels[x] ) , GAMMA ) );
	}
}
////////////////////////
// StreamingInputTile //
////////////////////////
template< class Real , int Channels >
StreamingInputTile< Real , Channels >::StreamingInputTile( void )
{
	rc=1;
	data=NULL;
	w=h=0;
	sX=sY=0;
	rows=NULL;
}
template< class Real , int Channels >
StreamingInputTile< Real , Channels >::~StreamingInputTile( void )
{
	if(data)	delete data;
	if(rows)
	{
		for(int i=0;i<rc;i++)	if(rows[i])	delete[] rows[i];
		delete[] rows;
	}
	data=NULL;
	rows=NULL;
	rc=0;
};

template< class Real , int Channels >
void StreamingInputTile< Real , Channels >::Init( const char* tName , int rowCount , int startX , int startY , bool gammaDecode )
{
	memset( average , 0 , sizeof(double)*Channels );
	strcpy(tileName,tName);
	rc=rowCount;
	sX=startX;
	sY=startY;
	if(rowCount<1)	exit(0);
	rows = new Real*[rc];
	data = GetReadStream< Real >( tileName , w , h , gammaDecode , false , NULL );
	for(int i=-rc;i<0;i++)	init(i);
}

template< class Real , int Channels >
void StreamingInputTile< Real , Channels >::init( int idx )
{
	if( idx+rc-1==sY )
	{
		if(!data)	exit(0);
		for(int i=0;i<rc;i++)	rows[i] = new Real[w*Channels];
		for(int i=0;i<rc-1;i++)
		{
			memcpy( rows[i] , (*data)[i] , w*sizeof(Real)*Channels );
			data->advance();
		}
	}
	else if(idx==sY+h)
	{
		delete data;
		for(int i=0;i<rc;i++) delete[] rows[i];
		delete[] rows;
		data=NULL;
		rows=NULL;
	}
	if(idx>=sY && data)
	{
		memcpy( rows[(idx+rc-1-sY)%rc] , (*data)[idx+rc-1-sY] , w*Channels*sizeof(Real) );
		data->advance();
	}
}
template< class Real , int Channels >
Color< Real , Channels > StreamingInputTile< Real , Channels >::operator() ( int i , int j )
{
	Color< Real , Channels > clr;
	Real* data = &rows[(j-sY)%rc][Channels*(i-sX)];
	for( int c=0 ; c<Channels ; c++ ) clr[c] = data[c];
	return clr;
}
template< class Real , int Channels >
Color< Real , Channels > StreamingInputTile< Real , Channels >::operator() ( int i , int j , Real& a )
{
	Color< Real , Channels > clr;
	Real* data=&rows[(j-sY)%rc][Channels*(i-sX)];
	for( int c=0 ; c<Channels ; c++ ) clr[c] = data[c];
	return clr;
}
template< class Real , int Channels > int StreamingInputTile< Real , Channels >::width ( void ) const { return  w; }
template< class Real , int Channels > int StreamingInputTile< Real , Channels >::height( void ) const { return  h; }
template< class Real , int Channels > int StreamingInputTile< Real , Channels >::startX( void ) const { return sX; }
template< class Real , int Channels > int StreamingInputTile< Real , Channels >::startY( void ) const { return sY; }


///////////////////
// RGridOfImages //
///////////////////
template< class Real , int Channels >
void RGridOfImages< Real , Channels >::GetImageSize( char* fileName , int& width , int& height )
{
	int c , r;
	FILE* fp = fopen(fileName,"r");
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
	char imageName[1024];
	if( fscanf( fp , "Columns: %d " , &c ) != 1) fprintf( stderr , "Failed to read in columns in RGridOfImages\n" ) , exit( 0 );
	if( fscanf( fp , "Rows: %d " , &r )     !=1) fprintf( stderr , "Failed to read in rows in RGridOfImages\n" )    , exit( 0 );
	width = height = 0;
	int* widths  = new int[c];
	int* heights = new int[r];
	for ( int j = 0 ; j < r ; j++ )
		for ( int i = 0 ; i < c ; i++ )
		{
			if( fscanf( fp , " %s " , imageName )!=1 )	fprintf(stderr,"Failed to read in image [%d][%d] in RGridOfImages\n",i,j) , exit(0);
			int ww , hh;
			GetReadSize( imageName , ww , hh );
			if( !j ) widths[i]  = ww	,	width += ww;
			if( !i ) heights[j] = hh	,	height += hh;
			if( widths[i]!=ww )		fprintf( stderr , "Inconsistent widths in column [%d]: %d != %d\n" , i ,  widths[i] , ww );
			if( heights[j]!=hh )	fprintf( stderr , "Inconsistent heights in row [%d]: %d != %d\n"   , j , heights[j] , hh );
		}
	fclose(fp);
	delete[] widths;
	delete[] heights;

}
template< class Real , int Channels >
RGridOfImages< Real , Channels >::RGridOfImages( char* fileName , int& width , int& height , bool gammaDecode , MultiStreamIOServer* ioServer , int threads )
{
	FILE* fp = fopen(fileName,"r");
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
	char imageName[1024];
	if( fscanf( fp , "Columns: %d " , &_c )!=1 ) fprintf( stderr , "Failed to read in columns in RGridOfImages\n" ) , exit(0);
	if( fscanf( fp , "Rows: %d "    , &_r )!=1 ) fprintf( stderr , "Failed to read in rows in RGridOfImages\n" )    , exit(0);

	int w,h;
	_gammaDecode = gammaDecode;
	_threads = threads;
	_widths  = new int[_c];
	_heights = new int[_r];
	_grid = new StreamingGrid**[_c];
	_fileNames = new char**[_c];
	_server = ioServer;
	for ( int i = 0 ; i < _c ; i++ )
	{
		_grid[i] = new StreamingGrid*[_r];
		_fileNames[i] = new char*[_r];
	}
	_w = _h = 0;
	for ( int j = 0 ; j < _r ; j++ )
		for ( int i = 0 ; i < _c ; i++ )
		{
			if(fscanf(fp," %s ",imageName)!=1)	fprintf(stderr,"Failed to read in image [%d][%d] in RGridOfImages\n",i,j) , exit(0);
			_fileNames[i][j]=new char[strlen(imageName)+1];
			strcpy( _fileNames[i][j] , imageName );
			GetReadSize( _fileNames[i][j] , w , h );
			if( !j ) _widths[i]  = w	,	_w += w;
			if( !i ) _heights[j] = h	,	_h += h;
			if( _widths[i]!=w )  fprintf( stderr , "Inconsistent widths in column [%d]: %d != %d\n" , i , _widths[i]  , w ); 
			if( _heights[j]!=h ) fprintf( stderr , "Inconsistent heights in row   [%d]: %d != %d\n" , j , _heights[j] , h );
		}
	fclose(fp);
	_buffer = AllocPointer< Real >( _w * Channels );
	_current = -1;
	advance();
	width = _w;
	height = _h;
}
template< class Real , int Channels >
RGridOfImages< Real , Channels >::~RGridOfImages( void )
{
	delete[] _widths;
	delete[] _heights;
	_widths = _heights = NULL;
	for (int i = 0 ; i < _c ; i++ )
	{
		for ( int j = 0 ; j <_r ; j++ )
		{
			if(_grid[i][j])			delete _grid[i][j];
			if(_fileNames[i][j])	delete _fileNames[i][j];
		}
		delete[] _grid[i];
		delete[] _fileNames[i];
	}
	delete[] _grid;
	delete[] _fileNames;
	FreePointer( _buffer );
	_grid = NULL;
	_fileNames = NULL;
}
template< class Real , int Channels > int  RGridOfImages< Real , Channels >::rows( void ) const { return _h; }
template< class Real , int Channels > int  RGridOfImages< Real , Channels >::rowSize( void ) const { return _w*Channels*sizeof(Real); }
template< class Real , int Channels >
Pointer( byte ) RGridOfImages< Real , Channels >::operator [] (int idx)
{
	return ( Pointer( byte ) )( _buffer );
}
template< class Real , int Channels >
void RGridOfImages< Real , Channels >::advance( void )
{
	_current++;
	if( _current<_h )
	{
		int row;
		int x = 0 , y = _current;
		for ( row = 0 ; y >= _heights[row] ; y -= _heights[row++] ) ;

#pragma omp parallel for num_threads( _threads )
		for( int col=0 ; col<_c ; col++ )
		{
			int x = 0;
			for( int _col=0 ; _col<col ; _col++ ) x += _widths[_col];
			if( !y )
			{
				int w , h;
				_grid[col][row] = GetReadStream< Real , Channels >( _fileNames[col][row] , w , h , _gammaDecode , _server , 0 );
				_grid[col][row]->SetServer( _server );
			}
			if( !_grid[col][row] )	fprintf(stderr,"Error: Attempting to read from NULL tile\n") , exit(0);
			Pointer( Real ) subRow = ( Pointer( Real ) )(*_grid[col][row])[y];

			memcpy( _buffer+Channels*x , subRow , sizeof(Real)*_widths[col]*Channels );

			_grid[col][row]->advance();

			if( y>=_grid[col][row]->rows()-1 )
			{
				delete _grid[col][row];
				_grid[col][row]=NULL;
			}
		}
	}
}
///////////////////
// WGridOfImages //
///////////////////
template< class Real , int Channels >
WGridOfImages< Real , Channels >::WGridOfImages( char* fileName , int width , int height , bool gammaEncode , int quality , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , int threads , bool hdr )
{
	char* header = GetFileHeader( fileName );
	_init( fileName , header , DefaultOutputTileExtension , width , height , gammaEncode , quality , tileWidth , tileHeight , ioServer , threads , hdr );
	delete[] header;
}
template< class Real , int Channels >
WGridOfImages< Real , Channels >::WGridOfImages( char* fileName , const char* header , const char* ext , int width , int height , bool gammaEncode , int quality,int tileWidth,int tileHeight , MultiStreamIOServer* ioServer , int threads , bool hdr )
{
	_init(fileName,header,ext,width,height , gammaEncode , quality,tileWidth,tileHeight , ioServer , threads , hdr );
}
template< class Real , int Channels >
void WGridOfImages< Real , Channels >::_init( char* fileName , const char* header , const char* ext , int width , int height , bool gammaEncode , int quality,int tileWidth,int tileHeight , MultiStreamIOServer* ioServer , int threads , bool hdr )
{
	FILE* fp = fopen( fileName , "w" );
	_gammaEncode = gammaEncode;
	_threads = threads;
	_w = width;
	_h = height;
	_c = (width + tileWidth -1) / tileWidth;
	_r = (height + tileHeight -1) / tileHeight;
	_quality = quality;
	_hdr = hdr;
	_server = ioServer;
	char imageName[1024];
	fprintf( fp , "Columns: %d\n" , _c ) , fflush( fp );
	fprintf( fp , "Rows: %d\n" , _r )    , fflush( fp );

	_widths  = new int[_c];
	_heights = new int[_r];
	for ( int i = 0 ; i < _c ; i++ )	_widths[i]  = tileWidth;
	for ( int j = 0 ; j < _r ; j++ )	_heights[j] = tileHeight;
	_widths[_c-1]  = width  - (_c-1)*tileWidth;
	_heights[_r-1] = height - (_r-1)*tileHeight;

	_grid = new StreamingGrid**[_c];
	_fileNames = new char**[_c];
	for ( int i = 0 ; i < _c ; i++ )
	{
		_grid[i] = new StreamingGrid*[_r];
		_fileNames[i] = new char*[_r];
	}
	for ( int j = 0 ; j < _r ; j++ )
		for ( int i = 0 ; i < _c ; i++ )
		{
			sprintf( imageName , "%s.%d.%d.%s" , header , i , j , ext );
			_fileNames[i][j]=new char[strlen(imageName)+1];
			strcpy(_fileNames[i][j],imageName);
			fprintf(fp,"%s\n",imageName) , fflush( fp );
			_grid[i][j] = NULL;
		}
	fclose(fp);
	_buffer = AllocPointer< Real >( _w * Channels );
	_current = 0;
}
template< class Real , int Channels >
WGridOfImages< Real , Channels >::~WGridOfImages( void )
{
	delete[] _widths;
	delete[] _heights;
	_widths = _heights = NULL;
	for (int i = 0 ; i < _c ; i++ )
	{
		for ( int j = 0 ; j <_r ; j++ )
		{
			if(_grid[i][j])			delete _grid[i][j];
			if(_fileNames[i][j])	delete _fileNames[i][j];
		}
		delete[] _grid[i];
	}
	delete[] _grid;
	FreePointer( _buffer );
	_grid = NULL;
	_fileNames = NULL;
}
template< class Real , int Channels > int  WGridOfImages< Real , Channels >::rows( void ) const { return _h; }
template< class Real , int Channels > int  WGridOfImages< Real , Channels >::rowSize( void ) const { return _w*Channels*sizeof(Real); }
template< class Real , int Channels > Pointer( byte ) WGridOfImages< Real , Channels >::operator [] (int idx){ return ( Pointer( byte ) )( _buffer ); }

template< class Real , int Channels >
void WGridOfImages< Real , Channels >::advance( void )
{
	if( _current<_h )
	{
		int row;
		int x = 0 , y = _current;
		for( row=0 ; y>=_heights[row] ; y-=_heights[row++] ) ;

#pragma omp parallel for num_threads( _threads )
		for( int col=0 ; col<_c ; col++ )
		{
			int x=0;
			for( int _col=0 ; _col<col ; _col++ ) x += _widths[_col];

			if( y==0 ) _grid[col][row] = GetWriteStream< Real , Channels >( _fileNames[col][row] , _widths[col] , _heights[row] , _gammaEncode , _quality , _server , _hdr , 0 );
			if( !_grid[col][row] ) fprintf( stderr , "Error: Attempting to write to NULL tile\n" ) , exit(0);
			Pointer( Real ) subRow = ( Pointer( Real ) )(*_grid[col][row])[y];
			memcpy( subRow , _buffer+Channels*x , sizeof(Real)*_widths[col]*Channels );

			_grid[col][row]->advance();

			if ( y>=_grid[col][row]->rows()-1 )
			{
				delete _grid[col][row];
				_grid[col][row] = NULL;
			}
		}
	}
	_current++;
}



///////////////////
// WImageClipper //
///////////////////
template< class Real , int Channels >
WImageClipper< Real , Channels >::WImageClipper( char* fileName , int width , int height , bool gammaEncode , int outOffX , int outOffY , int outWidth , int outHeight , int quality , MultiStreamIOServer* ioServer , bool hdr , int threads )
{
	deleteGrid = true;
	_hdr = hdr;
	_w = width;
	_h = height;
	_ox = outOffX;
	_oy = outOffY;
	_ow = outWidth;
	_oh = outHeight;

	if(_ox<0 || _oy<0 || _ox+_ow>_w || _oy+_oh>_h)
	{
		fprintf(stderr,"Clip-region not a subset of the domain in WImageClipper: [%d,%d) x [%d,%d) !< [0,%d) x [0,%d)\n",_ox,_ox+_ow,_oy,_oy+_oh,_w,_h);
		exit(0);
	}
	grid = GetWriteStream< Real , Channels >( fileName , _ow , _oh , gammaEncode , quality , ioServer , hdr , threads );
	buffer = AllocPointer< Real >( _w * Channels );
	current = 0;
}
template< class Real , int Channels >
WImageClipper< Real , Channels >::WImageClipper(StreamingGrid* outGrid,int width,int height , int outOffX,int outOffY,int outWidth,int outHeight , bool hdr )
{
	deleteGrid = false;
	_hdr = hdr;
	_w = width;
	_h = height;
	_ox = outOffX;
	_oy = outOffY;
	_ow = outWidth;
	_oh = outHeight;

	if(_ox<0 || _oy<0 || _ox+_ow>_w || _oy+_oh>_h)
	{
		fprintf(stderr,"Clip-region not a subset of the domain in WImageClipper: [%d,%d) x [%d,%d) !< [0,%d) x [0,%d)\n",_ox,_ox+_ow,_oy,_oy+_oh,_w,_h);
		exit(0);
	}

	grid = outGrid;
	buffer = AllocPointer< Real >( _w * Channels );
	current = 0;
}
template< class Real , int Channels >
WImageClipper< Real , Channels >::~WImageClipper(void)
{
	if( deleteGrid ) delete grid;
	FreePointer( buffer );
}
template< class Real , int Channels > int  WImageClipper< Real , Channels >::rows( void ) const { return _h;}
template< class Real , int Channels > int  WImageClipper< Real , Channels >::rowSize( void ) const { return sizeof(Real)*_w*Channels;}

template< class Real , int Channels >
Pointer( byte ) WImageClipper< Real , Channels >::operator [] (int idx)
{
	return ( Pointer( byte ) )( buffer );
}
template< class Real , int Channels >
void WImageClipper< Real , Channels >::advance( void )
{

	if(current>=_oy && current<_oy+_oh)
	{
		Pointer( Real ) row = ( Pointer( Real ) )(*grid)[current-_oy];

		memcpy( row , buffer+Channels*_ox , sizeof(Real)*_ow*Channels );
		grid->advance();
	}
	current++;
}
///////////////////////
// RChannelExtractor //
///////////////////////
template< class Real , int Channels >
RChannelExtractor< Real , Channels >::RChannelExtractor( char* fileName , int& width , int& height , bool gammaDecode , int channel )
{
	_c = channel;
	grid = GetReadStream< Real , Channels >( fileName , _w , _h , gammaDecode );
	width = _w;
	height = _h;

	buffer = new Real[_w];
	current = -1;
	advance();
}
template< class Real , int Channels >
RChannelExtractor< Real , Channels >::~RChannelExtractor( void ){ delete[] buffer; }
template< class Real , int Channels > int  RChannelExtractor< Real , Channels >::rows( void ) const { return _h; }
template< class Real , int Channels > int  RChannelExtractor< Real , Channels >::rowSize( void ) const { return _w*sizeof( Real ); }
template< class Real , int Channels >
Pointer( byte ) RChannelExtractor< Real , Channels >::operator [] (int idx)
{
	return ( Pointer( byte ) )( buffer , rowSize() );
}

template< class Real , int Channels >
void RChannelExtractor< Real , Channels >::advance(void)
{
	current++;
	if(current<_h)
	{
		Real* row = (Real*)(*grid)[current];
		for(int x=0 ; x<_w ; x++)	buffer[x]=row[Channels*x+_c];
		grid->advance();
	}
}

///////////////////
// RImageSampler //
///////////////////
template< class Real , int Channels >
RImageSampler< Real , Channels >::RImageSampler( char* fileName , int width , int  height , bool nearest , bool gammaDecode , MultiStreamIOServer* ioServer , int threads )
{	
	_grid = GetReadStream< Real , Channels >( fileName , _inW , _inH , gammaDecode , ioServer , threads );
	_outW = width;
	_outH = height;
	_nearest = nearest;

	_outBuffer = AllocPointer< Real >( Channels*_outW );
	if( !nearest )
	{
		_inBuffers[0] = AllocPointer< Real >( Channels*_inW );
		_inBuffers[1] = AllocPointer< Real >( Channels*_inW );
	}
	else _inBuffers[0] = _inBuffers[1] = NullPointer< Real >( );
	_inCurrent = 0 , _outCurrent = -1;
	advance();
}
template< class Real , int Channels >
RImageSampler< Real , Channels >::~RImageSampler( void )
{
	FreePointer( _outBuffer );
	FreePointer( _inBuffers[0] );
	FreePointer( _inBuffers[1] );
	delete _grid;
}
template< class Real , int Channels > int  RImageSampler< Real , Channels >::rows( void ) const { return _outH; }
template< class Real , int Channels > int  RImageSampler< Real , Channels >::rowSize(void) const { return _outW * sizeof( Real ) * Channels; }
template< class Real , int Channels >
Pointer( byte ) RImageSampler< Real , Channels >::operator [] ( int idx )
{
	return ( Pointer( byte ) ) _outBuffer;
}
template< class Real , int Channels >
void RImageSampler< Real , Channels >::advance( void )
{
	_outCurrent++;
	double y = (double)( (long long)( _outCurrent ) * (long long)( _inH-1 ) ) / ( _outH-1 );
	int yy;
	if( _nearest )	// Align the input grid to the row we will be reading from
	{
		yy = int( y+0.5 );
		while( _inCurrent<yy && _inCurrent<_inH ) _grid->advance( ) , _inCurrent++;
	}
	else	// Align the input grid to the row we will be reading from and read in the row into the appropriate buffer
	{
		// This needs to be fixed!!!
		yy = int( y );
		while( _inCurrent<yy && _inCurrent<_inH ) _grid->advance( ) , _inCurrent++;
		if( _inCurrent==yy && _inCurrent<_inH )
		{
			Pointer( Real ) row = ( Pointer( Real ) )(*_grid)[_inCurrent];
			memcpy( _inBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * Channels );
			_grid->advance( );
			_inCurrent++;
		}
		if( _inCurrent==yy+1 && _inCurrent<_inH )
		{
			Pointer( Real ) row = ( Pointer( Real ) )(*_grid)[_inCurrent];
			memcpy( _inBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * Channels );
			_grid->advance( );
			_inCurrent++;
		}
	}
	if( _outCurrent<_outH )
	{
		if( _nearest )
		{
			Pointer( Real ) row = ( Pointer( Real ) )(*_grid)[_inCurrent];
			for( int i=0 ; i<_outW ; i++ )
			{
//				float x = float( i * ( _inW-1 ) ) / ( _outW-1 );
				double x = float( (long long)( i ) * (long long)( _inW-1 ) ) / ( _outW-1 );
				int xx = int( x+0.5 );
				for( int c=0 ; c<Channels ; c++ ) _outBuffer[i*Channels+c] = row[Channels*xx+c];
			}
		}
		else
		{
			int y1 = yy;
			int y2 = y1+1;
//			float dy = y-yy;
			double dy = y-yy;
			for( int i=0 ; i<_outW ; i++ )
			{
//				float x = float( i * ( _inW-1 ) ) / ( _outW-1 );
				double x = (double)( (long long)( i ) * (long long)( _inW-1 ) ) / ( _outW-1 );
				int x1 = int( x );
				int x2 = x1 + 1;
//				float dx = x-x1;
				double dx = x-x1;
				if( x2>=_inW ) x2 = x1;
				for( int c=0 ; c<Channels ; c++ )
					_outBuffer[i*Channels+c] =
					_inBuffers[y1&1][x1*Channels+c] * Real(1.0-dx) * Real(1.0-dy) +
					_inBuffers[y1&1][x2*Channels+c] * Real(    dx) * Real(1.0-dy) +
					_inBuffers[y2&1][x1*Channels+c] * Real(1.0-dx) * Real(    dy) +
					_inBuffers[y2&1][x2*Channels+c] * Real(    dx) * Real(    dy) ;
			}
		}
	}
}
/////////////////////////
// RMaskedImageSampler //
/////////////////////////
template< class PixelType , class LabelType , int PChannels , int LChannels >
RMaskedImageSampler< PixelType , LabelType , PChannels , LChannels >::RMaskedImageSampler( char* pixelName , char* labelName , int width , int  height , bool gammaDecode , MultiStreamIOServer* ioServer , int threads )
{
	_pixelGrid = GetReadStream< PixelType , PChannels >( pixelName , _inW , _inH , gammaDecode , ioServer , threads );
	_outW = width;
	_outH = height;
	{
		int w , h;
		_labelGrid = GetReadStream< LabelType , LChannels >( labelName , w , h , false , ioServer , threads );
		if( w!=_inW || h!=_inH ) fprintf( stderr , "[ERROR] Pixel and label dimensions don't match: %d x %d != %d x %d\n" , _inW , _inH , w , h ) , exit( 0 );
	}

	_outBuffer = AllocPointer< PixelType >( PChannels*_outW );
	_pixelBuffers[0] = AllocPointer< PixelType >( PChannels*_inW );
	_pixelBuffers[1] = AllocPointer< PixelType >( PChannels*_inW );
	_labelBuffers[0] = AllocPointer< LabelType >( LChannels*_inW );
	_labelBuffers[1] = AllocPointer< LabelType >( LChannels*_inW );
	_inCurrent = 0 , _outCurrent = -1;
	advance();
}
template< class PixelType , class LabelType , int PChannels , int LChannels >
RMaskedImageSampler< PixelType , LabelType , PChannels , LChannels >::~RMaskedImageSampler( void )
{
	FreePointer( _outBuffer );
	FreePointer( _pixelBuffers[0] );
	FreePointer( _pixelBuffers[1] );
	FreePointer( _labelBuffers[0] );
	FreePointer( _labelBuffers[1] );
	delete _pixelGrid;
	delete _labelGrid;
}
template< class PixelType , class LabelType , int PChannels , int LChannels > int  RMaskedImageSampler< PixelType , LabelType , PChannels , LChannels >::rows( void ) const { return _outH; }
template< class PixelType , class LabelType , int PChannels , int LChannels > int  RMaskedImageSampler< PixelType , LabelType , PChannels , LChannels >::rowSize(void) const { return _outW * sizeof( PixelType ) * PChannels; }
template< class PixelType , class LabelType , int PChannels , int LChannels >
Pointer( byte ) RMaskedImageSampler< PixelType , LabelType , PChannels , LChannels >::operator [] ( int idx )
{
	return ( Pointer( byte ) ) _outBuffer;
}
template< class PixelType , class LabelType , int PChannels , int LChannels >
void RMaskedImageSampler< PixelType , LabelType , PChannels , LChannels >::advance( void )
{
	_outCurrent++;
	double y = (double)( (long long)( _outCurrent ) * (long long)( _inH-1 ) ) / ( _outH-1 );
	int yy = (int)floor( y+0.5 ) , y0 = (int)floor( y ) , y1 = (int)floor( y ) + 1;

	{
		// This needs to be fixed!!!
		while( _inCurrent<y0 && _inCurrent<_inH ) _pixelGrid->advance() , _labelGrid->advance() , _inCurrent++;
		if( _inCurrent==y0 && _inCurrent<_inH )
		{
			Pointer( PixelType ) pixelRow = ( Pointer( PixelType ) )(*_pixelGrid)[_inCurrent];
			Pointer( LabelType ) labelRow = ( Pointer( LabelType ) )(*_labelGrid)[_inCurrent];
			memcpy( _pixelBuffers[_inCurrent&1] , pixelRow , sizeof( PixelType ) * _inW * PChannels );
			memcpy( _labelBuffers[_inCurrent&1] , labelRow , sizeof( LabelType ) * _inW * LChannels );
			_pixelGrid->advance() , _labelGrid->advance();
			_inCurrent++;
		}
		if( _inCurrent==y1 && _inCurrent<_inH )
		{
			Pointer( PixelType ) pixelRow = ( Pointer( PixelType ) )(*_pixelGrid)[_inCurrent];
			Pointer( LabelType ) labelRow = ( Pointer( LabelType ) )(*_labelGrid)[_inCurrent];
			memcpy( _pixelBuffers[_inCurrent&1] , pixelRow , sizeof( PixelType ) * _inW * PChannels );
			memcpy( _labelBuffers[_inCurrent&1] , labelRow , sizeof( LabelType ) * _inW * LChannels );
			_pixelGrid->advance() , _labelGrid->advance();
			_inCurrent++;
		}
	}
	if( _outCurrent<_outH )
	{
		{
			double dy = y-y0;
			Color< LabelType , LChannels > label;
			double weights[2][2];
			for( int i=0 ; i<_outW ; i++ )
			{
				double x = (double)( (long long)( i ) * (long long)( _inW-1 ) ) / ( _outW-1 );
				int xx = (int)floor( x+0.5 ) , x0 = (int)floor( x ) , x1 = (int)floor( x ) + 1;
				double dx = x-x0;
				if( x1>=_inW ) x1 = x0;

				for( int c=0 ; c<LChannels ; c++ ) label[c] = _labelBuffers[yy&1][xx*LChannels+c];
				weights[0][0] = label == Color< LabelType , LChannels >( _labelBuffers[y0&1]+x0*LChannels ) ? ( (1.-dx) * (1.-dy) ) : 0.;
				weights[0][1] = label == Color< LabelType , LChannels >( _labelBuffers[y1&1]+x0*LChannels ) ? ( (1.-dx) * (   dy) ) : 0.;
				weights[1][0] = label == Color< LabelType , LChannels >( _labelBuffers[y0&1]+x1*LChannels ) ? ( (   dx) * (1.-dy) ) : 0.;
				weights[1][1] = label == Color< LabelType , LChannels >( _labelBuffers[y1&1]+x1*LChannels ) ? ( (   dx) * (   dy) ) : 0.;

				double sum = weights[0][0] + weights[0][1] + weights[1][0] + weights[1][1];
				weights[0][0] /= sum , weights[0][1] /= sum , weights[1][0] /= sum , weights[1][1] /= sum;
				for( int c=0 ; c<PChannels ; c++ )
					_outBuffer[i*PChannels+c] =
					_pixelBuffers[y0&1][x0*PChannels+c] * (PixelType) weights[0][0] +
					_pixelBuffers[y0&1][x1*PChannels+c] * (PixelType) weights[1][0] +
					_pixelBuffers[y1&1][x0*PChannels+c] * (PixelType) weights[0][1] +
					_pixelBuffers[y1&1][x1*PChannels+c] * (PixelType) weights[1][1] ;
			}
		}
	}
}
///////////////////
// WImageSampler //
///////////////////
template< class Real , int Channels >
WImageSampler< Real , Channels >::WImageSampler( char* fileName , int inWidth , int  inHeight , int outWidth , int outHeight , bool nearest , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , int threads )
{
	_grid = GetWriteStream< Real , Channels >( fileName , outWidth , outHeight , gammaEncode , quality , ioServer , threads );
	_inW = inWidth;
	_inH = inHeight;
	_outW = outWidth;
	_outH = outHeight;
	_nearest = nearest;

	_inBuffers[0] = AllocPointer< Real >( Channels * _inW );
	_inBuffers[1] = AllocPointer< Real >( Channels * _inW );

	_inCurrent = 0 , _outCurrent = 0;
}
template< class Real , int Channels >
WImageSampler< Real , Channels >::~WImageSampler( void )
{
	FreePointer( _inBuffers[0] );
	FreePointer( _inBuffers[1] );
	delete _grid;
}
template< class Real , int Channels > int  WImageSampler< Real , Channels >::rows( void ) const { return _inH; }
template< class Real , int Channels > int  WImageSampler< Real , Channels >::rowSize( void ) const { return _inW * sizeof( Real ) * Channels; }
template< class Real , int Channels >
Pointer( byte ) WImageSampler< Real , Channels >::operator [] ( int idx )
{
	return ( Pointer( byte ) ) _inBuffers[idx&1];
}
template< class Real , int Channels >
void WImageSampler< Real , Channels >::advance( void )
{
	long long start = (_inCurrent-1) * (_outH-1 ) , end = _inCurrent * (_outH-1);

	while( (long long)( _outCurrent ) * (long long)( _inH-1 )>=start && (long long)( _outCurrent ) * (long long)( _inH-1 )<=end )
	{
		double y = (double)( (long long)( _outCurrent ) * (long long)( _inH-1 ) ) / ( _outH-1 );
		if( _nearest )
		{
			int yy = int( y+0.5 );
			Pointer( Real )  inRow = _inBuffers[yy%2];
			Pointer( Real ) outRow = ( Pointer( Real ) )(*_grid)[_outCurrent];
			for( int i=0 ; i<_outW ; i++ )
			{
				double x = (double)( (long long)( i ) * (long long)( _inW-1 ) ) / ( _outW-1 );
				int xx = int( x+0.5 );
				for( int c=0 ; c<Channels ; c++ ) outRow[i*Channels+c] = inRow[Channels*xx+c];
			}
		}
		else
		{
			Pointer( Real ) inRow1 = _inBuffers[(_inCurrent-1)&1];
			Pointer( Real ) inRow2 = _inBuffers[(_inCurrent  )&1];
			Pointer( Real ) outRow = ( Pointer( Real ) )(*_grid)[_outCurrent];
			double dy = y - (_inCurrent-1);
			for( int i=0 ; i<_outW ; i++ )
			{
				double x = float( (long long)( i ) * (long long)( _inW-1 ) ) / ( _outW-1 );
				double dx = x - int(x);
				int x1 = int(x);
				int x2 = x1 + 1;
				if( x2>=_inW ) x2 = _inW-1;
				double value[Channels];
				for( int c=0 ; c<Channels ; c++ ) value[c] = inRow1[Channels*x1+c]*(1.-dx)*(1.-dy) + inRow1[Channels*x2+c]*(dx)*(1.-dy) + inRow2[Channels*x1+c]*(1.-dx)*(dy) + inRow2[Channels*x2+c]*(dx)*(dy);
				for( int c=0 ; c<Channels ; c++ ) outRow[i*Channels+c] = value[c];
			}
		}
		_outCurrent++;
		_grid->advance( );
	}
	_inCurrent++;
}

////////////////////////////////
#include <Util/CmdLineParser.h>
void GetReadSize( char* fileName , int& width , int& height )
{
	char* ext = GetFileExtension( fileName );
	if     ( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) )	JPGGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "bmp" ) )									BMPGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "kro" ) )									KROGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "png" ) )									PNGGetImageSize( fileName , width , height );
#ifndef NO_TIFF_SUPPORT
	else if( !strcasecmp( ext , "tiff" ) || !strcasecmp( ext , "tif" ) )	TIFGetImageSize( fileName , width , height );
#endif // !NO_TIFF_SUPPORT
	else if( !strcasecmp( ext , "pgm" ) || !strcasecmp( ext , "ppm" ) )		PGMGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "iGrid" ) )									RGridOfImages< unsigned char , 3 >::GetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "int" ) )									RAWint32GetImageSize ( fileName , width , height );
	else if( !strcasecmp( ext , "int16" ) )									RAWint16GetImageSize ( fileName , width , height );
	else if( !strcasecmp( ext , "half" ) )									RAWhalfGetImageSize  ( fileName , width , height );
	else if( !strcasecmp( ext , "float" ) )									RAWfloatGetImageSize ( fileName , width , height );
	else if( !strcasecmp( ext , "double" ) )								RAWdoubleGetImageSize( fileName , width , height );
	else fprintf( stderr , "Unsupported input file extension: %s (%s)\n" , ext , fileName );
	delete[] ext;
}
template< class Real , int Channels >
StreamingGrid* GetReadStream( char* fileName , int& width , int& height , bool gammaDecode , MultiStreamIOServer* ioServer , int threads  )
{
	StreamingGrid* data=NULL;
	char* ext = GetFileExtension( fileName );
	if     ( !strcasecmp(ext,"jpeg") || !strcasecmp(ext,"jpg") )			data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::JPGReadInfo()       , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "kro" ) )									data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::KROReadInfo()       , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "bmp" ) )									data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::BMPReadInfo()       , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "png" ) )									data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::PNGReadInfo()       , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "pfm" ) )									data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::PFMReadInfo()       , fileName , width , height , gammaDecode , ioServer );
#ifndef NO_TIFF_SUPPORT
	else if( !strcasecmp( ext , "tiff" ) || !strcasecmp( ext , "tif" ) )	data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::TIFReadInfo()       , fileName , width , height , gammaDecode , ioServer );
#endif // !NO_TIFF_SUPPORT
	else if( !strcasecmp( ext , "pgm" ) || !strcasecmp( ext , "ppm" ) )		data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::PGMReadInfo()       , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "int"    ) )								data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWint32ReadInfo()  , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "int16"  ) )								data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWint16ReadInfo()  , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "half"   ) )								data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWhalfReadInfo()   , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "float"  ) )								data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWfloatReadInfo()  , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "double" ) )								data = new ReadImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWdoubleReadInfo() , fileName , width , height , gammaDecode , ioServer );
	else if( !strcasecmp( ext , "iGrid" ) )									data = new RGridOfImages< Real , Channels >( fileName , width , height , gammaDecode , ioServer , threads );
	else fprintf( stderr , "Unsupported input file extension: %s (%s)\n" , ext , fileName );
	delete[] ext;
	return data;
}
template< class Real , int Channels >
StreamingGrid* GetWriteStream( char* fileName , const int& width , const int& height , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , bool hdr , int threads )
{
	StreamingGrid* data = NULL;
	char* ext = GetFileExtension( fileName );

	if( hdr )
	{
		if     ( !strcasecmp( ext , "jpg" ) || !strcasecmp( ext , "jpeg" ) )	data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::JPGHDRWriteInfo()    , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "bmp" ) )									data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::BMPHDRWriteInfo()    , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "kro" ) )									data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::KROHDRWriteInfo()    , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "png" ) )									data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::PNGHDRWriteInfo()    , fileName , width , height , gammaEncode , quality , true , ioServer );
#ifndef NO_TIFF_SUPPORT
		else if( !strcasecmp( ext , "tiff" ) || !strcasecmp( ext , "tif" ) )	data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::TIFHDRWriteInfo()    , fileName , width , height , gammaEncode , quality , true , ioServer );
#endif // !NO_TIFF_SUPPORT
		else if( !strcasecmp( ext , "int"    ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWint32WriteInfo()  , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "int16"  ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWint16WriteInfo()  , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "half"   ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWhalfWriteInfo()   , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "float"  ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWfloatWriteInfo()  , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "double" ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWdoubleWriteInfo() , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "iGrid") )									data = new WGridOfImages< Real , Channels >( fileName , width , height , gammaEncode , quality , DefaultOutputTileWidth , DefaultOutputTileHeight , ioServer , threads , hdr );
		else fprintf( stderr , "Unsupported output file extension: %s\n" , ext );
	}
	else
	{
		if     ( !strcasecmp( ext , "jpg" ) || !strcasecmp( ext , "jpeg" ) )	data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::JPGWriteInfo()       , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "bmp" ) )									data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::BMPWriteInfo()       , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "kro" ) )									data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::KROWriteInfo()       , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "png" ) )									data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::PNGWriteInfo()       , fileName , width , height , gammaEncode , quality , true , ioServer );
#ifndef NO_TIFF_SUPPORT
		else if( !strcasecmp(ext,"tiff") || !strcasecmp( ext , "tif" ) )		data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::TIFWriteInfo()       , fileName , width , height , gammaEncode , quality , true , ioServer );
#endif // !NO_TIFF_SUPPORT
		else if( !strcasecmp( ext , "int"    ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWint32WriteInfo()  , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "int16"  ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWint16WriteInfo()  , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "half"   ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWhalfWriteInfo()   , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "float"  ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWfloatWriteInfo()  , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "double" ) )								data = new WriteImageStream< Real , Channels >( ImageReadWriteInfo< Real , Channels >::RAWdoubleWriteInfo() , fileName , width , height , gammaEncode , quality , true , ioServer );
		else if( !strcasecmp( ext , "iGrid") )									data = new WGridOfImages< Real , Channels >( fileName , width , height , gammaEncode , quality ,  DefaultOutputTileWidth , DefaultOutputTileHeight , ioServer , threads , hdr );
		else fprintf( stderr , "Unsupported output file extension: %s\n" , ext );
	}
	delete[] ext;
	return data;
}
