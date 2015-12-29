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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "XPlatform.h"
#include "CmdLineParser.h"


#ifdef WIN32
int strcasecmp(char* c1,char* c2){return _stricmp(c1,c2);}
#endif

cmdLineReadable::cmdLineReadable(const char* name)
{
	set=false;
	this->name=new char[strlen(name)+1];
	strcpy(this->name,name);
}
cmdLineReadable::~cmdLineReadable(void)
{
	if(name) delete[] name;
	name=NULL;
}
int cmdLineReadable::read(char**,int){
	set=true;
	return 0;
}
void cmdLineReadable::write( FILE* fp ) const
{
	fprintf( fp , "--%s" , name );
}
void cmdLineReadable::print( int offset )
{
	for( int i=0 ; i<offset ; i++ ) printf( " " );
	printf( "--%s\n" , name ) , fflush( stdout );
}
void cmdLineReadable::writeValue( char* str )
{
	str[0] = 0;
}

////////////////
// cmdLineInt //
////////////////
cmdLineInt::cmdLineInt(const char* name) : cmdLineReadable(name) {value=0;}
cmdLineInt::cmdLineInt(const char* name,const int& v) : cmdLineReadable(name) {value=v;}
int cmdLineInt::read(char** argv,int argc){
	if(argc>0)
	{
		value=atoi(argv[0]);
		set=true;
		return 1;
	}
	else{return 0;}
}
void cmdLineInt::write( FILE* fp ) const
{
	fprintf( fp , "--%s %d" , name , value );
}
void cmdLineInt::print( int offset )
{
	for( int i=0 ; i<offset ; i++ ) printf( " " );
	printf( "--%s (int) = %d\n" , name , value ) , fflush( stdout );
}
void cmdLineInt::writeValue(char* str)
{
	sprintf(str,"%d",value);
}

/////////////////
// cmdLineInts //
/////////////////
cmdLineInts::cmdLineInts(const char* name) : cmdLineReadable(name)
{
	values = NULL;
	count = 0;
}
cmdLineInts::~cmdLineInts(void)
{
	if(values)	delete[] values;
	values = NULL;
	count = 0;
}
int cmdLineInts::read(char** argv,int argc)
{
	if(argc>0)
	{
		count=atoi(argv[0]);
		if(count <= 0 || argc <= count)	return 1;
		values = new int[count];
		if(!values)	return 0;
		for ( int i = 0 ; i < count ; i++ )	values[i] = atoi ( argv[i+1] );
		set=true;
		return count + 1;
	}
	else return 0;
}
void cmdLineInts::write( FILE* fp ) const
{
	fprintf( fp , "--%s" , name );
	for( int i=0 ; i<count ; i++ ) fprintf( fp , " %d" , values[i] );
}
void cmdLineInts::print( int offset )
{
	for( int i=0 ; i<offset ; i++ ) printf( " " );
	printf( "--%s (int[%d])\n" , name , count ) , fflush( stdout );
}
void cmdLineInts::writeValue( char* str )
{
	char* temp=str;
	for( int i=0 ; i<count ; i++ )
	{
		sprintf( temp , "%d " , values[i] );
		temp = str+strlen(str);
	}
}

//////////////////
// cmdLineFloat //
//////////////////
cmdLineFloat::cmdLineFloat(const char* name) : cmdLineReadable(name) {value=0;}
cmdLineFloat::cmdLineFloat(const char* name, const float& v) : cmdLineReadable(name) {value=v;}
int cmdLineFloat::read(char** argv,int argc){
	if(argc>0){
		value=(float)atof(argv[0]);
		set=true;
		return 1;
	}
	else{return 0;}
}
void cmdLineFloat::write( FILE* fp ) const
{
	fprintf( fp , "--%s %f" , name , value );
}
void cmdLineFloat::print( int offset )
{
	for( int i=0 ; i<offset ; i++ ) printf( " " );
	printf( "--%s (float) = %f\n" , name , value ) , fflush( stdout );
}
void cmdLineFloat::writeValue(char* str)
{
	sprintf(str,"%f",value);
}

///////////////////
// cmdLineString //
///////////////////
cmdLineString::cmdLineString( const char* name , const char* value ) : cmdLineReadable(name)
{
	if( value ) this->value = new char[ strlen(value)+1 ] , strcpy( this->value , value );
	else        this->value = NULL;
}
cmdLineString::~cmdLineString(void)
{
	if(value)	delete[] value;
	value=NULL;
}
int cmdLineString::read(char** argv,int argc){
	if(argc>0)
	{
		value=new char[strlen(argv[0])+1];
		strcpy(value,argv[0]);
		set=true;
		return 1;
	}
	else{return 0;}
}
void cmdLineString::write( FILE* fp ) const
{
	fprintf( fp , "--%s %s" , name , value );
}
void cmdLineString::print( int offset )
{
	for( int i=0 ; i<offset ; i++ ) printf( " " );
	printf( "--%s (char*) = %s\n" , name , value ) , fflush( stdout );
}
void cmdLineString::writeValue( char* str )
{
	sprintf( str , "%s" , value );
}

////////////////////
// cmdLineStrings //
////////////////////
cmdLineStrings::cmdLineStrings(const char* name) : cmdLineReadable(name)
{
	values = NULL;
	count = 0;
}
cmdLineStrings::~cmdLineStrings(void)
{
	if(values)
	{
		for(int i=0;i<count;i++)	if(values[i])	delete[] values[i];
		delete[] values;
	}
	values = NULL;
	count = 0;
}
int cmdLineStrings::read(char** argv,int argc)
{
	if(argc>0)
	{
		count=atoi(argv[0]);
		if(count <= 0 || argc <= count)	return 1;
		values = new char*[count];
		for ( int i = 0 ; i < count ; i++ )
		{
			values[i] = new char[strlen(argv[i+1])+1];
			strcpy(values[i],argv[i+1]);
		}
		set=true;
		return count + 1;
	}
	else return 0;
}
void cmdLineStrings::write( FILE* fp ) const
{
	fprintf( fp ,"--%s" , name );
	for( int  i=0 ; i<count ; i++ ) fprintf( fp , " %s" , values[i] );
}
void cmdLineStrings::print( int offset )
{
	for( int i=0 ; i<offset ; i++ ) printf( " " );
	printf( "--%s (char*[%d])\n" , name , count ) , fflush( stdout );
}
void cmdLineStrings::writeValue( char* str )
{
	char* temp=str;
	for( int i=0 ; i<count ; i++ )
	{
		sprintf( temp , "%s " , values[i] );
		temp = str+strlen(str);
	}
}

char* GetFileExtension( char* fileName )
{
	char* fileNameCopy;
	char* ext = NULL;
	char* temp;

	IOServer::SystemLock lock;

	fileNameCopy = new char[strlen(fileName)+1];
	assert( fileNameCopy );
	strcpy( fileNameCopy , fileName );
	temp = strtok( fileNameCopy , "." );
	while( temp!=NULL )
	{
		if( ext!=NULL ) delete[] ext;
		ext = new char[strlen(temp)+1];
		assert( ext );
		strcpy( ext , temp );
		temp = strtok( NULL , "." );
	}
	delete[] fileNameCopy;
	return ext;
}
char* GetFileHeader(char* fileName)
{
	char* header=NULL;

	header=new char[strlen(fileName)+1];
	assert(header);
	strcpy(header,fileName);
	int i = int( strlen(header) ) - 1;
	while ( i>=0 && header[i]!='.' )	i--;
	header[i] = 0;
	return header;
}

char* GetLocalFileName(char* fileName)
{
	char* fileNameCopy;
	char* name=NULL;
	char* temp;

	IOServer::SystemLock loc;

	fileNameCopy = new char[strlen(fileName)+1];
	assert( fileNameCopy );
	strcpy( fileNameCopy , fileName );
	char tokens[2];
	sprintf( tokens , "%c" , FileSeparator );
	temp = strtok( fileNameCopy , tokens );
	while( temp!=NULL )
	{
		if( name!=NULL ) delete[] name;
		name = new char[strlen(temp)+1];
		assert( name );
		strcpy( name , temp );
		temp = strtok( NULL , tokens );
	}
	delete[] fileNameCopy;
	return name;
}

bool cmdLineParse( int argc , char **argv , int num , cmdLineReadable** readable , bool dumpError )
{
	bool success = true;
	int i,j;
	while (argc > 0)
	{
		if (argv[0][0] == '-' && argv[0][1]=='-')
		{
			for(i=0;i<num;i++)
			{
				if (!strcmp(&argv[0][2],readable[i]->name))
				{
					argv++, argc--;
					j=readable[i]->read(argv,argc);
					argv+=j,argc-=j;
					break;
				}
			}
			if(i==num){
				if( dumpError )
				{
					fprintf(stderr, "invalid option: %s\n",*argv);
					fprintf(stderr, "possible options are:\n");
					for(i=0;i<num;i++)	fprintf(stderr, "  %s\n",readable[i]->name);
					success = false;
				}
				argv++, argc--;
			}
		}
		else
		{
			if(dumpError)
			{
				fprintf(stderr, "invalid option: %s\n", *argv);
				fprintf(stderr, "  options must start with a \'--\'\n");
			}
			argv++, argc--;
		}
	}
	return success;
}
