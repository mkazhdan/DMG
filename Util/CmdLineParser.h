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

#ifndef CMD_LINE_PARSER_INCLUDED
#define CMD_LINE_PARSER_INCLUDED
#include <stdarg.h>
#include <string.h>


#ifdef WIN32
#ifndef strcasecmp
#define strcasecmp _stricmp
#endif // strcasecmp
#endif // WIN32

class cmdLineReadable{
public:
	bool set;
	char* name;
	cmdLineReadable(const char* name);
	virtual ~cmdLineReadable(void);
	virtual int read(char** argv,int argc);
	virtual void write( FILE* fp ) const;
	virtual void print( int offset=0 );
	virtual void writeValue( char* str );
};

class cmdLineInt : public cmdLineReadable {
public:
	int value;
	cmdLineInt(const char* name);
	cmdLineInt(const char* name,const int& v);
	int read(char** argv,int argc);
	void write( FILE* fp ) const;
	void print( int offset=0 );
	void writeValue( char* str );
};
class cmdLineInts : public cmdLineReadable
{
public:
	int count;
	int* values;
	cmdLineInts(const char* name);
	~cmdLineInts(void);
	int read( char** argv , int argc );
	void write( FILE* fp ) const;
	void print( int offset=0 );
	void writeValue( char* str );
};

template<int Dim>
class cmdLineIntArray : public cmdLineReadable {
public:
	int values[Dim];
	cmdLineIntArray(const char* name);
	cmdLineIntArray(const char* name,const int v[Dim]);
	int read(char** argv,int argc);
	void write( FILE* fp ) const;
	void print( int offset=0 );
	void writeValue( char* str );
};

class cmdLineFloat : public cmdLineReadable {
public:
	float value;
	cmdLineFloat(const char* name);
	cmdLineFloat(const char* name,const float& f);
	int read(char** argv,int argc);
	void write( FILE* fp ) const;
	void print( int offset=0 );
	void writeValue( char* str );
};
template<int Dim>
class cmdLineFloatArray : public cmdLineReadable {
public:
	double values[Dim];
	cmdLineFloatArray(const char* name);
	cmdLineFloatArray(const char* name,const float f[Dim]);
	int read(char** argv,int argc);
	void write( FILE* fp ) const;
	void print( int offset=0 );
	void writeValue( char* str );
};
class cmdLineString : public cmdLineReadable {
public:
	char* value;
	cmdLineString( const char* name , const char* value=NULL );
	~cmdLineString();
	int read(char** argv,int argc);
	void write( FILE* fp ) const;
	void print( int offset=0 );
	void writeValue( char* str );
};
class cmdLineStrings : public cmdLineReadable {
public:
	int count;
	char** values;
	cmdLineStrings(const char* name);
	~cmdLineStrings();
	int read(char** argv,int argc);
	void write( FILE* fp ) const;
	void print( int offset=0 );
	void writeValue( char* str );
};


// This reads the arguments in argc, matches them against "names" and sets
// the values of "r" appropriately. Parameters start with "--"
bool cmdLineParse( int argc , char **argv , int num , cmdLineReadable** r , bool dumpError=true );

char* GetFileExtension(char* fileName);
char* GetFileHeader(char* fileName);
char* GetLocalFileName(char* fileName);

#include "Util/CmdLineParser.inl"
#endif // CMD_LINE_PARSER_INCLUDED
