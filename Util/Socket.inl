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

template<class C>
bool ReceiveOnSocket( Socket& s , Pointer( C ) data , size_t dataSize )
{
	unsigned long long rec=0;
	while( rec!=dataSize )
	{
		int tmp = socket_receive( s , ( ( Pointer( char ) ) data) + rec , dataSize-rec );
		if( tmp<=0 )
		{
			IOServer::StderrLock lock;
			if( !tmp )	fprintfId( stderr , "Connection Closed...\n" );
			else		fprintfId( stderr , "socket_receive from client failed (%d): %s\n" , s , LastSocketError() );
			return false;
		}
		rec+=tmp;
	}
	return true;
}
template<class C>
bool SendOnSocket( Socket& s , ConstPointer( C ) data , size_t dataSize )
{
	if( socket_send( s , ( ConstPointer( char ) )data , dataSize )<0 )
	{
		IOServer::StderrLock lock;
		fprintfId( stderr , "socket_send to client failed (%d): %s\n" , s , LastSocketError());
		return false;
	}
	return true;
}
template<class C>
void ReceiveOnSocket( Socket& s , Pointer( C ) data , size_t dataSize , const char* errorMessage , ... )
{
	unsigned long long rec=0;
	while( rec!=dataSize )
	{
		int tmp = socket_receive( s , ( ( Pointer( char ) ) data) + rec , dataSize-rec );
		if( tmp<=0 )
		{
			IOServer::StderrLock lock;
			if(!tmp)	fprintfId( stderr , "Connection Closed...\n" );
			else		fprintfId( stderr , "socket_receive from client failed (%d): %s\n" , s , LastSocketError() );
			{
				fprintf( stderr , "\t" );
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit(0);
		}
		rec+=tmp;
	}
}
template<class C>
void SendOnSocket( Socket& s , ConstPointer( C ) data , size_t dataSize , const char* errorMessage , ... )
{
	if( socket_send( s , ( ConstPointer( char ) )data , dataSize )<0 )
	{
		IOServer::StderrLock lock;
		fprintfId( stderr , "socket_send to client failed (%d): %s\n" , s , LastSocketError());
		{
			fprintf( stderr , "\t" );
			va_list args;
			va_start( args , errorMessage );
			vfprintf( stderr , errorMessage , args );
			va_end( args );
			fprintf( stderr , "\n" );
		}
		exit(0);
	}
}
