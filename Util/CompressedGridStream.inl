
/////////////////////////////////////
// CompressedBufferedStreamingGrid //
/////////////////////////////////////
CompressedBufferedStreamingGrid::CompressedBufferedStreamingGrid(StreamingGrid* sg)
{
	data=NULL;
	current=0;
	this->sg=sg;
}
CompressedBufferedStreamingGrid::~CompressedBufferedStreamingGrid(void)
{
	if(data)	free(data);
	data=NULL;
	sg=NULL;
}
int BufferedStreamingGrid::rows		(void) const	{return sg->rows();}
int BufferedStreamingGrid::rowSize	(void) const	{return sg->rowSize();}
void BufferedStreamingGrid::reset	(bool read,int minWindowSize)
{
	this->read=read;
	win=minWindowSize<rows() ? minWindowSize : rows();

	if(data)	free(data), data=NULL;
	data=malloc(rowSize()*win);
	if(!data)
	{
		fprintf(stderr,"Failed to allocate memory for BufferedStreamingGrid\n");
		exit(0);
	}
	sg->reset(read,1);
	if(read)
	{
		for(int w=0;w<win;w++)
		{
			void* row=(*sg)[w];
			memcpy((void*)(long long(data)+w*rowSize()),row,rowSize());
			sg->advance();
		}
	}
	current=0;
}
void BufferedStreamingGrid::advance(void)
{
	if(read)
	{
		current++;
		if(current+win-1<rows())
		{
			sg->advance();
			memcpy((void*)(long long(data)+((current+win-1)%win)*rowSize()),(*sg)[current+win-1],rowSize());
		}
	}
	else
	{
		if(current-win+1>=0)
		{
			memcpy((*sg)[current-win+1],(void*)(long long(data)+((current-win+1)%win)*rowSize()),rowSize());
			sg->advance();
		}
		current++;
	}
}
void* BufferedStreamingGrid::operator[]	(int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx<0 || idx<current || idx>=rows() || idx>=current+win)
		fprintf(stderr,"BufferedStreamingGrid: Index out of bounds: %d\t[%d, %d]\t%d x %d\n",idx,current,current+win,rowSize(),rows()), exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return (void*)(long long(data)+(idx%win)*rowSize());
}
