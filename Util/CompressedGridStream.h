#ifndef COMPRESSED_GRID_STREAM_INCLUDED
#define COMPRESSED_GRID_STREAM_INCLUDED
#include <Util/ZLIB/zlib.h>
#include "GridStream.h"
class CompressedBufferedStreamingGrid : public StreamingGrid
{
	bool read;
	void* data;
	int current,win;
	StreamingGrid* sg;
public:
	BufferedStreamingGrid	(StreamingGrid* sg);
	~BufferedStreamingGrid	(void);
	int		rows			(void) const;
	int		rowSize			(void) const;
	void*	operator[]		(int idx);
	void	advance			(void);
	void	reset			(bool read,int minWindowSize);
};

#include "CompressedGridStream.inl"
#endif // COMPRESSED_GRID_STREAM_INCLUDED