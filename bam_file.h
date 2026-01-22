//----------------------------------------------------------------
#ifndef BAM_FILE_H
#define	BAM_FILE_H
//----------------------------------------------------------------
#include <stdexcept>
#include "sam.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class BamFile
{
private:

    bool        initialized;

    uint16_t    excludeFlags;

    htsFile *   handle;
    bam_hdr_t * header;
    hts_idx_t * index;

    hts_itr_t * readIterator;

public:

    BamFile(void);
    BamFile(const string & fileName,bool countDuplicates=false,bool countSecondary=false);
    ~BamFile(void);

    void Open(const string & fileName,bool countDuplicates=false,bool countSecondary=false);
    void Close(void);

    const string & GetFileName(void);

    void SetRegion(const string & chr,hts_pos_t pos,hts_pos_t len);
    void SetRegions(size_t nRegions,char ** regions);

    inline int Read(bam1_t * read)
    {
        int ret; for(ret=sam_read1(handle,header,read);ret>=0 && (read->core.flag&excludeFlags)!=0;ret=sam_read1(handle,header,read));
        if(ret<-1) throw runtime_error(string("Error: Could not read: ")+handle->fn);
        return ret;
    }

    inline int ReadRegion(bam1_t * read)
    {
        int ret; for(ret=sam_itr_next(handle,readIterator,read);ret>=0 && (read->core.flag&excludeFlags)!=0;ret=sam_itr_next(handle,readIterator,read));
        if(ret<-1) throw runtime_error(string("Error: Could not read region: ")+handle->fn);
        return ret;
    }
};
//----------------------------------------------------------------
#endif // BAM_FILE_H
