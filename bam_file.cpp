//----------------------------------------------------------------
// Name        : bam_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : Open a single bam file to read regions
//----------------------------------------------------------------
#include "bam_file.h"
//----------------------------------------------------------------
void BamFile::Close(void)
{
    if(this->readIterator!=nullptr){sam_itr_destroy(this->readIterator);readIterator=nullptr;}

    if(initialized)
    {
        initialized=false;
        hts_idx_destroy(index); sam_hdr_destroy(header); sam_close(handle);
    }
}
//----------------------------------------------------------------
BamFile::BamFile(void) : initialized(false),readIterator(nullptr) {}
BamFile::BamFile(const string & filename,bool countDuplicates,bool countSecondary) : initialized(false),readIterator(nullptr) {Open(filename,countDuplicates,countSecondary);}
BamFile::~BamFile(void) {Close();}
//----------------------------------------------------------------
void BamFile::Open(const string & filename,bool countDuplicates,bool countSecondary)
{
    if(initialized && filename==handle->fn) return;

    //----------------------------------------------------------------
    //Open bam file
    //----------------------------------------------------------------

    htsFile * handle=sam_open(filename.c_str(),"r");
    if(handle==nullptr)  throw runtime_error(string("Error: Could not open bam file: ")+filename);

    //----------------------------------------------------------------
    //Read header
    //----------------------------------------------------------------

    bam_hdr_t * header=sam_hdr_read(handle);

    if(handle==nullptr)
    {
        sam_close(handle);
        throw runtime_error(string("Error: Could not read bam header: ")+filename);
    }

    //----------------------------------------------------------------
    //Load index
    //----------------------------------------------------------------

    hts_idx_t * index=sam_index_load(handle,filename.c_str());

    if(index==nullptr)
    {
        sam_hdr_destroy(header); sam_close(handle);
        throw runtime_error(string("Error: Could not read bam index: ")+filename);
    }

    //----------------------------------------------------------------
    //Update object
    //----------------------------------------------------------------

    Close();

    this->handle=handle;
    this->header=header;
    this->index=index;

    excludeFlags=BAM_FQCFAIL;
    if(countDuplicates==false) excludeFlags|=BAM_FDUP;
    if(countSecondary==false) excludeFlags|=BAM_FSECONDARY;

    initialized=true;
}
//----------------------------------------------------------------
const string & BamFile::GetFileName(void)
{
    static string fileName;

    if(initialized==false)
    {
        fileName.clear();
        return fileName;
    }

    fileName=handle->fn;
    return fileName;
}
//----------------------------------------------------------------
void BamFile::SetRegion(const string & chr,hts_pos_t pos,hts_pos_t len)
{
    //----------------------------------------------------------------
    //Check if the bam file is initialized
    //----------------------------------------------------------------

    if(initialized==false) throw runtime_error("Error: No open bam file");

    //----------------------------------------------------------------
    //Create iterator
    //----------------------------------------------------------------

    string region(chr+string(":")+to_string(pos+1)+string("-")+to_string(pos+len));
    hts_itr_t * readIterator=sam_itr_querys(index,header,region.c_str());

    if(readIterator==nullptr) throw runtime_error(string("Error: Could not init read iterator: ") + region);

    if(this->readIterator!=nullptr) sam_itr_destroy(this->readIterator);

    this->readIterator=readIterator;
}
//----------------------------------------------------------------
void BamFile::SetRegions(size_t nRegions,char ** regions)
{
    //----------------------------------------------------------------
    //Check if the bam file is initialized
    //----------------------------------------------------------------

    if(initialized==false) throw runtime_error("Error: No open bam file");

    //----------------------------------------------------------------
    //Create iterator
    //----------------------------------------------------------------

    hts_itr_t * readIterator=sam_itr_regarray(index,header,regions,uint(nRegions));

    if(readIterator==nullptr) throw runtime_error("Error: Could not init read iterator for multiple regions");

    if(this->readIterator!=nullptr) sam_itr_destroy(this->readIterator);

    this->readIterator=readIterator;
}
//----------------------------------------------------------------
