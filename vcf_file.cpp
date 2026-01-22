//----------------------------------------------------------------
// Name        : vcf_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : Open a single vcf file
//----------------------------------------------------------------
#include <stdexcept>
#include <vcf.h>
#include "vcf_file.h"
//----------------------------------------------------------------
VCFFile::VCFFile(void){}
//----------------------------------------------------------------
VCFFile::VCFFile(const string & filename){ Open(filename); }
//----------------------------------------------------------------
void VCFFile::Open(const string & filename)
{
    //----------------------------------------------------------------
    //Open file
    //----------------------------------------------------------------

    htsFile * handle=bcf_open(filename.c_str(),"r");

    if(handle==nullptr) throw runtime_error(string("Error: Could not open vcf file: ")+filename);

    //----------------------------------------------------------------
    //Read header
    //----------------------------------------------------------------

    bcf_hdr_t * header=bcf_hdr_read(handle);

    if(header==nullptr)
    {
        bcf_close(handle);
        throw runtime_error(string("Error: Could not read header in vcf file: ")+filename);
    }

    //----------------------------------------------------------------
    //Parse sequece names
    //----------------------------------------------------------------

    int nSeqs; const char ** seqNames=bcf_hdr_seqnames(header,&nSeqs);

    if(seqNames==nullptr)
    {
        bcf_hdr_destroy(header);
        bcf_close(handle);
        throw runtime_error(string("Error: Could not parse seqeuence names: ")+filename);
    }

    //----------------------------------------------------------------
    //Iterate over vcf entries
    //----------------------------------------------------------------

    entries.clear();

    bcf1_t * vcfEntry=bcf_init();

    while(bcf_read(handle,header,vcfEntry)==0)
    {
        if(bcf_is_snp(vcfEntry)==false) continue;
        auto & entry=entries[seqNames[vcfEntry->rid]][vcfEntry->pos]; entry[2]='\0'; entry[3]='\0'; for(size_t i=0,n=min(vcfEntry->n_allele,4U);i<n;i++) entry[i]=vcfEntry->d.allele[i][0];
    }

    //----------------------------------------------------------------
    //Clean up
    //----------------------------------------------------------------

    bcf_destroy(vcfEntry);
    free(seqNames);
    bcf_hdr_destroy(header);
    bcf_close(handle);
}
//----------------------------------------------------------------
