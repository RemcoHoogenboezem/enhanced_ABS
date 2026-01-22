//----------------------------------------------------------------
#ifndef VEP_FILE_H
#define VEP_FILE_H
//----------------------------------------------------------------
#include <map>
#include "fasta_file.h"
#include "variant_entry.h"
//----------------------------------------------------------------
class VEPFile
{
private:

    void Clear(void);

public:

    vector<string>                                                  info;
    string                                                          header;
    map<string,map<hts_pos_t,vector<pair<hts_pos_t,VarEntry*> > > > entries;

    VEPFile(void);
    VEPFile(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples);

    ~VEPFile(void);

    void Open(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples);
    void Write(const vector<string> & bamFilenames);
};
//----------------------------------------------------------------
#endif // VEP_FILE_H
