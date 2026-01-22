//----------------------------------------------------------------
#ifndef AnnovarFileH
#define	AnnovarFileH
//----------------------------------------------------------------
#include <map>
#include "fasta_file.h"
#include "variant_entry.h"
//----------------------------------------------------------------
class AnnovarFile
{
private:

    void Clear(void);

public:

    string                                                          header;
    map<string,map<hts_pos_t,vector<pair<hts_pos_t,VarEntry*> > > > entries;

    AnnovarFile(void);
    AnnovarFile(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples);

    ~AnnovarFile(void);

    void Open(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples);
    void Write(const vector<string> & bamFilenames);
};
//----------------------------------------------------------------
#endif
