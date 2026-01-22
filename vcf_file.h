//----------------------------------------------------------------
#ifndef VCF_FILE_H
#define VCF_FILE_H
//----------------------------------------------------------------
#include "hts.h"
#include <string>
#include <map>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class VCFFile
{
private:
public:

    map<string,map<hts_pos_t,char[4]> > entries;

    VCFFile(void);
    VCFFile(const string & filename);

    void Open(const string & filename);
};
//----------------------------------------------------------------
#endif // VCF_FILE_H
