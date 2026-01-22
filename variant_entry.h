//----------------------------------------------------------------
#ifndef VARIANT_ENTRY_H
#define VARIANT_ENTRY_H
//----------------------------------------------------------------
#include <iostream>
#include <vector>
#include <stdint.h>
#include "hts.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class Statistics
{
private:
public:

    uint32_t totalDepth;
    uint32_t unknown;
    uint32_t ambiguous;
    uint32_t altDepth;
    uint32_t altBias[2];

    uint32_t hqDepth;
    uint32_t hqAltDepth;
    uint32_t hqAltBias[2];

    uint32_t mismatchesTooHigh;

    hts_pos_t posSum;
    hts_pos_t posSum2;

    float vafITD;
    float coeITD;

    Statistics(void) : totalDepth(0U),unknown(0U),ambiguous(0U),altDepth(0U),altBias{0U,0U},hqDepth(0U),hqAltDepth(0U),hqAltBias{0U,0U},mismatchesTooHigh(0U),posSum(0U),posSum2(0U),vafITD(0.0f),coeITD(0.0f){}

    static inline void PrintHeader(const string & fileName)
    {
        size_t pos=fileName.find_last_of("\\/")+1;
        string sampleName=fileName.substr(pos,fileName.find_last_of('.')-pos);

        cout    << '\t' << sampleName << ":"
                "\ttotal_depth"
                "\tunknown"
                "\tambiguous"
                "\tdepth"
                "\talt_depth"
                "\talt_freq"
                "\talt_bias"
                "\tHQ_depth"
                "\tHQ_alt_depth"
                "\tHQ_alt_freq"
                "\tHQ_alt_bias"
                "\tHQ_ratio"
                "\tmismatches_too_high"
                "\tstart_pos_var"
                "\tVAF_ITD"
                "\tCOE_ITD";
    }

    inline void PrintEntry(void) const
    {
        uint32_t depth=totalDepth-(unknown+ambiguous);
        uint32_t altReadDepth=max(altBias[0]+altBias[1],1U);

        cout    << "\t\t" << totalDepth
                << '\t' << unknown
                << '\t' << ambiguous
                << '\t' << depth
                << '\t' << altDepth
                << '\t' << double(altDepth) / double(max(depth,1U))
                << '\t' << double(altBias[0]) / double(altReadDepth)
                << '\t' << hqDepth
                << '\t' << hqAltDepth
                << '\t' << double(hqAltDepth) / double(max(hqDepth,1U))
                << '\t' << double(hqAltBias[0]) / double(max(hqAltBias[0]+hqAltBias[1],1U))
                << '\t' << double(hqDepth)/double(max(depth,1U))
                << '\t' << double(mismatchesTooHigh)/double(altReadDepth)
             << '\t' << double(posSum2)-double(posSum*posSum)/double(altReadDepth)
                << '\t' << vafITD
                << '\t' << coeITD;
    }
};
//----------------------------------------------------------------
enum VarType {SNV=0,DEL=1,INS=2,ITD=3,PTD=4};
inline char varStrings[][4]={"SNV","DEL","INS","ITD","PTD"};
//----------------------------------------------------------------
class VarEntry
{
private:
public:

    VarType varType;
    string line,ref,alt;
    vector<Statistics> statistics;

    VarEntry(size_t nBamFiles) { statistics.assign(nBamFiles,Statistics()); }
    VarEntry(size_t nBamFiles,char * line) : line(line) { statistics.assign(nBamFiles,Statistics()); }
    VarEntry(size_t nBamFiles,const string & line) : line(line) { statistics.assign(nBamFiles,Statistics()); }

    hts_pos_t AssessTandem(hts_pos_t pos1,const char * ref)
    {
        const char * alt=this->alt.c_str();
        hts_pos_t pos2; for(pos2=pos1;*ref==*alt;pos2++,ref++,alt++){} pos2=max(pos1+1,pos2-1);

        if(pos1+1==pos2) {varType=INS; return pos2;}
        if(*alt=='\0') {varType=ITD; return pos2;}
        varType=PTD; return pos2;
    }

};
//----------------------------------------------------------------
#endif // VARIANT_ENTRY_H
