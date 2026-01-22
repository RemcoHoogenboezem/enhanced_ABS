//----------------------------------------------------------------
// Name        : vep_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : Open a single VEP file and store the variants in human readable format
//----------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <cstring>
#include "vep_file.h"
//----------------------------------------------------------------
void VEPFile::Clear(void)
{
    for(const auto & chr : entries)
    {
        for(const auto & pos : chr.second)
        {
            for(const auto & entry : pos.second)
            {
                if(pos.first<=entry.first) delete entry.second;
            }
        }
    }

    entries.clear();
}
//----------------------------------------------------------------
VEPFile::VEPFile(void){}
//----------------------------------------------------------------
VEPFile::VEPFile(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples)
{
    Open(filename,fastaEntries,nSamples);
}
//----------------------------------------------------------------
VEPFile::~VEPFile(void)
{
    Clear();
}
//----------------------------------------------------------------
void VEPFile::Open(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples)
{
    //----------------------------------------------------------------
    //Open VEP file
    //----------------------------------------------------------------

    FILE * vepFile=fopen(filename.c_str(),"r");

    if(vepFile==nullptr) throw runtime_error(string("Error: Could not open vep file: ")+filename);

    fclose(vepFile);

    vepFile=popen((string("bzcat -f ")+filename+string(" | zcat -f")).c_str(),"r");

    //----------------------------------------------------------------
    //Read info + header
    //----------------------------------------------------------------

    #define LINE_LEN 16777216
    char * line=(char*)malloc(LINE_LEN);

    vector<string> info;

    while(fgets(line,LINE_LEN,vepFile)!=nullptr && line[0]=='#' && line[1]=='#')
    {
        line[strcspn(line,"\n")]='\0';
        info.push_back(line);
    }

    if(line[0]!='#' || line[1]=='#')
    {
        free(line);pclose(vepFile);
        throw runtime_error(string("Error: No header line found in vep file: ")+filename);
    }

    line[strcspn(line,"\n")]='\0';

    string header(line);
    size_t nFields=count(header.begin(),header.end(),'\t');

    //----------------------------------------------------------------
    //Read entries and remove duplicate variants
    //----------------------------------------------------------------

    map<string,vector<vector<string>>> entries;
    {
        while(fgets(line,LINE_LEN,vepFile)!=nullptr)
        {
            if(line[0]=='\0' || line[0]=='\n' || line[0]=='#') continue;

            line[strcspn(line,"\n")]='\0';
            char * pLine=line; char * pVar=strsep(&pLine,"\t"); char * pChr=strsep(&pVar,"_"); char * pPos=strsep(&pVar,"_"); char * pRef=strsep(&pVar,"_"); char * pAlt=strsep(&pVar,"_");

            if(pVar!=nullptr || pAlt==nullptr)
            {
                free(line); pclose(vepFile);
                throw runtime_error(string("Error: Inconsistent number of columns in VEP entry: ")+filename);
            }

            string chr(pChr);

            if(fastaEntries.count(chr)==0)
            {
                free(line); pclose(vepFile);
                throw runtime_error(string("Error: Chromosome names in vep and fasta file must agree: ")+filename);
            }

            const FastaEntry & fastaEntry=fastaEntries.at(chr);

            hts_pos_t pos=atoll(pPos)-1;

            if(pos<0 || (pos + hts_pos_t(strlen(pRef)))>fastaEntry.len)
            {
                free(line); pclose(vepFile);
                throw runtime_error(string("Error: Variant position outside fasta reference sequence: ")+filename);
            }

            string var(pVar);

            if(entries.count(var)==0) entries[var].assign(nFields,vector<string>());

            size_t field=0;for(char * pField=strsep(&pLine,"\t");pField!=nullptr;pField=strsep(&pLine,"\t"),field++)
            {
                auto & vec=entries[var][field]; auto vecEnd=vec.end();
                if(find(vec.begin(),vecEnd,pField)==vecEnd) vec.push_back(pField);
            }

            if(field<nFields)
            {
                free(line); pclose(vepFile);
                throw runtime_error(string("Error: Inconsistent number of columns in VEP entry: ")+filename);
            }
        }
    }

    pclose(vepFile);

    this->info=info;
    this->header=header;

    //----------------------------------------------------------------
    //Convert entries
    //----------------------------------------------------------------

    Clear();

    for(auto & entry : entries)
    {
        VarEntry * newEntry=new VarEntry(nSamples,entry.first);

        char * pVar=strcpy(line,entry.first.c_str());
        string chr(strsep(&pVar,"_")); hts_pos_t pos1=atoll(strsep(&pVar,"_"))-1; newEntry->ref=strsep(&pVar,"_"); newEntry->alt=strsep(&pVar,"_");

        for(auto & field : entry.second)
        {
            auto it=field.begin();
            auto itEnd=field.end();

            newEntry->line+=string("\t")+*it;
            for(it++;it!=itEnd;it++) newEntry->line+=string(";")+*it;
        }

        for(auto & c : newEntry->ref) c=toupper(c);
        for(auto & c : newEntry->alt) c=toupper(c);

        hts_pos_t refLen=hts_pos_t(newEntry->ref.length()); hts_pos_t altLen=hts_pos_t(newEntry->alt.length());

        if(refLen+altLen==2)
        {
            newEntry->varType=SNV;
            this->entries[chr][pos1].emplace_back(pos1,newEntry);
            continue;
        }

        if(refLen>=altLen)
        {
            newEntry->varType=DEL;
            hts_pos_t pos2=pos1+refLen-1;
            this->entries[chr][pos1].emplace_back(pos2,newEntry);
            this->entries[chr][pos2].emplace_back(pos1,newEntry);
            continue;
        }

        hts_pos_t pos2=newEntry->AssessTandem(pos1,fastaEntries.at(chr).seq+pos1);
        this->entries[chr][pos1].emplace_back(pos2,newEntry);
        this->entries[chr][pos2].emplace_back(pos1,newEntry);
    }

    free(line);
}
//----------------------------------------------------------------
void VEPFile::Write(const vector<string> & bamFilenames)
{
    //----------------------------------------------------------------
    //Write info lines
    //----------------------------------------------------------------

    for(const auto & line : info) cout << line << '\n';

    //----------------------------------------------------------------
    //Write header
    //----------------------------------------------------------------

    cout << header << "\tvar_type";
    size_t nBamFiles=bamFilenames.size(); for(size_t i=0;i<nBamFiles;i++) Statistics::PrintHeader(bamFilenames[i]);
    cout << '\n';

    //----------------------------------------------------------------
    //Write entries
    //----------------------------------------------------------------

    for(const auto & chr : entries)
    {
        for(const auto & pos : chr.second)
        {
            for(const auto & entry : pos.second)
            {
                if(pos.first<=entry.first)
                {
                    cout << entry.second->line << '\t' << varStrings[entry.second->varType];
                    const auto & statistics=entry.second->statistics; for(size_t i=0;i<nBamFiles;i++) statistics[i].PrintEntry();
                    cout << '\n';
                }
            }
        }
    }
}
//----------------------------------------------------------------

