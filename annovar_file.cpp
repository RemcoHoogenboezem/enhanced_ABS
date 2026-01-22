//----------------------------------------------------------------
// Name        : annovar_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : Open a single annovar file (Can be compressed) and store variants in upper case characters human readable
//----------------------------------------------------------------
#include <iostream>
#include <cstring>
#include <parasail.h>
#include "annovar_file.h"
//----------------------------------------------------------------
void AnnovarFile::Clear(void)
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
}
//----------------------------------------------------------------
AnnovarFile::AnnovarFile(void) {}
//----------------------------------------------------------------
AnnovarFile::AnnovarFile(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples)
{
    Open(filename,fastaEntries,nSamples);
}
//----------------------------------------------------------------
AnnovarFile::~AnnovarFile(void)
{
    Clear();
}
//----------------------------------------------------------------
void AnnovarFile::Open(const string & filename,const map<string,FastaEntry> & fastaEntries,size_t nSamples)
{
    //----------------------------------------------------------------
    //Open annovar file
    //----------------------------------------------------------------

    FILE * annovarFile=fopen(filename.c_str(),"r");

    if(annovarFile==nullptr) throw runtime_error(string("Error: Could not open annovar file: ")+filename);  //Test open annovar file

    fclose(annovarFile);

    annovarFile=popen((string("bzcat -f ")+filename+string(" | zcat -f")).c_str(),"r");

    //----------------------------------------------------------------
    //Read and check annovar header
    //----------------------------------------------------------------

    #define LINE_LEN 16777216   //16MB per line max
    char * line=(char*)malloc(LINE_LEN*sizeof(char));

    if(fgets(line,LINE_LEN,annovarFile)==nullptr)
    {
        free(line); pclose(annovarFile);
        throw runtime_error(string("Error: Could not read header from annovar file: ")+filename);
    }

    if(strcasestr(line,"Chr\tStart\tEnd\tRef\tAlt")==nullptr)
    {
        free(line); pclose(annovarFile);
        throw runtime_error(string("Error: Invalid annovar file header: ")+filename);
    }

    line[strcspn(line,"\n")]='\0';    //Remove new line
    string header(line);

    //----------------------------------------------------------------
    //Read annovar entries
    //----------------------------------------------------------------

    map<string,map<hts_pos_t,vector<pair<hts_pos_t,VarEntry*> > > > entries;

    while(fgets(line,LINE_LEN,annovarFile)!=nullptr)
    {
        if(line[0]=='\0' || line[0]=='\n') continue;

        line[strcspn(line,"\n")]='\0';
        VarEntry * entry=new VarEntry(nSamples,line);

        char * pLine=line; char * pChr=strsep(&pLine,"\t"); char * pPos=strsep(&pLine,"\t"); strsep(&pLine,"\t"); char * pRef=strsep(&pLine,"\t"); char * pAlt=strsep(&pLine,"\t");

        if(pAlt==nullptr)
        {
            free(line); pclose(annovarFile);
            throw runtime_error(string("Error: Insufficient number of columns in annovar file: ")+filename);
        }

        if(pRef[0]=='-' && pAlt[0]=='-')
        {
            free(line); pclose(annovarFile);
            throw runtime_error(string("Error: Invalid variant type in annovar file: ")+filename);
        }

        string chr(pChr);

        if(fastaEntries.count(chr)==0)
        {
            free(line); pclose(annovarFile);
            throw runtime_error(string("Error: Chromosome names in annovar and fasta file must agree: ")+filename);
        }

        const FastaEntry & fastaEntry=fastaEntries.at(chr);

        hts_pos_t pos1=atoll(pPos)-1;

        if(pos1<0 || (pos1 + hts_pos_t(strlen(pRef)))>fastaEntry.len)
        {
            free(line); pclose(annovarFile);
            throw runtime_error(string("Error: Variant position outside fasta reference sequence: ")+filename);
        }

        if(pAlt[0]=='-') pos1--;

        entry->ref=pRef; entry->alt=pAlt;

        if(pRef[0]=='-'){entry->ref=fastaEntry.seq[pos1]; entry->alt=entry->ref+entry->alt;}
        if(pAlt[0]=='-'){entry->alt=fastaEntry.seq[pos1]; entry->ref=entry->alt+entry->ref;}

        for(auto & c : entry->ref) c=toupper(c);
        for(auto & c : entry->alt) c=toupper(c);

        hts_pos_t refLen=hts_pos_t(entry->ref.length()); hts_pos_t altLen=hts_pos_t(entry->alt.length());

        if(refLen+altLen==2)
        {
            entry->varType=SNV;
            entries[chr][pos1].emplace_back(pos1,entry);
            continue;
        }

        if(refLen>=altLen)
        {
            entry->varType=DEL;
            hts_pos_t pos2=pos1+refLen;
            entries[chr][pos1].emplace_back(pos2,entry);
            entries[chr][pos2].emplace_back(pos1,entry);
            continue;
        }

        hts_pos_t pos2=entry->AssessTandem(pos1,fastaEntry.seq+pos1);
        entries[chr][pos1].emplace_back(pos2,entry);
        entries[chr][pos2].emplace_back(pos1,entry);
    }

    free(line); pclose(annovarFile);

    //----------------------------------------------------------------
    //Keep results in object
    //----------------------------------------------------------------

    Clear();

    this->header=header;
    this->entries=entries;
}
//----------------------------------------------------------------
void AnnovarFile::Write(const vector<string> & bamFilenames)
{
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
