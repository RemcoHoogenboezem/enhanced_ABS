//----------------------------------------------------------------
// Name        : annotate_bam_statistics.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <set>
#include <omp.h>
#include <getopt.h>
#include "vcf_file.h"
#include "annovar_file.h"
#include "vep_file.h"
#include "bam_file.h"
#include "progress_bar.h"
#include "ssw.h"
#include "annotate_bam_statistics.h"
//----------------------------------------------------------------
class Job
{
private:
public:

    size_t fileIndex;
    const char * chr;
    map<hts_pos_t,vector<pair<hts_pos_t,VarEntry*> > >::iterator pos;

    Job(size_t fileIndex,const char * chr,map<hts_pos_t,vector<pair<hts_pos_t,VarEntry*> > >::iterator pos) : fileIndex(fileIndex),chr(chr),pos(pos) {}
};
//----------------------------------------------------------------
//class Read;
//----------------------------------------------------------------
class Thread
{
private:

    bool countDuplicates;
    bool countSecondary;

    uint8_t minHQBaseScore;
    uint8_t minHQAlignmentScore;

    hts_pos_t minMatchLength;
    hts_pos_t pileupTolerance;

    size_t nBamFiles;
    size_t nMaxMismatches;

    double minAlignmentRate;

    const vector<string> & bamFileNames;

    const map<string,FastaEntry> & fastaEntries;
    const map<string,map<hts_pos_t,char[4]> > & vcfEntries;
    map<string,map<hts_pos_t,vector<pair<hts_pos_t,VarEntry*> > > > & varEntries;

    BamFile bamFile;

    void ProcessSNP(Job & job);
    void ProcessIndel(Job & job);

public:

    Thread(bool countDuplicates,bool countSecondary,uint8_t minHQBaseScore,uint8_t minHQAlignmentScore,hts_pos_t minMatchLength,hts_pos_t pileupTolerance,size_t nBamFiles,double minAlignmentRate,const vector<string> & bamFileNames,const map<string,FastaEntry> & fastaEntries,const map<string,map<hts_pos_t,char[4]> > & vcfEntries,map<string,map<hts_pos_t,vector<pair<hts_pos_t,VarEntry*> > > > & varEntries)
        : countDuplicates(countDuplicates),countSecondary(countSecondary),minHQBaseScore(minHQBaseScore),minHQAlignmentScore(minHQAlignmentScore),minMatchLength(minMatchLength),pileupTolerance(pileupTolerance),nBamFiles(nBamFiles),minAlignmentRate(minAlignmentRate),bamFileNames(bamFileNames),fastaEntries(fastaEntries),vcfEntries(vcfEntries),varEntries(varEntries) {}

    void RunJob(Job & job);
};
//----------------------------------------------------------------
void Thread::ProcessSNP(Job & job)
{
    //----------------------------------------------------------------
    //Get fasta sequence and vcf entries by chromosome
    //----------------------------------------------------------------

    const char * fastaSeq=fastaEntries.at(job.chr).seq;
    const map<hts_pos_t,char[4]> * vcfEntriesByChr=vcfEntries.count(job.chr)!=0 ? &vcfEntries.at(job.chr) : nullptr;

    //----------------------------------------------------------------
    //Set pileup position SNP
    //----------------------------------------------------------------

    hts_pos_t snpPos=job.pos->first;
    bamFile.SetRegion(job.chr,snpPos,1);

    //----------------------------------------------------------------
    //Read class SNP
    //----------------------------------------------------------------

    class Read
    {
    private:
    public:

        bool empty;
        bool mismatchesTooHigh;

        uint8_t alignmentScore;
        char base;
        uint8_t baseScore;

        hts_pos_t pos;

        Read(void) : empty(true){}
        void Init(bool mismatchesTooHigh,uint8_t alignmentScore,char base,uint8_t baseScore,hts_pos_t pos){empty=false;this->mismatchesTooHigh=mismatchesTooHigh;this->alignmentScore=alignmentScore;this->base=base;this->baseScore=baseScore;this->pos=pos;}
    };

    //----------------------------------------------------------------
    //Iterate over reads and gather fragments
    //----------------------------------------------------------------

    map<string,Read[2]> fragments;

    static const char bamSeq2ASCII[]="NACNGNNNTNNNNNNN";

    bam1_t * read=bam_init1(); while(bamFile.ReadRegion(read)>=0)
    {
        hts_pos_t refPos=read->core.pos;
        size_t queryPos=0;

        hts_pos_t matchLength=0;
        size_t nMismatches=0;

        char base; uint8_t baseScore;

        for(uint32_t * cigar=bam_get_cigar(read),*cigarEnd=cigar+size_t(read->core.n_cigar);cigar!=cigarEnd;cigar++)
        {
            hts_pos_t opLen=bam_cigar_oplen(*cigar);

            switch(bam_cigar_op(*cigar))
            {
            //Consumes reference only
            case BAM_CDEL:
            case BAM_CREF_SKIP:

                if(refPos<=snpPos && snpPos<refPos+opLen) goto NEXT_READ;   //Skip if snpPos is inside deletion or refskip
                refPos+=opLen;
                continue;

            //Consumes query only
            case BAM_CINS:
            case BAM_CSOFT_CLIP:

                queryPos+=opLen;
                continue;

            //Consumes reference and query
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:

                matchLength+=opLen;

                uint8_t * bamSeq=bam_get_seq(read);

                for(size_t i=refPos,iEnd=refPos+opLen,j=queryPos;i<iEnd;i++,j++)
                {
                    if(vcfEntriesByChr==nullptr || vcfEntriesByChr->count(i)==0)
                    {
                        nMismatches+=fastaSeq[i]!=bamSeq2ASCII[bam_seqi(bamSeq,j)];
                        continue;
                    }

                    auto vcfBases=vcfEntriesByChr->at(i); char bamBase=bamSeq2ASCII[bam_seqi(bamSeq,j)];
                    nMismatches+=vcfBases[0]!=bamBase && vcfBases[1]!=bamBase && vcfBases[2]!=bamBase && vcfBases[3]!=bamBase;
                }

                if(refPos<=snpPos && snpPos<refPos+opLen)
                {
                    queryPos+=snpPos-refPos;
                    base=bamSeq2ASCII[bam_seqi(bamSeq,queryPos)]; baseScore=bam_get_qual(read)[queryPos];
                }

                refPos+=opLen;
                queryPos+=opLen;
                continue;
            }
        }

        if(matchLength<minMatchLength) continue;

        fragments[bam_get_qname(read)][bam_is_rev(read)].Init(nMismatches>nMaxMismatches,read->core.qual,base,baseScore,read->core.pos);

    NEXT_READ:;
    }

    bam_destroy1(read);

    //----------------------------------------------------------------
    //Compute statistics
    //----------------------------------------------------------------

    auto & entries=job.pos->second;
    size_t fileIndex=job.fileIndex;
    size_t nFragments=fragments.size();
    for(auto & entry : entries) entry.second->statistics[fileIndex].totalDepth=nFragments;

    for(auto & fragment : fragments)
    {
        auto & fwRead=fragment.second[0]; auto & rvRead=fragment.second[1];

        if(fwRead.empty==false)
        {
            bool fwReadIsHQ=fwRead.alignmentScore>=minHQAlignmentScore && fwRead.baseScore>=minHQBaseScore;

            if(rvRead.empty==false)
            {
                bool rvReadIsHQ=rvRead.alignmentScore>=minHQAlignmentScore && rvRead.baseScore>=minHQBaseScore;
                bool fragmentIsHQ=fwReadIsHQ || rvReadIsHQ;

                for(auto & entry : entries)
                {
                    auto varEntry=entry.second;
                    bool fwReadIsAlt=fwRead.base==varEntry->alt[0];
                    bool rvReadIsAlt=rvRead.base==varEntry->alt[0];

                    auto & statistics=varEntry->statistics[fileIndex];
                    statistics.hqDepth+=fragmentIsHQ;
                    statistics.altBias[0]+=fwReadIsAlt;
                    statistics.altBias[1]+=rvReadIsAlt;
                    statistics.hqAltBias[0]+=fwReadIsHQ&&fwReadIsAlt;
                    statistics.hqAltBias[1]+=rvReadIsHQ&&rvReadIsAlt;

                    if(fwRead.baseScore>rvRead.baseScore){ statistics.altDepth+=fwReadIsAlt; statistics.hqAltDepth+=fwReadIsHQ&&fwReadIsAlt; continue;}
                    if(fwRead.baseScore<rvRead.baseScore){ statistics.altDepth+=rvReadIsAlt; statistics.hqAltDepth+=rvReadIsHQ&&rvReadIsAlt; continue;}
                    statistics.altDepth+=fwReadIsAlt||rvReadIsAlt; statistics.hqAltDepth+=(fwReadIsHQ&&fwReadIsAlt)||(rvReadIsHQ&&rvReadIsAlt);

                    statistics.mismatchesTooHigh+=fwReadIsAlt&&fwRead.mismatchesTooHigh;
                    statistics.mismatchesTooHigh+=rvReadIsAlt&&rvRead.mismatchesTooHigh;
                }

                continue;
            }

            for(auto & entry : entries)
            {
                auto varEntry=entry.second;
                bool fwReadIsAlt=fragment.second[0].base==varEntry->alt[0];

                auto & statistics=varEntry->statistics[fileIndex];
                statistics.hqDepth+=fwReadIsHQ;
                statistics.altBias[0]+=fwReadIsAlt;
                statistics.hqAltBias[0]+=fwReadIsHQ&&fwReadIsAlt;
                statistics.altDepth+=fwReadIsAlt;
                statistics.hqAltDepth+=fwReadIsHQ&&fwReadIsAlt;

                statistics.mismatchesTooHigh+=fwReadIsAlt&&fwRead.mismatchesTooHigh;
            }

            continue;
        }

        bool rvReadIsHQ=fragment.second[1].alignmentScore>=minHQAlignmentScore && fragment.second[1].baseScore>=minHQBaseScore;

        for(auto & entry : entries)
        {
            auto varEntry=entry.second;
            bool rvReadIsAlt=fragment.second[1].base==varEntry->alt[0];

            auto & statistics=varEntry->statistics[fileIndex];
            statistics.hqDepth+=rvReadIsHQ;
            statistics.altBias[1]+=rvReadIsAlt;
            statistics.hqAltBias[1]+=rvReadIsHQ&&rvReadIsAlt;
            statistics.altDepth+=rvReadIsAlt;
            statistics.hqAltDepth+=rvReadIsHQ&&rvReadIsAlt;

            statistics.mismatchesTooHigh+=rvReadIsAlt&&rvRead.mismatchesTooHigh;
        }
    }
}
//----------------------------------------------------------------
void Thread::ProcessIndel(Job & job)
{
    //----------------------------------------------------------------
    //Calculate and set pileup regions
    //----------------------------------------------------------------

    hts_pos_t jobPos1=job.pos->first; hts_pos_t jobPos2=0;

    auto fastaEntry=fastaEntries.at(job.chr); hts_pos_t fastaEntryEnd=fastaEntry.len-1;
    {
        map<hts_pos_t,hts_pos_t> intervals;

        intervals[max(jobPos1-pileupTolerance,0L)]=min(jobPos1+pileupTolerance,fastaEntryEnd);
        for(auto  entry : job.pos->second) {hts_pos_t pos2=entry.first; jobPos2=max(jobPos2,pos2); intervals[max(pos2-pileupTolerance,0L)]=min(pos2+pileupTolerance,fastaEntryEnd);}

        vector<pair<hts_pos_t,hts_pos_t> > mergedIntervals; mergedIntervals.reserve(16);

        auto it=intervals.begin(); auto itEnd=intervals.end();
        hts_pos_t intervalBegin=it->first; hts_pos_t intervalEnd=it->second;

        for(it++;it!=itEnd;it++)
        {
            if(it->first<=intervalEnd+1)
            {
                intervalEnd=max(intervalEnd,it->second);
                continue;
            }

            mergedIntervals.emplace_back(intervalBegin,intervalEnd);
            intervalBegin=it->first; intervalEnd=it->second;
        }

        mergedIntervals.emplace_back(intervalBegin,intervalEnd);

        size_t nRegions=mergedIntervals.size(); char ** regions=(char**)malloc(nRegions*sizeof(char*));

        char region[1024]; const char * chr=job.chr;

        for(size_t i=0;i<nRegions;i++)
        {
            regions[i]=strcpy((char*)malloc(sprintf(region,"%s:%li-%li",chr,mergedIntervals[i].first,mergedIntervals[i].second)+1),region);
        }

        bamFile.SetRegions(nRegions,regions);

        for(size_t i=0;i<nRegions;i++) {free(regions[i]);} free(regions);
    }

    //----------------------------------------------------------------
    //Gather fragments and determine reference intervals
    //----------------------------------------------------------------

    class Reference
    {
    private:
    public:

        bool inJobSet;
        SSW ssw;
        VarEntry * varEntry;

        Reference(bool inJobSet,const string & refSeq,VarEntry * varEntry) : inJobSet(inJobSet),ssw(refSeq),varEntry(varEntry) {}
    };

    class Read
    {
    private:
    public:

        bool isHQ;
        bool isUnknown;
        bool isAmbiguous;

        int maxPosScore;

        hts_pos_t begin,end;
        char * query;

        Read(void) : query(nullptr) {}
        ~Read(void) {if(query!=nullptr) free(query);}

        void Init(const bam1_t * bamRead,uint8_t minHQAlignmentScore,uint8_t minHQBaseScore,hts_pos_t & readBegin,hts_pos_t & readEnd)
        {
            isHQ=bamRead->core.qual>=minHQAlignmentScore;

            readBegin=readEnd=bamRead->core.pos;

            for(uint32_t *cigarBegin=bam_get_cigar(bamRead),*cigarEnd=cigarBegin+bamRead->core.n_cigar-1,*cigar=cigarBegin;cigar<=cigarEnd;cigar++)
            {
                hts_pos_t opLen=hts_pos_t(bam_cigar_oplen(*cigar));

                switch(bam_cigar_op(*cigar))
                {
                case BAM_CSOFT_CLIP:

                    if(cigar==cigarBegin) {readBegin-=opLen; continue;}
                    if(cigar==cigarEnd) readEnd+=opLen;
                    continue;

                case BAM_CDEL:
                case BAM_CREF_SKIP:
                case BAM_CMATCH:
                case BAM_CEQUAL:
                case BAM_CDIFF:

                    readEnd+=opLen;
                    continue;
                }
            }

            readEnd--;

            static const char bamSeq2ASCII[][2]={ {'n','N'},{'a','A'},{'c','C'},{'n','N'},{'g','G'},{'n','N'},{'n','N'},{'n','N'},{'t','T'},{'n','N'},{'n','N'},{'n','N'},{'n','N'},{'n','N'},{'n','N'},{'n','N'} };

            maxPosScore=0;

            begin=readBegin; end=readEnd;

            size_t queryLen=size_t(bamRead->core.l_qseq); uint8_t * seq=bam_get_seq(bamRead); uint8_t * qual=bam_get_qual(bamRead); uint8_t * qualEnd=qual+(queryLen&~1UL); char * pQuery=query=(char*)malloc((queryLen+1)*sizeof(char));

            for(;qual<qualEnd;seq++,qual+=2,pQuery+=2)
            {
                uint8_t packedBases=seq[0]; size_t isHQ;
                isHQ=size_t(qual[0]>=minHQBaseScore); maxPosScore+=1+isHQ; pQuery[0]=bamSeq2ASCII[size_t(packedBases>>4)][isHQ];
                isHQ=size_t(qual[1]>=minHQBaseScore); maxPosScore+=1+isHQ; pQuery[1]=bamSeq2ASCII[size_t(packedBases&15)][isHQ];
            }

            if(queryLen&1){ size_t isHQ=size_t(qual[0]>=minHQBaseScore); maxPosScore+=1+isHQ; pQuery[0]=bamSeq2ASCII[size_t(seq[0]>>4)][isHQ]; }

            query[queryLen]='\0';
        }
    };

    class Score
    {
    private:
    public:

        bool inJobSet;
        int reads[2];

        Score(void) : reads{0,0} {}
    };

    class Fragment
    {
    private:
    public:

        Read reads[2];
        map<VarEntry*,Score> scores;

        void CountUnknown(size_t fileIndex)
        {
            for(auto itScore=scores.begin(),scoresEnd=scores.end();itScore!=scoresEnd;itScore++)
            {
                if(itScore->second.inJobSet==false) continue;
                auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.unknown++;
            }
        }

        void CountRead(size_t fileIndex,Reference * maxRef,size_t strand)
        {
            uint32_t isHQ=uint32_t(reads[strand].isHQ);

            for(auto itScore=scores.begin(),scoresEnd=scores.end();itScore!=scoresEnd;itScore++)
            {
                if(itScore->second.inJobSet==false) continue;
                auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.hqDepth+=isHQ;
            }

            if(maxRef->inJobSet==false) return;
            auto & statistics=maxRef->varEntry->statistics[fileIndex]; statistics.altDepth++; statistics.altBias[strand]++; statistics.hqAltDepth+=isHQ; statistics.hqAltBias[strand]+=isHQ;
        }

        void CountFragment(size_t fileIndex,Reference * maxRef)
        {
            uint32_t isHQ=uint32_t(reads[0].isHQ)|uint32_t(reads[1].isHQ);

            for(auto itScore=scores.begin(),scoresEnd=scores.end();itScore!=scoresEnd;itScore++)
            {
                if(itScore->second.inJobSet==false) continue;
                auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.hqDepth+=isHQ;
            }

            if(maxRef->inJobSet==false) return;
            auto & statistics=maxRef->varEntry->statistics[fileIndex]; statistics.altDepth++; statistics.altBias[0]++; statistics.altBias[1]++; statistics.hqAltDepth+=isHQ; statistics.hqAltBias[0]+=isHQ; statistics.hqAltBias[1]+=isHQ;
        }

        uint32_t CountAmbigousRead(size_t fileIndex,size_t strand)
        {
            int maxScore=0; size_t nMaxScore=0; map<VarEntry*,Score>::iterator itMaxScore;

            auto scoresEnd=scores.end(); for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
            {
                int score=itScore->second.reads[strand];
                if(score>maxScore){ maxScore=score; nMaxScore=0;}
                if(score==maxScore && itScore->first->statistics[fileIndex].altDepth>0 && ++nMaxScore==1) itMaxScore=itScore;
            }

            if(nMaxScore!=1)
            {
                for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
                {
                    if(itScore->second.inJobSet==false) continue;
                    auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.ambiguous+=itScore->second.reads[strand]==maxScore;
                }
                return 0U;
            }

            uint32_t isHQ=uint32_t(reads[strand].isHQ);

            for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
            {
                if(itScore->second.inJobSet==false) continue;
                auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.hqDepth+=isHQ;
            }

            if(itMaxScore->second.inJobSet==true){ auto & statistics=itMaxScore->first->statistics[fileIndex]; statistics.altDepth++; statistics.altBias[strand]++; statistics.hqAltDepth+=isHQ; statistics.hqAltBias[strand]+=isHQ; }
            return 1U;
        }

        uint32_t CountAmbigousFragment(size_t fileIndex)
        {
            int fwMaxScore=0,rvMaxScore=0; size_t fwNMaxScore=0,rvNMaxScore=0; map<VarEntry*,Score>::iterator fwITMaxScore,rvITMaxScore;

            auto scoresEnd=scores.end(); for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
            {
                int fwScore=itScore->second.reads[0],rvScore=itScore->second.reads[1];
                if(fwScore>fwMaxScore){ fwMaxScore=fwScore; fwNMaxScore=0;}
                if(rvScore>rvMaxScore){ rvMaxScore=rvScore; rvNMaxScore=0;}

                if(itScore->first->statistics[fileIndex].altDepth>0)
                {
                    if(fwScore==fwMaxScore && ++fwNMaxScore==1) fwITMaxScore=itScore;
                    if(rvScore==rvMaxScore && ++rvNMaxScore==1) rvITMaxScore=itScore;
                }
            }

            uint32_t isHQ=uint32_t(reads[0].isHQ)|uint32_t(reads[1].isHQ);

            if(fwNMaxScore!=1)
            {
                if(rvNMaxScore!=1)
                {
                    for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
                    {
                        if(itScore->second.inJobSet==false) continue;
                        auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; auto reads=itScore->second.reads; statistics.ambiguous+=reads[0]==fwMaxScore || reads[1]==rvMaxScore;
                    }

                    return 0U;
                }

                for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
                {
                    if(itScore->second.inJobSet==false) continue;
                    auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.hqDepth+=isHQ;
                }

                if(rvITMaxScore->second.inJobSet==true){ auto & statistics=rvITMaxScore->first->statistics[fileIndex]; statistics.altDepth++; statistics.altBias[1]++; statistics.hqAltDepth+=isHQ; statistics.hqAltBias[1]+=isHQ; }
                return 1U;
            }

            if(rvNMaxScore!=1)
            {
                for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
                {
                    if(itScore->second.inJobSet==false) continue;
                    auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.hqDepth+=isHQ;
                }

                if(fwITMaxScore->second.inJobSet==true){ auto & statistics=fwITMaxScore->first->statistics[fileIndex]; statistics.altDepth++; statistics.altBias[0]++; statistics.hqAltDepth+=isHQ; statistics.hqAltBias[0]+=isHQ; }
                return 1U;
            }

            if(fwITMaxScore->first!=rvITMaxScore->first)
            {
                for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
                {
                    if(itScore->second.inJobSet==false) continue;
                    auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; auto reads=itScore->second.reads; statistics.ambiguous+=reads[0]==fwMaxScore || reads[1]==rvMaxScore;
                }
                return 0U;
            }

            for(auto itScore=scores.begin();itScore!=scoresEnd;itScore++)
            {
                if(itScore->second.inJobSet==false) continue;
                auto & statistics=itScore->first->statistics[fileIndex]; statistics.totalDepth++; statistics.hqDepth+=isHQ;
            }

            if(fwITMaxScore->second.inJobSet==true){auto & statistics=fwITMaxScore->first->statistics[fileIndex]; statistics.altDepth++; statistics.altBias[0]++; statistics.altBias[1]++; statistics.hqAltDepth+=isHQ; statistics.hqAltBias[0]+=isHQ; statistics.hqAltBias[1]+=isHQ;}
            return 1U;
        }
    };

    map<string,Fragment> fragments;
    vector<pair<hts_pos_t,hts_pos_t> > refIntervals; refIntervals.reserve(16);
    {
        bam1_t * bamRead=bam_init1();

        if(bamFile.ReadRegion(bamRead)>=0)
        {
            hts_pos_t intervalBegin,intervalEnd; fragments[bam_get_qname(bamRead)].reads[bam_is_rev(bamRead)].Init(bamRead,minHQAlignmentScore,minHQBaseScore,intervalBegin,intervalEnd);

            while(bamFile.ReadRegion(bamRead)>=0)
            {
                hts_pos_t bamReadBegin,bamReadEnd; fragments[bam_get_qname(bamRead)].reads[bam_is_rev(bamRead)].Init(bamRead,minHQAlignmentScore,minHQBaseScore,bamReadBegin,bamReadEnd);

                if(bamReadBegin<=intervalEnd+1)
                {
                    intervalEnd=max(intervalEnd,bamReadEnd);
                    continue;
                }

                refIntervals.emplace_back(intervalBegin,intervalEnd);
                intervalBegin=bamReadBegin; intervalEnd=bamReadEnd;
            }

            refIntervals.emplace_back(intervalBegin,intervalEnd);
        }

        bam_destroy1(bamRead);
    }

    if(fragments.size()==0) return;

    size_t nRefIntervals=refIntervals.size();

    //----------------------------------------------------------------
    //Create reference sequences
    //----------------------------------------------------------------

    vector<Reference> allReferences; allReferences.reserve(1024);
    multimap<hts_pos_t,Reference*> references;  //References by position
    {
        for(size_t i=0;i<nRefIntervals;i++)
        {
            hts_pos_t intervalBegin=refIntervals[i].first;
            hts_pos_t intervalEnd=refIntervals[i].second;

            string refSeq(fastaEntry.seq+intervalBegin,1+intervalEnd-intervalBegin);
            allReferences.emplace_back(true,refSeq,new VarEntry(nBamFiles));
        }

        set<VarEntry*> jobSet; for(auto entry : job.pos->second) if(jobPos1<=entry.first) jobSet.insert(entry.second);   //Jobset contains only the in job variants in normal orientation

        for(size_t i=0;i<nRefIntervals;i++)
        {
            hts_pos_t intervalBegin=refIntervals[i].first;
            hts_pos_t intervalEnd=refIntervals[i].second;

            //----------------------------------------------------------------
            //Iterate over variants in the interval
            //----------------------------------------------------------------

            auto & varEntriesByChr=varEntries.at(job.chr);

            for(auto var=varEntriesByChr.lower_bound(intervalBegin),varEnd=varEntriesByChr.upper_bound(intervalEnd);var!=varEnd;var++)
            {
                hts_pos_t pos1=var->first;

                for(auto & entry : var->second)
                {
                    hts_pos_t pos2=entry.first;

                    auto varEntry=entry.second;

                    bool inJobSet=jobSet.count(varEntry)==1; if(i>0 && inJobSet==true) continue;

                    hts_pos_t refSize=hts_pos_t(varEntry->ref.size()); string alt(varEntry->alt); hts_pos_t altSize=hts_pos_t(alt.size());

                    //----------------------------------------------------------------
                    //Variant is SNV
                    //----------------------------------------------------------------

                    if(varEntry->varType==SNV)
                    {
                        if(references.count(pos1)==0) references.emplace(pos1,&allReferences[i]);

                        string altSeq(fastaEntry.seq+intervalBegin,1+intervalEnd-intervalBegin); altSeq[pos1-intervalBegin]=alt[0];
                        references.emplace(pos1,&allReferences.emplace_back(inJobSet,altSeq,varEntry));

                        continue;
                    }

                    //----------------------------------------------------------------
                    //Variant is multibase substitution or deletion
                    //----------------------------------------------------------------

                    if(varEntry->varType==DEL)
                    {
                        if(pos1<pos2)
                        {
                            if(pos1==intervalEnd) continue; //Deletion has no effect

                            auto pReference=&allReferences[i];
                            if(references.count(pos1)==0) references.emplace(pos1,pReference);
                            if(references.count(pos2)==0) references.emplace(pos2,pReference);

                            hts_pos_t padBegin=min(pos1,intervalBegin)-max(min(pos2-1,jobPos1-1)-max(pos1,intervalBegin)+1L-altSize,0L);
                            hts_pos_t padEnd=max(pos2-1,intervalEnd)+max(min(pos2-1,intervalEnd)-max(pos1,jobPos2)+1L-hts_pos_t(pos1>=jobPos2)*altSize,0L);

                            string altSeq(fastaEntry.seq+padBegin,1+padEnd-padBegin); altSeq.replace(pos1-padBegin,refSize,alt);
                            pReference=&allReferences.emplace_back(inJobSet,altSeq,varEntry);
                            references.emplace(pos1,pReference); references.emplace(pos2,pReference);

                            continue;
                        }

                        if(pos1==intervalBegin || pos2>=intervalBegin) continue;    //Deletion has no effect or already in reference set

                        auto pReference=&allReferences[i];
                        if(references.count(pos1)==0) references.emplace(pos1,pReference);
                        if(references.count(pos2)==0) references.emplace(pos2,pReference);

                        hts_pos_t padBegin=pos2-max(min(pos1-1,jobPos1-1)-intervalBegin+1L-altSize,0L);
                        hts_pos_t padEnd=max(pos1-1,intervalEnd)+max(min(pos1-1,intervalEnd)-max(pos2,jobPos2)+1L-hts_pos_t(pos1>=jobPos2)*altSize,0L);

                        string altSeq(fastaEntry.seq+padBegin,1+padEnd-padBegin); altSeq.replace(pos1-padBegin,refSize,alt);
                        pReference=&allReferences.emplace_back(inJobSet,altSeq,varEntry);
                        references.emplace(pos1,pReference); references.emplace(pos2,pReference);

                        continue;
                    }

                    //----------------------------------------------------------------
                    //Variant is insertion or tandem duplication
                    //----------------------------------------------------------------

                    if(pos1<pos2)
                    {
                        if(pos1==intervalEnd) continue; //Insertion has no effect

                        auto pReference=&allReferences[i];
                        if(references.count(pos1)==0) references.emplace(pos1,pReference);
                        if(references.count(pos2)==0) references.emplace(pos2,pReference);

                        string altSeq(fastaEntry.seq+intervalBegin,1+intervalEnd-intervalBegin); altSeq.replace(pos1-intervalBegin,refSize,alt);
                        pReference=&allReferences.emplace_back(inJobSet,altSeq,varEntry);
                        references.emplace(pos1,pReference); references.emplace(pos2,pReference);

                        continue;
                    }
                }
            }
        }
    }

    //----------------------------------------------------------------
    //Align and count first pass
    //----------------------------------------------------------------

    uint32_t assigned=0;

    auto fragmentsEnd=fragments.end();

    for(auto itFragment=fragments.begin();itFragment!=fragmentsEnd;itFragment++)
    {
        Reference *fwMaxRef=nullptr,*rvMaxRef=nullptr;

        auto & fragment=itFragment->second;
        auto & fwRead=fragment.reads[0]; auto & rvRead=fragment.reads[1];

        if(fwRead.query!=nullptr)
        {
            int maxScore=0; size_t nMaxScore=0;

            for(auto itReference=references.lower_bound(fwRead.begin),referencesEnd=references.upper_bound(fwRead.end);itReference!=referencesEnd;itReference++)
            {
                auto reference=itReference->second; auto & scorePair=fragment.scores[reference->varEntry]; auto & score=scorePair.reads[0]; if(score!=0) continue;

                scorePair.inJobSet=reference->inJobSet;
                score=reference->ssw.Align(fwRead.query);

                if(score==maxScore) {nMaxScore++; continue;}
                if(score>maxScore) {fwMaxRef=reference; maxScore=score; nMaxScore=1;}
            }

            fwRead.isUnknown=double(maxScore)/double(fwRead.maxPosScore)<minAlignmentRate;
            fwRead.isAmbiguous=nMaxScore>1;

            if(rvRead.query==nullptr)
            {
                if(fwRead.isUnknown) {fragment.CountUnknown(job.fileIndex); continue;}
                if(fwRead.isAmbiguous) continue;
                fragment.CountRead(job.fileIndex,fwMaxRef,0); assigned++; continue;
            }
        }

        int maxScore=0; size_t nMaxScore=0;

        for(auto itReference=references.lower_bound(rvRead.begin),referencesEnd=references.upper_bound(rvRead.end);itReference!=referencesEnd;itReference++)
        {
            auto reference=itReference->second; auto & scorePair=fragment.scores[reference->varEntry]; auto & score=scorePair.reads[1]; if(score!=0) continue;

            scorePair.inJobSet=reference->inJobSet;
            score=reference->ssw.Align(rvRead.query);

            if(score==maxScore) {nMaxScore++; continue;}
            if(score>maxScore) {rvMaxRef=reference; maxScore=score; nMaxScore=1;}
        }

        rvRead.isUnknown=double(maxScore)/double(rvRead.maxPosScore)<minAlignmentRate;
        rvRead.isAmbiguous=nMaxScore>1;

        if(fwRead.query==nullptr)
        {
            if(rvRead.isUnknown) {fragment.CountUnknown(job.fileIndex); continue;}
            if(rvRead.isAmbiguous) continue;
            fragment.CountRead(job.fileIndex,rvMaxRef,1); assigned++; continue;
        }

        if(rvRead.isUnknown)
        {
            if(fwRead.isUnknown) {fragment.CountUnknown(job.fileIndex); continue;}
            if(fwRead.isAmbiguous) continue;
            fragment.CountRead(job.fileIndex,fwMaxRef,0); assigned++; continue;
        }

        if(fwRead.isUnknown)
        {
            if(rvRead.isAmbiguous) continue;
            fragment.CountRead(job.fileIndex,rvMaxRef,1); assigned++; continue;
        }

        if(rvRead.isAmbiguous)
        {
            if(fwRead.isAmbiguous) continue;
            fragment.CountRead(job.fileIndex,fwMaxRef,0); assigned++; continue;
        }

        if(fwRead.isAmbiguous)
        {
            fragment.CountRead(job.fileIndex,rvMaxRef,1); assigned++; continue;
        }

        if(fwMaxRef->varEntry!=rvMaxRef->varEntry) continue;

        fragment.CountFragment(job.fileIndex,fwMaxRef); assigned++;
    }

    //----------------------------------------------------------------
    //Second pass count ambiguous
    //----------------------------------------------------------------

    for(auto itFragment=fragments.begin();itFragment!=fragmentsEnd;itFragment++)
    {
        auto & fragment=itFragment->second;
        auto & fwRead=fragment.reads[0]; auto & rvRead=fragment.reads[1];

        if(fwRead.query!=nullptr)
        {
            if(rvRead.query==nullptr)
            {
                if(fwRead.isUnknown) continue;
                if(fwRead.isAmbiguous) {assigned+=fragment.CountAmbigousRead(job.fileIndex,0);} continue;
            }
        }
        else
        {
            if(rvRead.isUnknown) continue;
            if(rvRead.isAmbiguous) {assigned+=fragment.CountAmbigousRead(job.fileIndex,1);} continue;
        }

        if(rvRead.isUnknown)
        {
            if(fwRead.isUnknown) continue;
            if(fwRead.isAmbiguous) {assigned+=fragment.CountAmbigousRead(job.fileIndex,0);} continue;
        }

        if(fwRead.isUnknown)
        {
            if(rvRead.isAmbiguous) {assigned+=fragment.CountAmbigousRead(job.fileIndex,1);} continue;
        }

        if(rvRead.isAmbiguous)
        {
            if(fwRead.isAmbiguous) assigned+=fragment.CountAmbigousFragment(job.fileIndex);
        }
    }

    //----------------------------------------------------------------
    //Reestimate ITD/PTD VAF
    //----------------------------------------------------------------

    for(auto entry=job.pos->second.begin(),entryEnd=job.pos->second.end();entry!=entryEnd;entry++)
    {
        auto varEntry=entry->second; auto varType=varEntry->varType;

        if(varType!=ITD && varType!=PTD) continue;

        auto & statistics=varEntry->statistics[job.fileIndex];

        double MUTi=double(statistics.altDepth);
        double nonMUTi=max(0.0,0.5*double(statistics.ambiguous)-MUTi);
        double fAssigned=double(assigned);

        statistics.vafITD=MUTi/(nonMUTi+fAssigned);
        statistics.coeITD=nonMUTi/(nonMUTi+fAssigned-MUTi);
    }

    //----------------------------------------------------------------
    //Clean up ref intervals
    //----------------------------------------------------------------

    for(size_t i=0;i<nRefIntervals;i++) delete allReferences[i].varEntry;

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
}
//----------------------------------------------------------------
void Thread::RunJob(Job & job)
{
    //----------------------------------------------------------------
    //Open bam file
    //----------------------------------------------------------------

    bamFile.Open(bamFileNames[job.fileIndex],countDuplicates,countSecondary);

    //----------------------------------------------------------------
    //Determine if the locus is second position only or SNP only
    //----------------------------------------------------------------

    bool isSNVOnly=true;

    for(auto & entry : job.pos->second)
    {
        if(entry.second->varType!=SNV) isSNVOnly=false;
    }

    //----------------------------------------------------------------
    //Variant is SNV
    //----------------------------------------------------------------

    if(isSNVOnly){ProcessSNP(job); return;}

    //----------------------------------------------------------------
    //Variant is insertion
    //----------------------------------------------------------------

    ProcessIndel(job);
}
//----------------------------------------------------------------
#define FASTA_FILE                      'f'
#define VCF_FILE                        'V'
#define ANNOVAR_FILE                    'a'
#define VEP_FILE                        'e'
#define BAM_FILES                       'b'
#define MIN_MATCH_LENGTH                'm'
#define PILEUP_TOLERANCE                'T'
#define MIN_HQ_BASE_SCORE               's'
#define MIN_HQ_ALIGNMENT_SCORE          'S'
#define MIN_ALIGNMENT_RATE              'r'
#define THREADS                         't'
#define COUNT_DUPLICATES                'd'
#define COUNT_SECONDARY                 'u'
#define VERBOSE                         'v'
#define HELP                            'h'
#define SHORT_OPTIONS                   "f:V:a:e:b:m:T:s:S:r:t:duvh"
//----------------------------------------------------------------
struct option longOptions[] =
{
    {"fasta-file",required_argument,nullptr,FASTA_FILE},
    {"vcf-file",required_argument,nullptr,VCF_FILE},
    {"annovar-file",required_argument,nullptr,ANNOVAR_FILE},
    {"vep-file",required_argument,nullptr,VEP_FILE},
    {"bam-files",required_argument,nullptr,BAM_FILES},
    {"min-match-length",required_argument,nullptr,MIN_MATCH_LENGTH},
    {"pileup-tolerance",required_argument,nullptr,PILEUP_TOLERANCE},
    {"min-hq-base-score",required_argument,nullptr,MIN_HQ_BASE_SCORE},
    {"min-hq-alignment-score",required_argument,nullptr,MIN_HQ_ALIGNMENT_SCORE},
    {"min-alignment-rate",required_argument,nullptr,MIN_ALIGNMENT_RATE},
    {"threads",required_argument,nullptr,THREADS},
    {"count-duplicates",no_argument,nullptr,COUNT_DUPLICATES},
    {"count-secondary",no_argument,nullptr,COUNT_SECONDARY},
    {"verbose",no_argument,nullptr,VERBOSE},
    {"help",no_argument,nullptr,HELP},
    {0, 0, 0, 0}
};
//----------------------------------------------------------------
void AnnotateBamStatistics::Tokenize(char * s,const char * d,vector<string> & v)
{
    v.clear(); for(char * t=strsep(&s,d);t!=nullptr;t=strsep(&s,d)) v.push_back(t);
}
//----------------------------------------------------------------
int AnnotateBamStatistics::Run(int argc,char * argv[])
{
    try
    {
        //----------------------------------------------------------------
        //Get input arguments
        //----------------------------------------------------------------

        bool showHelp=(argc==1);
        int option,optionIndex;

        bool verbose=true;
        bool countDuplicates=false;
        bool countSecondary=false;

        uint8_t minHQBaseScore=30;
        uint8_t minHQAlignmentScore=40;

        int nThreads=1;

        hts_pos_t minMatchLength=15;
        hts_pos_t pileupTolerance=5;

        double minAlignmentRate=0.90;

        string fastaFilename,vcfFilename,annovarFilename,vepFilename;
        vector<string> bamFilenames;

        verbose=false;

        while((option=getopt_long(argc,argv,SHORT_OPTIONS,longOptions,&optionIndex))>=0)
        {
            switch(option)
            {
            case FASTA_FILE: fastaFilename=string(optarg); break;
            case VCF_FILE: vcfFilename=string(optarg); break;
            case ANNOVAR_FILE: annovarFilename=string(optarg); break;
            case VEP_FILE: vepFilename=string(optarg); break;
            case BAM_FILES: Tokenize(optarg,",",bamFilenames); break;
            case MIN_MATCH_LENGTH: minMatchLength=atoll(optarg); break;
            case PILEUP_TOLERANCE: pileupTolerance=atoll(optarg); break;
            case MIN_HQ_BASE_SCORE: minHQBaseScore=atoi(optarg); break;
            case MIN_HQ_ALIGNMENT_SCORE: minHQAlignmentScore=atoi(optarg); break;
            case MIN_ALIGNMENT_RATE: minAlignmentRate=atof(optarg); break;
            case THREADS: nThreads=atoi(optarg); break;
            case COUNT_DUPLICATES: countDuplicates=true; break;
            case COUNT_SECONDARY: countSecondary=true; break;
            case VERBOSE: verbose=true; break;
            case HELP:
            default:

                showHelp=true;
                break;
            }
        }

        //----------------------------------------------------------------
        //Show help
        //----------------------------------------------------------------

        if(verbose) cerr << "Info: Show help?" << endl;

        if(showHelp)
        {
            cerr << "enhanced_ABS v2.0 [options] > output.txt (output always to std::out)"                                                              << endl;
            cerr                                                                                                                                        << endl;
            cerr << "-f --fasta-file <text>             Single fasta file(required)"                                                                    << endl;
            cerr << "-V --vcf-file <text>               Single vcf file (optional)"                                                                     << endl;
            cerr << "-a --annovar-file <text>           Single annovar file (required either annovar-file or vep-file)"                                 << endl;
            cerr << "-e --vep-file <text>               Single VEP file (required either annovar-file or vep-file)"                                     << endl;
            cerr << "-b --bam-files <text>              One or more bam files (required)"                                                               << endl;
            cerr << "-m --min-match-length <int>        Minimum match length (optional default=15)"                                                     << endl;
            cerr << "-T --pileup-tolerance <int>        Pileup tolerance number of extra bases around the variant position (optional default=5)"        << endl;
            cerr << "-s --min-hq-base-score <int>       Minimum base score for filtered statistics (optional,default=30)"                               << endl;
            cerr << "-S --min-hq-alignment-score <int>  Minimum alignment score for filtered statistics (optional,default=40)"                          << endl;
            cerr << "-r --min-alignment-rate <float>    Minimum Smith-Waterman alignment rate (optional,default=0.95)"                                  << endl;
            cerr << "-t --threads <int>                 Number of threads to use (optional default=1)"                                                  << endl;
            cerr << "-d --count-duplicates <void>       If specified duplicates fragments are used in the statistics (optional default=false)"          << endl;
            cerr << "-u --count-secondary <void>        If specified secondary fragments are used in the statistics (optional default=false)"           << endl;
            cerr << "-v --verbose <void>                If specified be verbose (optional default=false)"                                               << endl;
            cerr                                                                                                                                        << endl;

            return 0;
        }

        //----------------------------------------------------------------
        //Check input arguments
        //----------------------------------------------------------------

        if(verbose) cerr << "Info: Check input arguments" << endl;

        if(fastaFilename.empty())
            throw runtime_error("Error: Please specify a single fasta file (-h for help)");

        if(int(annovarFilename.empty())+int(vepFilename.empty())!=1)
            throw runtime_error("Error: Please specify either an annovar file or vep file not both! (-h for help)");

        size_t nBamFiles=bamFilenames.size();

        if(nBamFiles==0)
            throw runtime_error("Error: Please specify at least one bam file! (-h for help)");

        if(minMatchLength<10) minMatchLength=10;
        if(pileupTolerance<0) pileupTolerance=0;
        if(nThreads<1) nThreads=1;

        minAlignmentRate=min(max(minAlignmentRate,0.2),1.0);

        //----------------------------------------------------------------
        //Open fasta file
        //----------------------------------------------------------------

        if(verbose) cerr << "Info: Open fasta file" << endl;
        FastaFile fastaFile(fastaFilename);

        //----------------------------------------------------------------
        //Open vcf file if specified
        //----------------------------------------------------------------

        VCFFile vcfFile;
        if(vcfFilename.empty()==false)
        {
            if(verbose) cerr << "Info: Open vcf file" << endl;
            vcfFile.Open(vcfFilename);
        }

        //----------------------------------------------------------------
        //Use annovar file
        //----------------------------------------------------------------

        if(annovarFilename.empty()==false)
        {
            if(verbose) cerr << "Info: Open annovar file" << endl;
            AnnovarFile annovarFile(annovarFilename,fastaFile.entries,nBamFiles);

            //----------------------------------------------------------------
            //Create jobs
            //----------------------------------------------------------------

            if(verbose) cerr << "Info: Create jobs" << endl;

            vector<Job> jobs;

            for(size_t fileIndex=0;fileIndex<nBamFiles;fileIndex++)
            {
                for(auto entry=annovarFile.entries.begin(),entriesEnd=annovarFile.entries.end();entry!=entriesEnd;entry++)
                {
                    for(auto pos=entry->second.begin(),posEnd=entry->second.end();pos!=posEnd;pos++)
                    {
                        bool isPos2Only=true;

                        hts_pos_t pos1=pos->first;

                        for(auto varEntry : pos->second)
                        {
                            if(pos1<=varEntry.first) isPos2Only=false;
                        }

                        if(isPos2Only==true) continue;

                        jobs.emplace_back(fileIndex,entry->first.c_str(),pos);
                    }
                }
            }

            size_t nJobs=jobs.size();

            //----------------------------------------------------------------
            //Pileup variants
            //----------------------------------------------------------------

            if(verbose)
            {
                cerr << "Info: Pileup variants" << endl;
                cerr << progressBar[0] << flush;
            }

            bool errorOccured=false;
            size_t prevProgress=0;

            omp_set_num_threads(nThreads);

            #pragma omp parallel
            {
                Thread thread(countDuplicates,countSecondary,minHQBaseScore,minHQAlignmentScore,minMatchLength,pileupTolerance,nBamFiles,minAlignmentRate,bamFilenames,fastaFile.entries,vcfFile.entries,annovarFile.entries);

                #pragma omp for schedule(dynamic)
                for(size_t i=0;i<nJobs;i++)
                {
                    try
                    {
                        thread.RunJob(jobs[i]);

                        if(verbose)
                        {
                            size_t progress=(i*100)/nJobs;
                            #pragma omp critical
                            {
                                if(progress>prevProgress)
                                {
                                    prevProgress=progress;
                                    cerr << progressBar[progress] << flush;
                                }
                            }
                        }
                    }
                    catch(const runtime_error & error)
                    {
                        errorOccured=true;
                        cerr << error.what() << endl;
                    }
                }
            }

            if(errorOccured) return 1;
            if(verbose) cerr << progressBar[100] << endl;

            //----------------------------------------------------------------
            //Output resulting annovar file
            //----------------------------------------------------------------

            if(verbose) cerr << "Info: Output results" << endl;
            annovarFile.Write(bamFilenames);

            //----------------------------------------------------------------
            //Done
            //----------------------------------------------------------------

            if(verbose) cerr << "Info: Done" << endl;
            return 0;
        }

        //----------------------------------------------------------------
        //Open VEP file
        //----------------------------------------------------------------

        if(verbose) cerr << "Info: Open VEP file" << endl;
        VEPFile vepFile(vepFilename,fastaFile.entries,nBamFiles);

        //----------------------------------------------------------------
        //Create jobs
        //----------------------------------------------------------------

        if(verbose) cerr << "Info: Create jobs" << endl;

        vector<Job> jobs;

        for(size_t fileIndex=0;fileIndex<nBamFiles;fileIndex++)
        {
            for(auto entry=vepFile.entries.begin(),entriesEnd=vepFile.entries.end();entry!=entriesEnd;entry++)
            {
                for(auto pos=entry->second.begin(),posEnd=entry->second.end();pos!=posEnd;pos++)
                {
                    bool isPos2Only=true;

                    hts_pos_t pos1=pos->first;

                    for(auto varEntry : pos->second)
                    {
                        if(pos1<=varEntry.first) isPos2Only=false;
                    }

                    if(isPos2Only==true) continue;

                    jobs.emplace_back(fileIndex,entry->first.c_str(),pos);
                }
            }
        }

        size_t nJobs=jobs.size();

        //----------------------------------------------------------------
        //Pileup variants
        //----------------------------------------------------------------

        if(verbose)
        {
            cerr << "Info: Pileup variants" << endl;
            cerr << progressBar[0] << flush;
        }

        bool errorOccured=false;
        size_t prevProgress=0;

        omp_set_num_threads(nThreads);

        #pragma omp parallel
        {
            Thread thread(countDuplicates,countSecondary,minHQBaseScore,minHQAlignmentScore,minMatchLength,pileupTolerance,nBamFiles,minAlignmentRate,bamFilenames,fastaFile.entries,vcfFile.entries,vepFile.entries);

            #pragma omp for schedule(dynamic)
            for(size_t i=0;i<nJobs;i++)
            {
                try //Catch errors within the same thread
                {
                    thread.RunJob(jobs[i]);

                    if(verbose)
                    {
                        size_t progress=(i*100)/nJobs;
                        #pragma omp critical
                        {
                            if(progress>prevProgress)
                            {
                                prevProgress=progress;
                                cerr << progressBar[progress] << flush;
                            }
                        }
                    }
                }
                catch(const runtime_error & error)
                {
                    errorOccured=true;
                    cerr << error.what() << endl;
                }
            }
        }

        if(errorOccured) return 1;
        if(verbose) cerr << progressBar[100] << endl;

        //----------------------------------------------------------------
        //Output resulting annovar file
        //----------------------------------------------------------------

        if(verbose) cerr << "Info: Output results" << endl;
        vepFile.Write(bamFilenames);

        //----------------------------------------------------------------
        //Done
        //----------------------------------------------------------------

        if(verbose) cerr << "Info: Done" << endl;
        return 0;

    }
    catch(const runtime_error & error)
    {
        cerr << error.what() << endl;
        return 1;
    }
}
//----------------------------------------------------------------


