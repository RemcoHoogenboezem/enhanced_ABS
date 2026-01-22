//----------------------------------------------------------------
// Name        : main.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <cstring>
#include "annotate_bam_statistics.h"
//----------------------------------------------------------------
int main(void)
{
    char command[]="enhanced_ABS -f /home/remco/Desktop/test_enhanced_abs/hg38.fa -a /home/remco/Desktop/6241.txt -b /home/remco/Desktop/6241.bam"; char * pCommand=command;
    int argc=7; char * argv[16]; for(int i=0;i<argc;i++) argv[i]=strsep(&pCommand," ");

    return AnnotateBamStatistics::Run(argc,argv);
}
//----------------------------------------------------------------
