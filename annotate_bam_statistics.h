//----------------------------------------------------------------
#ifndef ANNOTATE_BAM_STATISTICS_H
#define ANNOTATE_BAM_STATISTICS_H
//----------------------------------------------------------------
#include <vector>
#include <string>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class AnnotateBamStatistics
{
private:

   static void Tokenize(char * s,const char * d,vector<string> & v);

public:

    static int Run(int argc,char * argv[]);
};
//----------------------------------------------------------------
#endif // ANNOTATE_BAM_STATISTICS_H
