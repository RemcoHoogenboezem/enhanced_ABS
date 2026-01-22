//----------------------------------------------------------------
#ifndef SSW_H
#define SSW_H
//----------------------------------------------------------------
#include <string>
#include <parasail.h>
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class SSW
{
private:

    parasail_profile_t * referenceProfile;

    void Clear(void);

public:

    SSW(void);
    SSW(const string & reference);
    SSW(size_t referenceLen,const char * reference);

    ~SSW(void);

    void Init(const string & reference);
    void Init(size_t referenceLen,const char * reference);

    int Align(const char * query);
    int Align(const string & query);
    int Align(size_t queryLen,const char * query);
};
//----------------------------------------------------------------
#endif // SSW_H
