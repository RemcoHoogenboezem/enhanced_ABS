//----------------------------------------------------------------
// Name        : ssw.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : Striped smith waterman implementation around parasail (Swap reference and query in this implementation for performance)
//----------------------------------------------------------------
#include <stdexcept>
#include <cstring>
#include "ssw.h"
//----------------------------------------------------------------
static const int gapOpen=3;
static const int gapExtension=1;
//----------------------------------------------------------------
static const int matrix[]=  //Substitution matrix (Rows=query Cols=reference)
{
/*       a   c   g   t   A   C   G   T   * */
/* a */  1, -1, -1, -1,  1, -1, -1, -1,  0,
/* c */ -1,  1, -1, -1, -1,  1, -1, -1,  0,
/* g */ -1, -1,  1, -1, -1, -1,  1, -1,  0,
/* t */ -1, -1, -1,  1, -1, -1, -1,  1,  0,
/* A */  2, -1, -1, -1,  2, -1, -1, -1,  0,
/* C */ -1,  2, -1, -1, -1,  2, -1, -1,  0,
/* G */ -1, -1,  2, -1, -1, -1,  2, -1,  0,
/* T */ -1, -1, -1,  2, -1, -1, -1,  2,  0,
/* * */  0,  0,  0,  0,  0,  0,  0,  0,  0
};
//----------------------------------------------------------------
static const int mapper[256]=   //Alphabet mapper
{
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	4,	8,	5,	8,	8,	8,	6,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	7,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	0,	8,	1,	8,	8,	8,	2,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	3,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,
    8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8
};
//----------------------------------------------------------------
static const parasail_matrix_t parasailMatrix=
{
    "parasailMatrix",
    matrix,
    mapper,
    9,
    2,
    -1,
    nullptr,
    PARASAIL_MATRIX_TYPE_SQUARE,
    9,
    "acgtACGT*",
    nullptr
};
//----------------------------------------------------------------
void SSW::Clear(void)
{
    if(referenceProfile!=nullptr) parasail_profile_free(referenceProfile);
}
//----------------------------------------------------------------
SSW::SSW(void) : referenceProfile(nullptr)
{
}
//----------------------------------------------------------------
SSW::SSW(const string & reference) : referenceProfile(nullptr)
{
    Init(reference);
}
//----------------------------------------------------------------
SSW::SSW(size_t referenceLen,const char * reference) : referenceProfile(nullptr)
{
    Init(referenceLen,reference);
}
//----------------------------------------------------------------
SSW::~SSW(void)
{
    Clear();
}
//----------------------------------------------------------------
void SSW::Init(const string & reference)
{
    Clear();
    referenceProfile=parasail_profile_create_16(reference.c_str(),int(reference.length()),&parasailMatrix);
}
//----------------------------------------------------------------
void SSW::Init(size_t referenceLen,const char * reference)
{
    Clear();
    referenceProfile=parasail_profile_create_16(reference,referenceLen,&parasailMatrix);
}
//----------------------------------------------------------------
int SSW::Align(const char * query)
{
    if(referenceProfile==nullptr) throw runtime_error("Error: Please initialze the SSW object first before aligning!");

    parasail_result_t * parasailResult=parasail_sw_striped_profile_16(referenceProfile,query,strlen(query),gapOpen,gapExtension);
    int score=parasailResult->score;
    parasail_result_free(parasailResult);

    return score;
}
//----------------------------------------------------------------
int SSW::Align(const string & query)
{
    if(referenceProfile==nullptr) throw runtime_error("Error: Please initialze the SSW object first before aligning!");

    parasail_result_t * parasailResult=parasail_sw_striped_profile_16(referenceProfile,query.c_str(),query.length(),gapOpen,gapExtension);
    int score=parasailResult->score;
    parasail_result_free(parasailResult);

    return score;
}
//----------------------------------------------------------------
int SSW::Align(size_t queryLen,const char * query)
{
    if(referenceProfile==nullptr) throw runtime_error("Error: Please initialze the SSW object first before aligning!");

    parasail_result_t * parasailResult=parasail_sw_striped_profile_16(referenceProfile,query,queryLen,gapOpen,gapExtension);
    int score=parasailResult->score;
    parasail_result_free(parasailResult);

    return score;
}
//----------------------------------------------------------------




