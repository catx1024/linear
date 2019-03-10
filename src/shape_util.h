#ifndef LINEAR_HEADER_SHAPE_EXTEND_H
#define LINEAR_HEADER_SHAPE_EXTEND_H
#include <utility>
#include <seqan/sequence.h>
using namespace seqan;

class LShape
{
public:
    typedef uint64_t THashValue;

    unsigned span;
    unsigned weight;
    THashValue hValue;      //hash value 
    THashValue crhValue;    //ReverseComplement
    THashValue XValue;     //minimizer 
    THashValue YValue;     //Y(h,x)
    THashValue strand;
    int leftChar;
    int x;

    LShape();
    LShape(unsigned);
    LShape(unsigned, unsigned);
};

typedef typename Iterator<String<Dna5> >::Type TIter_S;
uint64_t hashInit(LShape & me, TIter_S const &it);
uint64_t hashInit_hs(LShape & me, TIter_S const &it, int d = 0);
uint64_t hashNext(LShape & me, TIter_S const &it);
uint64_t hashNexth(LShape & me, TIter_S const &it);
uint64_t hashNext_hs(LShape &me, TIter_S const &it);
uint64_t hashPre_hs(LShape & me, TIter_S const &it);
uint64_t hashNextV(LShape & me, TIter_S const &it);
uint64_t hashNextX(LShape & me, TIter_S const &it);
#endif