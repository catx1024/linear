#ifndef LINEAR_HEADER_ALIGN_INT_H
#define LINEAR_HEADER_ALIGN_INT_H
#include <utility>
using namespace seqan

struct MiniValueBit
{
    enum{VALUEBIT = 64};
};

struct MiniHEX{               
    enum {HEX = 8 };
};

template <unsigned shapeLength>
struct MiniWeight{
    enum{ WEIGHT = shapeLength - MiniHEX::HEX};
};

template <unsigned TSPAN, unsigned TWEIGHT = MiniWeight<TSPAN>::WEIGHT, typename TSpec = void>
struct Minimizer;
typedef Minimizer<0, 0> SimpleMShape;

class LShape
{
public:
    typedef typename uint64_t THashValue;

    unsigned span;
    unsigned weight;
    THashValue hValue;      //hash value 
    THashValue crhValue;    //ReverseComplement
    THashValue XValue;     //minimizer 
    THashValue YValue;     //Y(h,x)
    THashValue strand;
    int  leftChar;
    int x;

    LShape();
    LShape(unsigned);
    LShape(unsigned, unsigned);
};

#endif