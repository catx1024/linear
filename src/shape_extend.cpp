#include "shape_extend.h"

using namespace seqan

const float _boundAlpha = 0.8;
unsigned d_span = 8;

LShape::LShape():
        span(25),
        weight(17),
        hValue(0),
        crhValue(0),
        XValue(0),
        YValue(0),
        strand(0),
        leftChar(0),
        x(0)
{}
LShape::LShape(unsigned v1):
        span(v1),
        weight(v1 - d_span),
        hValue(0),
        crhValue(0),
        XValue(0),
        YValue(0),
        strand(0),
        leftChar(0),
        x(0)
{}
LShape::LShape(unsigned v1, unsigned v2):
        span(v1),
        weight(v2),
        hValue(0),
        crhValue(0),
        XValue(0),
        YValue(0),
        strand(0),
        leftChar(0),
        x(0)
{}

inline unsigned weight(LShape const &me)
{
    return me.weight;
}

// ----------------------------------------------------------------------------
// hashInit() & hashNext()
// ----------------------------------------------------------------------------
typedef Iterator<String<Dna5> >::Type TIter; 

template <unsigned span> 
struct MASK
{
    static const uint64_t VALUE = (1ULL << span) - 1;
};

static const uint64_t COMP4 = 3;
static const int  ordC = 3;

inline void hashInit(LShape & me, TIter const &it)
{

    SEQAN_ASSERT_GT((unsigned)me.span, 0u);

    me.leftChar = 0;
    me.hValue = 0;
    for (unsigned i = 0; i < me.span - 1; ++i)
    {
        me.hValue = (me.hValue << 2) + ordValue((TValue)*(it + i));
    }
    me.x = 0;
}

inline uint64_t hashInit(LShape & me, TIter const &it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);

    me.leftChar = 0;
    me.hValue = 0;
    me.crhValue = 0;
    me.leftChar = 0;
    me.x = me.leftChar- ordC;
    uint64_t k =0, count = 0; //COMP for complemnet value;

    while (count < me.span)
    {
        if (ordValue(*(it + k + count)) == 4)
        {
            k += count + 1;
            count = 0;
        }
        else
            count++;
    }
    unsigned bit = 2;
    for (unsigned i = 0; i < me.span - 1; ++i)
    {
        uint64_t val = ordValue (*(it + k + i));
        me.x += (val << 1) - ordC;
        me.hValue = (me.hValue << 2) + val;
        me.crhValue += ((COMP4 - val) << bit);
        bit += 2;
        
    }
    return k;
}
/**
 *init for hashNexthS 
 */
inline uint64_t hashInit_hs(LShape & me, TIter const &it, int d)
{
    me.hValue = 0;
    for (unsigned i = d; i < me.span - 1 + d; ++i)
    {
        me.hValue = (me.hValue << 2) + ordValue (*(it + i));
    }
    me.hValue <<= (d << 1);
    return 0;
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 */ 
inline uint64_t hashNext(LShape & me, TIter const &it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned t, span = me.span << 1, weight = me.weight << 1;
    uint64_t v2 = ordValue((uint64_t)*(it + me.span - 1));
    uint64_t mask = (1ULL << (span - 2)) - 1;
    me.hValue=((me.hValue & mask) << 2 )+ v2;
    me.crhValue=((me.crhValue >> 2) & mask) + ((COMP4 - v2) << (span - 2));
    me.XValue = mask; 
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue(*(it));
    //printf("[debug]::hash %d\n", me.x);
    //v2 = (me.x > 0)?me.hValue:me.crhValue;
    if (me.x > 0)
    {
        v2 =me.hValue; 
        me.strand = 0;
    }
    else 
    {
        v2 = me.crhValue;
        me.strand = 1;
    }
    for (unsigned k = 64-span; k <= 64 - weight; k+=2)
    {
        v1 = v2 << k >> (64-weight);
        if(me.XValue > v1)
        {
            me.XValue=v1;
            t = k;
        }
    } 
    
    me.YValue = (v2 >> (64-t) << (64-t-weight)) +
            (v2 & ((1ULL<<(64-t-weight)) - 1)) + 
            (t << (span - weight - 1));
    //me.YValue = 0;
    return me.XValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate hValue;
 */ 
inline uint64_t hashNexth(LShape & me, TIter const &it)
{
    //typedef typename Size< Shape<TValue, TSpec> >::Type  TSize;
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t mask = (1ULL << (me.span * 2 - 2)) - 1;
    uint64_t v2 = ordValue((uint64_t)*(it + me.span - 1 ));
    me.hValue=((me.hValue & mask) << 2)+ v2;
    me.crhValue=((me.crhValue >> 2) & mask) + 
                ((COMP4 - v2) << ((me.span << 1) - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue(*(it));
    return me.x; 
}
/**
 * only calculate hash value for single strand
 * calculate hValue;
 */ 
inline uint64_t hashNext_hs(LShape &me, TIter const &it)
{
    uint64_t v2 = ordValue((uint64_t)*(it + me.span - 1 ));
    uint64_t mask = (1ULL << (me.span * 2 - 2)) - 1;
    me.hValue=((me.hValue & mask)<< 2) + v2;
    return me.hValue; 
}

inline uint64_t hashPre_hs(LShape & me, TIter const &it)
{
    uint64_t v2 = ordValue((uint64_t)*(it)) << ((me.span << 1)  - 2);
    uint64_t mask = (1ULL << (me.span * 2 - 2)) - 1;
    me.hValue=((me.hValue >> 2) & mask)+ v2;
    return me.hValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate hValue;
 */ 
inline uint64_t hashNextV(LShape & me, TIter const &it)
{
    //typedef typename Size< Shape<TValue, TSpec> >::Type  TSize;
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t  v2 = ordValue((uint64_t)*(it + me.span - 1 ));
    uint64_t mask = (1ULL << (me.span * 2 - 2)) - 1;
    me.hValue=((me.hValue & mask) << 2)+ v2;
    me.crhValue=((me.crhValue >> 2) & mask) + 
                ((COMP4 - v2) << (me.span * 2 - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue(*(it));
    me.strand = (me.x >> 63) & 1; //Note: me.x type is uint64_t
    return (me.x > 0)?me.hValue:me.crhValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate XValue
 */ 
inline uint64_t hashNextX(LShape & me, TIter const &it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned span = me.span << 1, weight = me.weight << 1;
    uint64_t t, v2;
    if (me.x > 0)
    {
        v2 =me.hValue; 
        me.strand = 0;
    }
    else 
    {
        v2 = me.crhValue;
        me.strand = 1;
    }
    me.XValue = (1ULL << span) - 1;
    for (unsigned k = 64-span; k <= 64 - weight; k+=2)
    {
        v1 = v2 << k >> (64-weight);
        if(me.XValue > v1)
        {
            me.XValue=v1;
            t = k;
        }
    } 
    me.YValue = (v2 >> (64-t) << (64-t-weight)) +
                (v2 & ((1ULL<<(64-t-weight)) - 1)) + 
                (t << (span - weight - 1));
                
    (void)it;
    return me.XValue; 
}
inline uint64_t getT(LShape & me)
{
    return (me.YValue >> ((me.span - me.weight) << 1));
}
inline uint64_t h2y(LShape & me, uint64_t h)
{
    uint64_t x = -1, v1, t=0;
    for (unsigned k = 64-(me.span << 1) ; k <= 64 - (me.weight << 1); k+=2)
    {
        v1 = h << k >> (64-(me.weight<<1));
        if(x > v1)
        { 
            x=v1;
            t=k;
        }
    } 
    return (h >> (64 - t) << (64 - t - (me.weight << 1))) + (h & (((uint64_t)1 << (64 - t - (me.weight << 1))) - 1))+(t << (((me.span - me.weight) << 1) - 1));
}

