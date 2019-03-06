#ifndef LINEAR_HEADER_PMP_FINDER_H
#define LINEAR_HEADER_PMP_FINDER_H

using namespace seqan;

struct CordBase
{
    //Cord(C): coordinates of the vertex of sliding windows
    //=|N/A[2]|strand[1]|cordEnd[1] genomeCord [40] |readCord [20bits]
    //cell [4] is the minimum length the window is allowed to slide in the alignment matrix.
    //genomeCord(gC or xC): = position in the genome >> cell_bit << cell_bit. the last cell_bit bits maybe set to 0
    //gC:= SA node = Seq num i1 [10] | Base num i2 [30]  
    //readCord(rC or yC): ~= position in the read >> cell_bit << cell_bit. the last cell_bit bits maybe set to 0 during process.
    //rC:= Base num [20]
    
    typedef unsigned Bit;
    typedef uint64_t Mask;
    typedef uint64_t Flag;
    typedef uint64_t CordType;
    typedef uint64_t CellType;
    typedef unsigned Size;
    
    Bit bit;
    uint64_t flagEnd;
    Mask mask;
    Mask maskx;
    Mask valueMask;
    Bit flag_bit;
    Flag flag_strand;
    Flag flag_end;
    Bit cell_bit;
    Size cell_size;
    Mask headFlag;
    Mask valueMask_dstr;
    
    CordBase():
        bit(20),
        flagEnd((1ULL << 60)),
        mask(0xfffff),
        maskx(0xffffffffff),
        valueMask((1ULL<< 60) - 1),
        flag_bit(61),
        flag_strand(1ULL << flag_bit),
        flag_end(0x1000000000000000),
        cell_bit(4),
        cell_size(16),
        headFlag((1ULL<<63)),
        valueMask_dstr(valueMask | flag_strand)
        {}
    
}_DefaultCordBase;

struct Cord
{
    typedef typename CordBase::CordType CordType;
    typedef typename CordBase::CellType CellType;
    typedef typename CordBase::Flag Flag;
    
    typedef String<CordType> CordString;
    typedef StringSet<CordString> CordSet;
    
    CordType getCordX(CordType const &, typename CordBase::Bit const &, typename CordBase::Mask const &) const;
    CordType getCordY(CordType const &, typename CordBase::Mask const &) const;
    CordType createCord(CordType const &, CordType const &, CordType const &, typename CordBase::Bit const &, typename CordBase::Bit const &) const ;
    CordType hit2Cord(PMRes::HitType const &, typename CordBase::Bit const &, typename CordBase::Mask const &, typename CordBase::Mask const &) const;
    CordType hit2Cord_dstr(PMRes::HitType const &, typename CordBase::Bit const &, typename CordBase::Mask const &, typename CordBase::Mask const &) const;
    CellType cord2Cell(CordType const &, typename CordBase::Bit const &) const;
    CordType cell2Cord(CellType const &, typename CordBase::Bit const &) const;
    void setCordEnd(CordType &, typename CordBase::Flag const &, typename CordBase::Flag const &);
    Flag getCordStrand(CordType const &, CordBase::Bit const &) const;
    Flag isCordEnd(CordType const &, CordBase::Flag const &)const;
    //void setHead(uint64_t &, uint64_t const &, uint64_t const & = _DefaultCordBase.headFlag);
    void setMaxLen(String<uint64_t> &, uint64_t const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t getMaxLen(String<uint64_t> const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & = _DefaultCordBase.bit); //add x and y

    bool isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultCordBase.flagEnd);
}_DefaultCord; 

inline uint64_t get_cord_x (uint64_t val);
inline uint64_t get_cord_y (uint64_t val);
inline uint64_t get_cord_strand (uint64_t val);
inline void cmpRevCord(uint64_t val1, uint64_t val2, uint64_t & cr_val1, uint64_t & cr_val2, uint64_t read_len);
inline uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y);

struct HitBase
{
    uint64_t bit;
    uint64_t bit2;
    uint64_t flag;
    uint64_t flag2;
    uint64_t mask;
    
    HitBase():
        bit(60),
        bit2(61),
        flag(1ULL<<bit),
        flag2(1ULL<<bit2),
        mask(flag - 1)
        {}
}_DefaultHitBase;
/**
 *   struct hit:
 *   extend the structure Cord;
 *   NA[2]|strand[1]|head[1]|genome pos[40]|read pos[20]
 *   NodeType: 1 Head, 0 Body
 */
struct Hit
{
    void setBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockBody(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    bool isBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void unsetBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockStrand(uint64_t &, uint64_t const &, 
                     uint64_t const & = _DefaultHitBase.flag2);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    unsigned getStrand(uint64_t const &, uint64_t const & = _DefaultHitBase.flag2);

}_DefaultHit;

#endif