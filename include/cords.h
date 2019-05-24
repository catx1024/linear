#ifndef LINEAR_HEADER_CORDS_H
#define LINEAR_HEADER_CORDS_H

#include <seqan/sequence.h>

using seqan::String;
using seqan::StringSet;
using seqan::CharString;
using namespace seqan;
/*
 * Cord(C): coordinates in the alignment matrix;
 * :=|N/A[2]|strand[1]|cordEnd[1] gC [40] |rC [20bits]
 * cell[4] the minimum length the window is allowed to slide in the alignment matrix.
 * genomeCord(gC or xC):= position in the genome >> cell_bit << cell_bit. the last ll_bit bits maybe set to 0
 * gC(genome coordinate):= SA node := Seq_id[10] | base_x i2[30]  
 * rC:= base_y[20]. Read coordinate
 */
typedef uint64_t CordType;
typedef String<CordType> CordsType;
typedef StringSet<CordsType > CordsSetType;
struct CordBase
{
    typedef unsigned Bit;
    typedef uint64_t Mask;
    typedef uint64_t Flag;
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
    unsigned bit_id;
    CordBase();
};
extern CordBase _DefaultCordBase;
struct Cord
{
    uint64_t getCordX(uint64_t const &, unsigned const & = _DefaultCordBase.bit, uint64_t const & = _DefaultCordBase.maskx) const;
    uint64_t getCordY(uint64_t const & cord, uint64_t const & mask = _DefaultCordBase.mask) const;
    uint64_t createCord(uint64_t const & x, 
                 uint64_t const & y, 
                 uint64_t const & strand,
                 unsigned const & bit = _DefaultCordBase.bit, 
                 unsigned const & bit2 = _DefaultCordBase.flag_bit) const;
    uint64_t hit2Cord(uint64_t const & hit, 
               unsigned const & bit = _DefaultCordBase.bit, 
               uint64_t const & mask = _DefaultCordBase.mask,
               uint64_t const & mask2 = _DefaultCordBase.valueMask
              ) const;
    uint64_t hit2Cord_dstr(uint64_t const & hit, 
               unsigned const & bit = _DefaultCordBase.bit, 
               uint64_t const & mask = _DefaultCordBase.mask,
               uint64_t const & mask2 = _DefaultCordBase.valueMask_dstr
              ) const;
    uint64_t cord2Cell(uint64_t const & cord, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const;
    uint64_t cell2Cord(uint64_t const & cell, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const;
    void setCordEnd(uint64_t & cord,
            typename CordBase::Flag const & end = _DefaultCordBase.flag_end);
    uint64_t getCordStrand(uint64_t const & cord,
            unsigned const & strand = _DefaultCordBase.flag_bit) const;
    uint64_t isCordEnd(uint64_t const & cord,
                typename CordBase::Flag const & end = _DefaultCordBase.flag_end) const;
    void setMaxLen(String<uint64_t> &, uint64_t const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t getMaxLen(String<uint64_t> const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & = _DefaultCordBase.bit); //add x and y
    bool isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultCordBase.flagEnd);
    uint64_t makeBlockEndVal(uint64_t, uint64_t const & = _DefaultCordBase.flagEnd);
};
extern Cord _DefaultCord; 
/**
 *   struct hit:
 *   extend the structure Cord;
 *   NA[2]|strand[1]|head[1]|genome pos[40]|read pos[20]
 *   NodeType: 1 Head, 0 Body
 */
struct HitBase
{
    uint64_t bit;
    uint64_t bit2;
    uint64_t flag;
    uint64_t flag2;
    uint64_t mask;
    HitBase();
};
extern HitBase _DefaultHitBase;
struct Hit
{
    void setBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockBody(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    bool isBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void unsetBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockStrand(uint64_t &, uint64_t const &, uint64_t const & = _DefaultHitBase.flag2);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    unsigned getStrand(uint64_t const &, uint64_t const & = _DefaultHitBase.flag2);
};
extern Hit _DefaultHit;

void cmpRevCord(uint64_t, uint64_t, uint64_t &, uint64_t &, uint64_t);
uint64_t get_cord_x (uint64_t);
uint64_t get_cord_y (uint64_t); 
uint64_t get_cord_strand (uint64_t);
uint64_t get_cord_id (uint64_t);
uint64_t shift_cord(uint64_t const &, int64_t, int64_t);
uint64_t create_id_x (uint64_t, uint64_t);
uint64_t create_cord (uint64_t, uint64_t, uint64_t, uint64_t);
uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y);
void set_cord_id (uint64_t & val, uint64_t id);
void set_cord_end (uint64_t &); 
void set_cord_block_end(uint64_t & val);
void print_cord(uint64_t, CharString = "");
int64_t atomic_inc_cord_y (int64_t & cord); // atomic cord++, return the new cord

#endif