#include <seqan/sequence.h>
#include "cord.h"
using namespace seqan;

CordBase::CordBase():
        bit(20),
        flagEnd(1ULL << 60),
        mask((1ULL << 20) - 1),
        maskx(0xffffffffff),
        valueMask((1ULL<< 60) - 1),
        flag_bit(61),
        flag_strand(1ULL << flag_bit),
        flag_end(0x1000000000000000),
        cell_bit(4),
        cell_size(16),
        headFlag((1ULL<<63)),
        valueMask_dstr(valueMask | flag_strand),
        bit_id (40)
{}
CordBase _DefaultCordBase;   
Cord _DefaultCord;
HitBase::HitBase():
        bit(60),
        bit2(61),
        flag(1ULL<<bit),
        flag2(1ULL<<bit2),
        mask(flag - 1)
{}
HitBase _DefaultHitBase;
Hit _DefaultHit;

uint64_t Cord::getCordX(uint64_t const & cord, 
               unsigned const & bit,
               uint64_t const & mask) const
{
    return (cord >> bit) & mask; 
}

uint64_t Cord::getCordY(uint64_t const & cord, 
               uint64_t const & mask) const 
{
    return cord & mask;
}

uint64_t Cord::createCord(uint64_t const & x, 
                 uint64_t const & y, 
                 uint64_t const & strand,
                 unsigned const & bit, 
                 unsigned const & bit2) const
{
    return (x << bit) + y + (strand << bit2);
}

uint64_t Cord::hit2Cord(uint64_t const & hit, 
               unsigned const & bit, 
               uint64_t const & mask,
               uint64_t const & mask2
              ) const
{
    return (hit + ((hit & mask) << bit)) & mask2;
}

uint64_t Cord::hit2Cord_dstr(uint64_t const & hit, 
               unsigned const & bit, 
               uint64_t const & mask,
               uint64_t const & mask2
              ) const
{
    return (hit + ((hit & mask) << bit)) & mask2;
}

uint64_t Cord::cord2Cell(uint64_t const & cord, 
                unsigned const & bit) const
{
    return cord >> bit;
}

 uint64_t Cord::cell2Cord(uint64_t const & cell, 
                unsigned const & bit) const
{
    return cell << bit;
}

void Cord::setCordEnd(uint64_t & cord,
            typename CordBase::Flag const & end)
{
    cord |= end;
}

typename CordBase::Flag Cord::getCordStrand(uint64_t const & cord,
            unsigned const & strand) const
{
    return (cord >> strand) & 1ULL;
}

typename CordBase::Flag Cord::isCordEnd(uint64_t const & cord,
                typename CordBase::Flag const & end) const
{
    return cord & end;
}

uint64_t Cord::shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & bit) //add x and y
{
    if (x < 0)
        return val - ((-x) << bit) + y;
    else
        return val + (x << bit) + y;
}

bool Cord::isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd)
{
    int64_t dx = _DefaultCord.getCordX(val2 - val1);
    int64_t dy = get_cord_y(val2 - val1);
    return (dx >= 0) && (dx < thd) && (dy >= 0) && (dy < thd);
}

bool Cord::isBlockEnd(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}


const unsigned _sa1_bit_ = 60;
const unsigned _sa2_bit_ = 40;
const uint64_t _sa1_mask_ = (1ULL << _sa1_bit_ - _sa2_bit_) - 1;
const uint64_t _sa2_mask_ = (1ULL << _sa2_bit_) - 1;
uint64_t _getSA_i1(uint64_t const & node)
{
    return (node >> _sa2_bit_) & _sa1_mask_;
}
uint64_t _getSA_i2(uint64_t const & node)
{
    return node & _sa2_mask_;
}
uint64_t get_cord_x (uint64_t val) {return _getSA_i2(_DefaultCord.getCordX(val));}
uint64_t get_cord_y (uint64_t val) {return _DefaultCord.getCordY(val);}
uint64_t get_cord_strand (uint64_t val) {return _DefaultCord.getCordStrand(val);}
uint64_t get_cord_id (uint64_t val) {return _getSA_i1(_DefaultCord.getCordX(val));}
void set_cord_y (uint64_t & cord, uint64_t y)
{
    cord &= ~(_DefaultCordBase.mask);
    cord += y;
}
void set_cord_end (uint64_t & val) {_DefaultCord.setCordEnd(val);}
uint64_t create_id_x(uint64_t const id, uint64_t const x)
{
    return (id << _DefaultCordBase.bit_id) + x;
}
uint64_t create_cord (uint64_t id, uint64_t cordx, uint64_t cordy, uint64_t strand)
{
    return _DefaultCord.createCord(create_id_x (id, cordx), cordy, strand);
}

void cmpRevCord(uint64_t val1, 
                    uint64_t val2,
                    uint64_t & cr_val1,
                    uint64_t & cr_val2,
                    uint64_t read_len)
{
    cr_val1 = (val1 - get_cord_y(val1) + read_len - get_cord_y(val2) - 1) ^ _DefaultCordBase.flag_strand;
    cr_val2 = (val2 - get_cord_y(val1) + read_len - get_cord_y(val2) - 1) ^ _DefaultCordBase.flag_strand;
}
uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y)
{
    return (val & (~_DefaultCordBase.valueMask)) + (x << _DefaultCordBase.bit) + y;
}

 void Hit::setBlockStart(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}

 void Hit::setBlockBody(uint64_t & val, uint64_t const & flag)
{
    val &= (~flag);
}

 bool Hit::isBlockStart(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}

 void Hit::setBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}

 void Hit::unsetBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val &= ~flag;
}

 void Hit::setBlockStrand(uint64_t & val, uint64_t const & strand, uint64_t const & flag)
{
    if (strand)
        val |= flag;
    else
        val &= ~flag;
}

 bool Hit::isBlockEnd(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}

 unsigned Hit::getStrand(uint64_t const & val, uint64_t const & flag)
{
    return (val & flag)?1:0;
}

void _printHit(unsigned j, unsigned id1, unsigned id2, String<uint64_t> & hit, unsigned len)
{
    unsigned end;
    for (unsigned k = 0; k < length(hit); k++)
    {
        if (_DefaultHit.isBlockEnd(hit[k]))
            end = 1;
        else
            end = 0;
        //printf("[printhit] %d %d %d %d %d\n", j, id1, id2, len, end);
    }
}

void _printHit(String<uint64_t>  & hit)
{
    for (unsigned k = 0; k < length(hit); k++)
    {
        std::cout << "[P]::_printHit() " 
              << _getSA_i1(_DefaultCord.getCordX(_DefaultCord.hit2Cord(hit[k]))) << " " 
              << _getSA_i2(_DefaultCord.getCordX(_DefaultCord.hit2Cord(hit[k]))) << " " 
              << get_cord_y(hit[k]) << "\n";
        if (_DefaultHit.isBlockEnd(hit[k]))
        {
            std::cout << "[P]::_printHit() end\n";
        }
    }
}
