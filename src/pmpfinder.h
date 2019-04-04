#ifndef SEQAN_HEADER_PMP_FINDER_H
#define SEQAN_HEADER_PMP_FINDER_H
#include <seqan/sequence.h>
using namespace seqan;

//WARNING:The length of read should be < 1MB;
extern const float band_width;
extern const unsigned cmask;
extern const unsigned cell_size;
extern const unsigned cell_num;
extern const unsigned window_width; //16*12
extern const unsigned window_delta;
extern const unsigned sup;
extern const unsigned med;
extern const unsigned inf;
extern const unsigned initx; 
extern const unsigned inity;
extern const unsigned scriptStep;
extern const unsigned scriptBit;
extern const unsigned scriptWindow; //script_length = 2^scriptWindow
extern const unsigned scriptWindow2;
extern const unsigned scriptWindow3;
extern const int scriptCount[5];
extern const int scriptMask;
extern const int scriptMask2;
extern const int scriptMask3;
extern const uint64_t hmask;
extern const unsigned windowThreshold; // 36;

typedef Iterator <String <Dna5> >::Type TIter5;
  
void cmpRevCord(uint64_t, uint64_t, uint64_t &, uint64_t &, uint64_t);
uint64_t get_cord_x (uint64_t);
uint64_t get_cord_y (uint64_t); 
uint64_t get_cord_strand (uint64_t);
uint64_t get_cord_id (uint64_t);
uint64_t shift_cord(uint64_t const &, int64_t, int64_t);
uint64_t create_id_x (uint64_t, uint64_t);
uint64_t create_cord (uint64_t, uint64_t, uint64_t, uint64_t);
uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y);
void set_cord_end (uint64_t &); 

void setCordsMaxLen(String<uint64_t> &, uint64_t);
uint64_t getCordsMaxLen(String<uint64_t> const &);

void createFeatures(TIter5 const &, TIter5 const &, String<short> & );
void createFeatures(StringSet<String<Dna5> > &, StringSet<String<short> > &, unsigned);
unsigned _windowDist(Iterator<String<short> >::Type const &, 
                     Iterator<String<short> >::Type const &);

bool path_dst(typename Iterator<String<uint64_t> >::Type, 
              typename Iterator<String<uint64_t> >::Type, 
              StringSet<String<short> > &,
              StringSet<String<short> > &, 
              String<uint64_t> &,
              float const & );

int extendPatch(StringSet<String<short> > & f1, 
                StringSet<String<short> > & f2, 
                String<uint64_t> & cords,
                int k,
                uint64_t cord1,
                uint64_t cord2,
                int revscomp_const,
                int overlap_size = window_width,
                int gap_size = window_width);
#endif