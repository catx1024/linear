#ifndef SEQAN_HEADER_INDEX_UTIL2_H
#define SEQAN_HEADER_INDEX_UTIL2_H

class DIndex 
{
    String <int> dir1;
    String <int> hs1;
    String <uint64_t> dir2;
    String <uint64_t> hs2;
    void * pt_dir[2];
    void * pt_hs[2];
    int pt;
    int shape_len;
public:
    DIndex();
    //getDir();
    //getHs();
    int fullSize();
}; 
int createDirctIndex(StringSet<String<Dna5> > & seqs, DIndex & index);


#endif