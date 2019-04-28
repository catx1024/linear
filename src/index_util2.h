#ifndef SEQAN_HEADER_INDEX_UTIL2_H
#define SEQAN_HEADER_INDEX_UTIL2_H

//using seqan::String;
//using seqan::StringSet;
using namespace seqan;

class DIndex 
{
    //String <int> dir;
    //String <int> hs;
    String <int64_t> dir;
    String <int64_t> hs;
    void * pt_dir[2];
    void * pt_hs[2];
    int pt;
    LShape shape;
public:
    DIndex();
    DIndex(unsigned); //shape_len
    String<int64_t> & getDir();
    String<int64_t> & getHs();
    LShape & getShape();
    int fullSize();
    int getShapeLen();
}; 
int createDIndex(StringSet<String<Dna5> > &, DIndex &, int64_t, int64_t);
int64_t queryHsStr(DIndex & index, int64_t xval);
int64_t queryHsEnd(DIndex & index, int64_t xval);


#endif