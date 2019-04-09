#include <seqan/sequence.h>
#include "shape.h" 
using namespace seqan;

class DIndex 
{
    String <int> dir1
    String <int> hs1
    String <uint64_t> dir2;
    String <uint64_t> hs2;
    void * pt_dir[2];
    void * pt_hs[2];
    int pt;
    int shape_len;
public:
    DIndex();
    getDir();
    getHs();
    int fullSize();
}
unsigned const dshape_len = 14;
int createDirctIndex(StringSet<String<Dna5> > & seqs, DIndex & index)
{
    int thd_step;
    LShape shape;
    resize (index.getDir(), index.fullSize());
    for (int i = 0; i < length(seqs); i++)
    {
        index.getDirPt()[index.getPt()] ->
        hashInit (shape, begin(seqs[i]));
        for (int j = 0; j < length(seqs); j++)
        {
            hashNext (shape, begin(seqs[i]) + j);
            ++dir[shape.hValue];
        }
    }
    for (int i = 0; i < length(index.dir); i++)
    {
        dir[i] += dir[i - 1];
    }
    for (int i = 0; i < length(seqs); i++)
    {
        hashInit (shape, begin(seqs[i]));
        for (int j = 0; j < length(seqs); j++)
        {
            hashNext (shape, begin(seqs[i]) + j);
            sa[dir[shape.hValue]] = createCord();
        }  
    }
}