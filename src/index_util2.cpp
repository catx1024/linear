#include <seqan/sequence.h>
#include "shape_extend.h" 
using namespace seqan;

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