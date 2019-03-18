#include <seqan/sequence.h>
#include "shape.h" 
using namespace seqan;

class DIndex 
{
public:
   String <uint64_t> dir;
   String <uint64_t> sa;
}
unsigned const dshape_len = 14;
int createDirctIndex(StringSet<String<Dna5> > & seqs, DIndex & index)
{
    int thd_step;
    LShape shape;
    resize (index.dir, (4 << dshape_len), 0);
    for (int i = 0; i < length(seqs); i++)
    {
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