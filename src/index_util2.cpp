#include <seqan/sequence.h>
#include "shape_extend.h" 
#include "index_util.h"
#include "index_util2.h"
#include "cord.h"
using namespace seqan;

unsigned dshape_len = 22;
DIndex::DIndex():
    shape(dshape_len)
{}
DIndex::DIndex (unsigned len):
    shape(len)
{}
LShape & DIndex::getShape()
{
    return shape;
}
int DIndex::fullSize()
{
    return 1 << shape.weight << shape.weight;
}
String<int64_t> & DIndex::getDir()
{
    return dir;
}
String<int64_t> & DIndex::getHs()
{
    return hs;
}
void DIndex::setHsValue (int64_t i_hs, int64_t gId, int64_t gPos, LShape & shape)
{
    hs[i_hs] = create_cord(gId, gPos, 0, shape.strand);
}

int createDIndex(StringSet<String<Dna5> > & seqs, DIndex & index)
{
    int thd_step;
    LShape & shape = index.getShape();
    String<int64_t> &  dir = index.getDir();
    String<int64_t> & hs = index.getHs();
    resize (index.getDir(), index.fullSize(), 0);
    for (int64_t i = 0; i < length(seqs); i++)
    {
        hashInit (shape, begin(seqs[i]));
        for (int j = 0; j < length(seqs[i]); j++)
        {
            hashNext (shape, begin(seqs[i]) + j);
            ++dir[shape.XValue];
        }
    }
    int64_t sum = 0; 
    for (int64_t i = 0; i < length(dir); i++)
    {
        sum += dir[i];
        dir[i] += sum - dir[i];
    }
    resize (index.getHs(), lengthSum(seqs));
    for (int64_t i = 0; i < length(seqs); i++)
    {
        hashInit (shape, begin(seqs[i]));
        for (int64_t j = 0; j < length(seqs[i]); j++)
        {
            hashNext (shape, begin(seqs[i]) + j);
            index.setHsValue(dir[shape.XValue], i, j, shape);
        }  
    }
}