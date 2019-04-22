#include <gtest/gtest.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "mapper.h"
#include "pmpfinder.h"
#include "cord.h"
#include "chain_map.h"
#include "gap.h"
#include "align_interface.h"
#include "args_parser.h"
#include "index_util2.h"

using ::testing::TestWithParam;
using ::testing::Values;

int testCreateDIndx(DIndex & index, StringSet<String<Dna5> > & seqs)
{
    createDIndex (seqs, index);
}
int testQueryDIndx (DIndex & index, String<Dna5> & seq)
{
    LShape & shape = index.getShape();
    hashInit(shape, begin(seq));
    int64_t pre_xval = ~0;
    bool find_f = false;
    int count_f = 0;
    for (int i = 0; i < length(seq); i++) 
    {
        hashNext(shape, begin(seq) + i);
        int64_t str = queryHsStr(index, shape.XValue);
        int64_t end = queryHsEnd(index, shape.XValue);
        for (int64_t j = str; j < end; j++)
        {
            int64_t val = index.getHs()[j];
            if (get_cord_x(val) == i)
            {
                //std::cout << i << " " << shape.XValue << "\n";
                find_f = true;
                break;
            }
        }
        if (!find_f)
        {
            //std::cerr << "err " << i << " " << shape.XValue << "\n";
            ++count_f;
        }
        find_f = false;
    }
    std::cout << "count_f " << count_f << "\n";
}
//TEST(DIndex, create) 
//{
//    createDIndex()
//    ASSERT_EQ(1,1) ;
//}

int main(int argc, char ** argv)
{
    std::cerr << "Linear unit test[]" << std::endl;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper mapper(options);
    DIndex index_d(22);
    testCreateDIndx(index_d, mapper.genomes());
    //testQueryDIndx(index_d, mapper.genomes()[0]);
    std::cout << "Lshape.len " << index_d.getShape().weight << "\n";


    //testing::InitGoogleTest(&argc, argv);
    //return RUN_ALL_TESTS();
    return 0;
}
