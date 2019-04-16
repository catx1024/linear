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

int testCreateDIndx(DIndex & index_d, Mapper & mapper)
{
    createDIndex (mapper.genomes(), index_d);
}

TEST(DIndex, create) 
{
    //createDIndex()
    ASSERT_EQ(1,1) ;
}

int main(int argc, char ** argv)
{
    std::cerr << "Linear unit test[]" << std::endl;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper mapper(options);
    DIndex index_d(14);
    testCreateDIndx(index_d, mapper);
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
    return 0;
}
