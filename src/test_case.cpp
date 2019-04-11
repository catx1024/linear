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

//using namespace seqan; 

TEST(SquareRootTest, PositiveNos) { 
    ASSERT_EQ(6, 30);
}

int main(int argc, char ** argv)
{
    std::cerr << "Encapsulated version: Mapping reads efficiently" << std::endl;
    //(void)argc;
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper mapper(options);
    //double t=sysTime();

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

 
    //g_test(mapper);

    //std::cerr << "results saved to " << options.getOutputPath() << "\n";
    
    return 0;
}
