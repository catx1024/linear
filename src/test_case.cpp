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

using ::testing::TestWithParam;
using ::testing::Values;
class FooTest :public ::testing::TestWithParam<char*> {
     public:
  ~FooTest() override {}
  void SetUp() override { val = GetParam(); }
  void TearDown() override {}
 protected:
    char* val;

};

TEST_P(FooTest, DoesBlah) {
  // Inside a test, access the test parameter with the GetParam() method
  // of the TestWithParam<T> class:
ASSERT_EQUAL(1,1);

}
INSTANTIATE_TEST_SUITE_P(DoesBlah,
                         FooTest,
                        ::testing::Values("meeny", "miny", "moe"));

int main(int argc, char ** argv)
{
    /*
    std::cerr << "Linear unit test[]" << std::endl;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
        */
    //Mapper mapper(options);
    testing::InitGoogleTest(&argc, argv);

 //FooTest foo;

    return RUN_ALL_TESTS();
    return 0;
}
