#include <seqan/arg_parse.h>
#include "args_parser.h"
#include "mapper.h"
#include "gap.h"
using namespace seqan; 

struct LinearParm : ParmBase
{
    //MapGapsParm map_gaps_parm;
    //MapGapsParm
    //LinearParm (int thD_tile_size, float thD_err_rate) :
    //    map_gaps_parm(thD_tile_size, thD_err_rate)
    //{}
};

int main(int argc, char const ** argv)
{
    double time = sysTime();
    //LinearParm parm(192, 0.2);

    std::ofstream of;
    //printParmBase2Json (of, parm.map_gaps_parm.map_g_anchor_parm.map_g_anchor2_parm_);
    //(void)argc;
    Options options;
    options.versions = "1.8.2";
    std::cerr << "["<< options.versions
              << "]\nEncapsulated: Mapping reads efficiently" << std::endl;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    Mapper mapper(options);
    map(mapper, options.p1);
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
