#include "base.h"
#include "args_parser.h"
#include "mapper.h"

int main(int argc, char ** argv)
{
    double time = sysTime();
    std::cerr << "[]\n";
    (void)argc;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    std::cerr << "Encapsulated version: Mapping reads efficiently" << std::endl;
    Mapper mapper(options);

    //omp_set_num_threads(mapper.thread());
    //map(mapper, options.p1);

    //mapper.print_vcf();
    std::cerr << "  Result Files: \033[1;31m" << options.oPath << "\033[0m" << std::endl;
    //std::cerr << "                \033[1;31m" << (mapper.getOutputPrefix() + ".gvf") << "\033[0m" << std::endl;
    //std::cerr << "                \033[1;31m" << (mapper.getOutputPrefix() + ".sam") << "\033[0m" << std::endl;
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}