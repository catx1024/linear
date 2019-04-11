#ifndef LINEAR_HEADER_ARGUMENT_PARSER_H
#define LINEAR_HEADER_ARGUMENT_PARSER_H
#include <seqan/arg_parse.h>
#include "base.h"

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char ** argv);

#endif