#ifndef PCA_OPTION_PARSER_H_
#define PCA_OPTION_PARSER_H_

#include <boost/program_options.hpp>
#include <iterator>
#include <iostream>
#include "bw_types.h"

void pca_option_parser(int argc, char *argv[], bw_options &options);

void pca_print_options(const bw_options &options);

#endif