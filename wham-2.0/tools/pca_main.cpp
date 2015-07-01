#include "pca_wham.h"

int main(int argc, char * argv[]) {
    bw_options opt;
    pca_option_parser(argc, argv, opt);
    
    PCA results(opt);
    return 0;
}