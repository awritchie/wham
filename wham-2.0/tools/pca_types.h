#ifndef PCA_TYPES_H_
#define PCA_TYPES_H_

#include <vector>
#include <string>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <fstream>
#include <iostream>

#ifndef M_PI
#define M_PI = atan(1.)*4.
#endif

#ifndef RAD2DEG
#define RAD2DEG 180./M_PI
#endif

#ifndef DEG2RAD
#define DEG2RAD M_PI/180.
#endif

struct pca_dat
{
    int experiment;
    int frame;
    int bin;
    double prob;
    std::vector<double> dat;
};


#endif
