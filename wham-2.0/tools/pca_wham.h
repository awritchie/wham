#ifndef PCA_WHAM_H_
#define PCA_WHAM_H_

#include "bw_types.h"
#include "pca_types.h"
#include "pca_option_parser.h"
#include "bw_hdf5.h"
#include "cpp_lapack.h"
#include "cpp_blas.h"
#include <random>

class PCA : public Interact_H5
{
protected :
    bw_options options;
    h5_dat prob;
    std::vector<pca_dat> frame, total_frame;
    std::vector<double> probs;
    int ntraj;
public :
    // constructor
    PCA(const bw_options &option);
    // 1) Read in all of the frames, bins, and probabilities
    void read_bins();
    void read_datasets();
    void prune_frames();
    void count_bins();
    // 2) Bootstrap sample frames
    void bootstrap(std::vector<pca_dat> &resample);
    // 4) Usual PCA routine
    void test_dataset(std::vector<std::vector<double> > &mat);
    void array_dataset(const std::vector<pca_dat> &resample, std::vector<std::vector<double> > &mat);
    void mean_vector(const std::vector<std::vector<double> > &mat, std::vector<double> &means);
    void scatter(const std::vector<std::vector<double> > &mat, std::vector<double> &scatter_mat);
    void transformation_matrix(const std::vector<double> &evecs, int Ndim, std::vector<double> &W);
    void project_pca(const std::vector<std::vector<double> > &dataset, const std::vector<double> &W);
    void do_pca(const std::vector<std::vector<double> > &dataset);
    // deconstructor
    ~PCA(){};
};

#endif