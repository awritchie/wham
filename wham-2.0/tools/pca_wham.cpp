#include "pca_wham.h"

PCA::PCA(const bw_options &option) : Interact_H5(option)
{
    options = option;
    // Read in the probability for each bin
    h5_get_dataset("/Ensemble",prob);
    // Read the bin values for each frame.  Also tells us how many trajectories we have
    read_bins();
    // Read all of the datasets selected with --dataset option
    read_datasets();
    // Remove frames missing data (likely due to clustering)
    prune_frames();
    // Count number of times bins are visited and modify (pca_dat) frames appropriately
    count_bins();
    // Make a bootstrapped dataset
    for (int i=0;i<options.nbootstrap;i++) {
        std::vector<pca_dat> sample;
        bootstrap(sample);
        // For each bootstrapped dataset, do PCA
        std::vector<std::vector<double> > dataset;
//        test_dataset(dataset);
        array_dataset(sample,dataset);
        do_pca(dataset);
        if (i > -1 ) {
            break;
        }
    }
    return;
}

void PCA::read_bins()
{
    ntraj = 0;
    int in_traj = 0;
    while (in_traj == 0) {
        std::stringstream dataset;
        dataset << "/Trajectories/Traj-" << ntraj << "/Bin";
        std::vector<std::vector<float> > dat;
        in_traj = h5_get_dataset(dataset.str(),dat);
        if (in_traj == 0) {
            // Rearrange data to be useful here.  In bins, dat[0][x] is binvalue of frame x, dat[1][x] (if it exists, it may not) is frame number of frame x.  However, if dat[1][x] exists, then dat[1][x] = x.
            for (int j=0; j<(int)dat[0].size(); j++) {
                pca_dat framej;
                framej.experiment = ntraj;
                framej.frame = j;
                framej.bin = (int) dat[0][j];
                framej.prob = (double) prob.bin_prob[framej.bin];
                total_frame.push_back(framej);
            }
            ntraj++;
        }
    }
    return;
}

void PCA::read_datasets()
{
    for (int i=0; i<ntraj; i++) {
        for (int j=0; j<(int)options.datnames.size(); j++) {
            std::stringstream dataset;
            dataset << "/Trajectories/Traj-" << i << "/" << options.datnames[j];
            std::vector<std::vector<float> > dat;
            if (h5_get_dataset(dataset.str(),dat) == 1) {
                std::cerr << "\nERROR: Cannot open " << dataset.str() << std::endl;
                std::exit(1);
            }
            // dat[0][x] = value of frame x.  dat[1][x] = frame number of frame x.
            for (int k=0; k<(int)dat[0].size(); k++) {
                int total_frame_number = i*prob.span[1]+dat[1][k];
                total_frame[total_frame_number].dat.push_back(dat[0][k]);
            }
        }
    }
    return;
}

void PCA::prune_frames()
{
    // Start at the end so that the index is not incorrect due to erasing
    // Faster to make a new vector than to delete members.
    /*for (int i=(int)frame.size()-1; i>=0; i--){
        if ( (int)frame[i].dat.size() != (int)options.datnames.size() ) {
            frame.erase(frame.begin() + i);
        }
    }*/
    for (int i=0; i<(int)total_frame.size(); i++){
        if ( (int)total_frame[i].dat.size() == (int)options.datnames.size() ) {
            frame.push_back(total_frame[i]);
        }
    }
    total_frame.clear();
    return;
}

void PCA::count_bins()
{
    int nbins = prob.bin_prob.size();
    std::vector<int> counts (nbins,0);
    for (int i=0; i<(int)frame.size(); i++) {
        counts[frame[i].bin]++;
    }
    // Divide the probability of being in a bin by the number of times that bin was visited
    for (int i=0; i<(int)frame.size(); i++) {
        frame[i].prob /= counts[frame[i].bin];
    }
    probs = std::vector<double> ((int)frame.size(),0);
    for (int i=0; i<(int)frame.size(); i++) {
        probs[i] = frame[i].prob;
    }
    return;
}

void PCA::bootstrap(std::vector<pca_dat> &resample)
{
    // Handle the random number generator seed
    unsigned seed;
    if (options.seed != NULL) {
        seed = options.seed;
        options.seed++;
    }
    else {
        std::random_device rd;
        seed = rd();
    }
    std::mt19937 gen(seed);
    // Sample from a discrete distribution
    std::discrete_distribution<> distribution (probs.begin(),probs.end());
    std::map<int, int> m;
    // Each resample is has nframes
    int nframes = (int)frame.size();
    for (int n=0;n<nframes;n++){
        ++m[distribution(gen)];
    }
    
    // Use the frame selections to generate the resample
    resample = std::vector<pca_dat> (nframes);
    int i = 0;
    for (auto p : m) {
        for (int j=0; j<p.second; j++) {
            resample[i] = frame[p.first];
            i++;
        }
    }
    return;
}

void PCA::test_dataset(std::vector<std::vector<double> > &mat)
{
    mat = std::vector<std::vector<double> > { {-0.693101,-1.437665,0.569708},
                                            {1.959635,1.034623,1.085524},
                                            {0.486191,0.204207,0.991877},
                                            {-1.487539,0.251570,0.189390},
                                            {-1.182513,-0.074457,0.837460},
                                            {0.656417,0.196768,1.682814},
                                            {1.589478,0.137858,0.293443},
                                            {-0.238769,0.084736,0.907351},
                                            {0.765011,1.408654,-0.012985},
                                            {0.022413,-0.839636,-0.047308},
                                            {1.293217,-0.141615,0.585553},
                                            {-0.559727,-0.554367,0.624575},
                                            {-1.378358,-0.448479,-0.848946},
                                            {2.308664,-0.026863,-0.582588},
                                            {-0.613424,1.338016,0.383870},
                                            {0.442568,-0.519120,0.470392},
                                            {-0.267398,-0.412677,0.717908},
                                            {-0.729059,-0.414493,0.313460},
                                            {0.360823,1.163593,-1.426253},
                                            {-0.227367,-0.982915,1.221052},
                                            {1.528242,-0.877785,0.437466},
                                            {0.033554,2.009362,2.643489},
                                            {1.196838,0.411106,1.607066},
                                            {-0.586405,0.049565,0.010722},
                                            {1.765401,1.886084,2.223172},
                                            {1.127229,-0.105636,1.605168},
                                            {0.132893,-1.110334,1.373132},
                                            {2.311513,-1.745904,1.510791},
                                            {1.062413,0.875206,1.425266},
                                            {1.911766,1.268031,1.830525},
                                            {-0.697130,0.971538,-0.708931},
                                            {0.066936,1.920977,-0.459258},
                                            {0.039514,1.669457,2.470940},
                                            {-1.296557,-0.337910,1.960989},
                                            {0.549928,0.822802,1.369806},
                                            {1.178806,1.687805,-0.357288},
                                            {1.855768,-0.135498,1.966774},
                                            {1.614855,0.608582,0.603608},
                                            {2.479533,0.357979,-0.279483},
                                            {1.448398,1.881475,1.393421} };
    int j = 3;
    for (int i=3;i<=(int)options.datnames.size();i++){
        std::cout << options.datnames[j] << std::endl;
        options.datnames.erase(options.datnames.begin() + j);
    }
    return;
}

void PCA::array_dataset(const std::vector<pca_dat> &resample, std::vector<std::vector<double> > &mat)
{
    int nrows = (int)resample.size();
    int ncols = (int)options.datnames.size();
    mat = std::vector<std::vector<double> > (nrows,std::vector<double> (ncols,3));
    std::vector<double> pca_mat (ncols*nrows,0);
    for (int i=0; i<nrows; i++) {
        for (int j=0; j<ncols; j++) {
            mat[i][j] = resample[i].dat[j];
        }
    }
    return;
}
void PCA::scatter(const std::vector<std::vector<double> > &mat, std::vector<double> &scatter_mat)
{
    int nrows = (int)mat.size();
    int ncols = (int)mat[0].size();
    // Compute the mean vector
    std::vector<double> means;
    mean_vector(mat,means);
    // Construct the scatter matrix
    scatter_mat = std::vector<double> (ncols*ncols,0);
    for (int i=0;i<nrows;i++) {
        double difference[ncols];
        for (int j=0;j<ncols;j++){
            difference[j] = mat[i][j] - means[j];
        }
        dgemm('T','N',ncols,ncols,1,1,difference,1,difference,1,1,&scatter_mat[0],ncols);
    }
    if (options.bVerbose) {
        std::cout << "Scatter Matrix:\n";
        for (int i=0;i<ncols;i++){
            std::cout << "\t";
            for (int j=0;j<ncols;j++){
                int n=j*ncols+i;
                std::cout << scatter_mat[n] << " ";
            }
            std::cout << "\n";
        }
    }
    return;
}

void PCA::mean_vector(const std::vector<std::vector<double> > &mat, std::vector<double> &means)
{
    int nframes = (int) mat.size();
    int ndat = (int) mat[0].size();
    means = std::vector<double> (ndat,0.0f);
    std::vector<double> vars (ndat,0.0f);
    for (int i=0;i<nframes;i++) {
        for (int j=0; j<ndat; j++) {
            means[j] += mat[i][j]/nframes;
            vars[j] += mat[i][j]*mat[i][j]/nframes;
        }
    }
    if (options.bVerbose || true ) {
        std::cout << "Average Vector:\n\t";
        for (int i=0;i<ndat;i++){
            std::cout << means[i] << " ";
        }
        std::cout << "\n";
        std::cout << "Standard Deviation Vector:\n\t";
        for (int i=0;i<ndat;i++){
            vars[i] = sqrt(vars[i] - means[i]*means[i]);
            std::cout << vars[i] << " ";
        }
        std::cout << "\n";
    }
    return;
}

void PCA::do_pca(const std::vector<std::vector<double> > &dataset)
{
    int ncols = (int) dataset[0].size();
    // Compute the scatter matrix
    std::vector<double> scatter_mat;
    scatter(dataset,scatter_mat);
    // Compute the eigenvalues and eigenvectors
    // S=eigenvalues, U=eigenvectors
    std::vector<double> U(ncols*ncols,0.0f), S(ncols,0.0f), VT(ncols*ncols,0.0f);
    int info = 0;
    int lwork = MAX(3*ncols+ncols,5*ncols);
    std::vector<double> WORK(lwork,0.0f);
    dgesvd('A','A',ncols,ncols,&scatter_mat[0],ncols,&S[0],&U[0],ncols,&VT[0],ncols,&WORK[0],lwork,info);
    if (info != 0) {
        std::cerr << "\nERROR! Could not compute SVD of scatter matrix.  May require additional linear algebra methods to be programmed!";
        std::exit(1);
        return;
    }
    // This is all standard out stuff
    for (int i=0;i<MIN(3,ncols);i++){
        std::cout << "EIGENVALUE: " << S[i]/S[0] << ",\t[ ";
        double maxv = 0.0f;
        int maxind = 0;
        for (int j=0;j<ncols;j++){
            int n = j + i*ncols;
            std::cout << U[n] << " ";
            if (fabs(U[n]) > fabs(maxv)) {
                maxv = U[n];
                maxind = j;
            }
        }
        std::cout << " ]\n";
        std::cout << i << " : " << maxv << " " << maxind << " " << options.datnames[maxind] << std::endl;
    }
    // Make the N x ncols dimensional transformation matrix
    std::vector<double> W;
    transformation_matrix(U,2,W);
    // Project onto W
    std::vector<std::vector<double> > fulldataset;
    array_dataset(frame,fulldataset);
    project_pca(fulldataset,W);
    return;
}

void PCA::transformation_matrix(const std::vector<double> &evecs, int Ndim, std::vector<double> &W)
{
    int ncols = (int)options.datnames.size();
    W = std::vector<double> (Ndim*ncols);
    for (int i=0;i<ncols;i++){
        for (int j=0;j<Ndim;j++){
            int n = i*Ndim+j;
            int n2= j*ncols+i;
            W[n] = evecs[n2];
        }
    }
    return;
}

void PCA::project_pca(const std::vector<std::vector<double> > &dataset, const std::vector<double> &W)
{
    int nrows = (int)dataset.size();
    int ncols = (int)dataset[0].size();
    std::vector<double> flat_dataset (nrows*ncols,0.0f);
    int Ndim = (int)W.size()/ncols;
    std::vector<double> projection (nrows*Ndim,0.0f);
    for (int i=0;i<nrows;i++){
        for (int j=0;j<ncols;j++){
            int n = j + i * ncols;
            flat_dataset[n] = dataset[i][j];
        }
    }
    std::vector<double> ww ((int)W.size(),0.0f);
    for (int i=0;i<(int)W.size();i++){
        ww[i] = W[i];
    }
    dgemm('N','N',Ndim,nrows,ncols,1,&ww[0],Ndim,&flat_dataset[0],ncols,1,&projection[0],Ndim);
    std::cout << "Results:\n";
    for (int i=0;i<nrows;i++){
        for (int j=0;j<Ndim;j++){
            std::cout << projection[i*Ndim+j] << " ";
        }
        std::cout << "\n";
    }
    return;
}