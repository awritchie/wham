#include "pca_option_parser.h"

void pca_option_parser(int argc, char *argv[], bw_options &options)
{
    namespace po = boost::program_options;
    
    const char *description =
    {
        "\tSomething about doing PCA.\n\tThe --h5file is generated from bin/wham.\n"
    };
    options.seed = NULL;
    options.doseed = false;
    options.nbootstrap = 1000;
    options.bVerbose = false;
    options.isangle = false;
    options.doWrite = true;
    try
    {
        po::options_description desc(description);
        desc.add_options()
            ("help,h","")
            ("h5file,p",po::value<std::string>(&options.hdf5file)->required(),
                "HDF5 file for the system being examined.  Required")
            ("datatype,d",po::value<std::vector<std::string> >(&options.datnames)->required(),
                "Data type names to read data from HDF5 file.  Do not use spaces.  Required")
            ("nowrite,w",
                "Do not save results to HDF5 file.")
            ("seed,s",po::value<int>(&options.seed),
                "Use a random seed to use for bootstrapping.")
            ("angle,a",
                "Use if the property is an angle and subject to periodic boundary conditions.  ALL PROPERTIES IN THE GIVEN FILE MUST BE ANGLES OR MUST NOT BE ANGLES.  DO NOT MIX ANGLES WITH NON-ANGLES")
            ("nbootstrap,b",po::value<double>(&options.nbootstrap),
                "Number of samples for bootstrap.  Default: 1000")
            ("bVerbose,v",
                "Print extra information")
        ;
        try
        {
            std::cerr << desc;
            
            po::variables_map vm;
            po::store(po::parse_command_line(argc,argv,desc), vm);
            po::notify(vm);
            
            if (vm.count("nowrite"))
            {
                options.doWrite = false;
            }
            if (vm.count("bVerbose"))
            {
                options.bVerbose = true;
            }
            if (vm.count("seed"))
            {
                std::cerr << options.seed << std::endl;
            }
            if (options.nbootstrap > 0)
            {
                if (options.nbootstrap > 1e15)
                {
                    std::cerr << "\n" << std::string(70, '*') << "\n";
                    std::cerr << "Cowardly refusing to perform more than 1e15 bootstrap steps." << std::endl;
                    std::cerr << std::string(70, '*') << "\n";
                    options.nbootstrap = 1e15;
                }
                options.bootstrap = true;
            }
            if (vm.count("angle"))
            {
                options.isangle = true;
            }
            /*if (vm.count("seed"))
            {
                options.doseed = true;
            }*/
            if (vm.count("help"))
            {
                std::exit(1);
            }
        }
        catch(boost::program_options::required_option& e)
        {
            std::cerr << "\nERROR: " << e.what() << std::endl;
            std::exit(1);
        }
    }
    catch(...)
    {
        std::cerr << "\nException of unknown type(s)\n";
        std::exit(1);
    }
    pca_print_options(options);
    return;
}

void pca_print_options(const bw_options &options)
{
    std::cerr << std::string(70, '-') << "\n";
    std::cerr << "Using data contained in " << options.hdf5file << std::endl;
    std::cerr << std::endl;
    std::cerr << "Looking at:\n\t";
    for (int i=0; i<(int)options.datnames.size(); i++) {
        if ( i != 0 ) { std::cerr << ", "; }
        std::cerr << options.datnames[i];
    }
    std::cerr << "\n";
    if (options.isangle)
    {
        std::cerr << "Taking into account angle periodicity..." << std::endl;
    }
    if (! options.doWrite )
    {
        std::cerr << "Will not write to " << options.hdf5file << "." << std::endl;
    }
    if (options.bootstrap)
    {
        std::cerr << "Will bootstrap " << options.nbootstrap << " times." << std::endl;
    }
    if (options.doseed)
    {
        std::cerr << "Using " << options.seed << " as seed" << std::endl;
    }
    std::cerr << std::string(70, '-') << "\n";
    return;
}

