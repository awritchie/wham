#include "bw_hdf5.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

Interact_H5::Interact_H5(const bw_options &option) : options(option)
{
    options = option;

    std::stringstream ss;
    char outname[1024];
    sprintf(outname,"%s",options.hdf5file.c_str());
    ss << outname;
    ss >> filename;
    const H5std_string FILE_NAME(outname);
    try
    {
        Exception::dontPrint();
        file = new H5File(FILE_NAME, H5F_ACC_RDWR);
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}

void Interact_H5::h5_write_dat(const bw_datfile &filedat)
{
    try
    {
        Exception::dontPrint();
        DataSpace *dataspace;
        DataSet *dataset;
        char name[1024];
        int ncol = filedat.dat[0].size();
        // For each column
        for (int i=0; i<ncol; i++)
        {
            int datsize = filedat.dat.size();
            float values[datsize][2];
            for (int j=0; j<datsize; j++)
            {
                values[j][0] = filedat.dat[j][i];
                values[j][1] = filedat.frameN[j];
            }
            /* Write data */
            sprintf(name,"/Trajectories/Traj-%i/%s",filedat.experiment,options.datnames[i].c_str());
            /* 
             See if it exists already.
             If the dataspace already exists and is the same size as we need,
             just replace the data currently there.  Otherwise, unlink the 
             dataset and make a new one (which will grow the file in size).
            */
            try
            {
                dataset = new DataSet(file->openDataSet(name));
                dataspace = new DataSpace(dataset->getSpace());
                int rank = dataspace->getSimpleExtentNdims();
                hsize_t current_dim[rank];
                dataspace->getSimpleExtentDims(current_dim,NULL);
                if (rank == 2 && current_dim[0] == datsize) {
                    std::cout << "Modifying " << name << " in place\n";
                    dataset->write(values,PredType::NATIVE_FLOAT);
                }
                else {
                    /* Unlink if the dataset already exists */
                    try {
                        file->unlink(name);
                        std::cout << "Unlinking " << name << std::endl;
                        throw std::out_of_range("Unlinking");
                    }
                    catch( FileIException unlink_error ) {
                        std::cout << "\nFailed to unlink " << name << std::endl;
                        std::exit(1);
                    }
                }
            }
            catch(...)
            {
                std::cout << "Writing " << name << std::endl;
                hsize_t dim[2] = { datsize,2 };
                dataspace = new DataSpace(2,dim);
                dataset = new DataSet(file->createDataSet(name,PredType::NATIVE_FLOAT,*dataspace));
                dataset->write(values,PredType::NATIVE_FLOAT);
                /* Attribute */
                StrType str_type(PredType::C_S1,1024);
                hsize_t attrsize = 1;
                char unit_attr[attrsize][1024];
                hsize_t attr_dim[1] = { attrsize };
                sprintf(unit_attr[0],"%s",options.datunits[i].c_str());
                DataSpace attr_dataspace(1,attr_dim);
                Attribute unit_attribute = dataset->createAttribute("Units",str_type,attr_dataspace);
                unit_attribute.write(str_type,unit_attr);
        
                dataset->close();
                delete dataset;
                delete dataspace;
            }
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}

int Interact_H5::h5_get_dataset(const std::string &path, h5_dat &data)
{
    try
    {
        Group grp = file->openGroup(path.c_str());
        // find out which frames this covered
//        hsize_t att_dim[1] = { 2 };
        int attr_data[2];
        Attribute span = grp.openAttribute("First frame read, last frame read");
        span.read(PredType::NATIVE_INT,attr_data);
        data.span.push_back(attr_data[0]);
        data.span.push_back(attr_data[1]);
    }
    catch(...)
    {
        return -1;
    }
    char bin_name[1024],prob_name[1024];
    sprintf(bin_name,"%s/BinCounts",path.c_str());
    sprintf(prob_name,"%s/Probability",path.c_str());
    DataSet *dataset;
    DataSpace *dataspace;
    try
    {
        Exception::dontPrint();
        // read the bin counts
        dataset = new DataSet(file->openDataSet(bin_name));
        dataspace = new DataSpace(dataset->getSpace());
        int rank = dataspace->getSimpleExtentNdims();
        hsize_t dims[rank];
        int ndims = dataspace->getSimpleExtentDims(dims,NULL);
        int bin_out[dims[0]];
        dataset->read(bin_out, PredType::NATIVE_INT);
        delete dataspace;
        delete dataset;
        // read the probabilities
        dataset = new DataSet(file->openDataSet(prob_name));
        dataspace = new DataSpace(dataset->getSpace());
        int prank = dataspace->getSimpleExtentNdims();
        hsize_t pdims[prank];
        int pndims = dataspace->getSimpleExtentDims(pdims,NULL);
        float prob_out[pdims[0]][pdims[1]];
        dataset->read(prob_out, PredType::NATIVE_FLOAT);
        delete dataspace;
        delete dataset; 
        // Make sure the sizes match
        if (pdims[0] != dims[0])
        {
            std::cerr << "\nERROR! " << bin_name << " and " << prob_name << " do not have the same number of bins!\n";
            std::exit(1);
        }
        /* 
           The probability of an individual bin is prob_out. But we 
           divide by the number of times that bin is visited to get 
           the probability of seeing a single frame in that bin.
        */
        for (int i=0;i<dims[0];i++)
        {
            if (bin_out[i] > 0)
            {
                // There are pdims[1] degrees of freedom, so we want the index of ndof - 1
                //data.bin_prob.push_back(prob_out[i][pdims[1]-1]/(float)bin_out[i]);
                /* Don't divide by the number of times the bin is visited yet */
                data.bin_prob.push_back(prob_out[i][pdims[1]-1]);
            }
            else
            {
                data.bin_prob.push_back(0);
            }
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return 0;
}

int Interact_H5::h5_get_dataset(const std::string name, std::vector<std::vector<float> > &full_dat)
{
    DataSet *dataset;
    DataSpace *dataspace;
    try
    {
        Exception::dontPrint();
        // read the bin counts
        dataset = new DataSet(file->openDataSet(name.c_str()));
        dataspace = new DataSpace(dataset->getSpace());
        // how many dimensions is the data set
        int rank = dataspace->getSimpleExtentNdims();
        // make array of the number of members in each dimension
        hsize_t dims[rank];
        int ndims = dataspace->getSimpleExtentDims(dims,NULL);
  
        // let's open the data
        for (int i=0; i<rank; i++) {
            hsize_t col_dims[1];
            hsize_t count[2] = { dims[0], 1 };
            hsize_t offset[2] = { 0, i };
            col_dims[0] = dims[0];
            std::vector<float> column_dat (dims[0],0);
            int rankc = 1;
            DataSpace mspace(rankc,col_dims);
            dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
            dataset->read(&column_dat[0],PredType::NATIVE_FLOAT,mspace,*dataspace);
            full_dat.push_back(column_dat);
        }
        
        // pretty print out
        if (options.bVerbose ) {
            std::cout << name << " dataset rank = " << rank << ", dimensions ";
            for (int i=0; i<rank; i++) {
                if (i==0) { std::cout << dims[i]; }
                else { std::cout << " x " << dims[i]; }
            }
            std::cout << std::endl;
        }
        
        delete dataspace;
        delete dataset;
        return 0;
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        return 1;
        error.printError();
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        return 1;
        error.printError();
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        return 1;
        error.printError();
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        return 1;
        error.printError();
    }
    catch(...)
    {
        std::cerr << "\nUNKNOWN ERROR" << std::endl;
        return 1;
        std::exit(1);
    }
    return 1;
}


void Interact_H5::h5_bin_assignments(std::vector<bw_datfile> &list)
{
    /* Obtain the bin number for each frame */
    char name[1024];
    DataSet *dataset;
    DataSpace *dataspace;
    try
    {
        Exception::dontPrint();
        for (int i=0; i<(int)list.size(); i++)
        {
            sprintf(name,"/Trajectories/Traj-%i/Bin",list[i].experiment);
            //std::cout << name << std::endl;
            dataset = new DataSet(file->openDataSet(name));
            dataspace = new DataSpace(dataset->getSpace());
            int rank = dataspace->getSimpleExtentNdims();
            hsize_t dims[rank];
            int ndims = dataspace->getSimpleExtentDims(dims,NULL);
            int bin_out[dims[0]];
            dataset->read(bin_out,PredType::NATIVE_INT);
            delete dataspace;
            delete dataset;
            for (int n=0; n<dims[0]; n++)
            {
                list[i].bin.push_back(bin_out[n]);
            }
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}

void Interact_H5::h5_write_prob(const std::vector<h5_dat> &probs)
{
    try
    {
        Exception::dontPrint();
        DataSpace *dataspace;
        DataSet *dataset;
        char name[1024];
        // For each column
        int ncol = probs[0].avg.size();
        int nconv = probs.size();
        for (int i=0; i<ncol; i++)
        {
            float averages[nconv][2];
            for (int j=0; j<nconv; j++)
            {
                averages[j][0] = probs[j].avg[i];
                averages[j][1] = probs[j].stdev[i];
            }

            /* Write data */
            sprintf(name,"/Ensemble/%s",options.datnames[i].c_str());
            /*
             See if it exists already.
             If the dataspace already exists and is the same size as we need,
             just replace the data currently there.  Otherwise, unlink the
             dataset and make a new one (which will grow the file in size).             
             */
            
            try
            {
                dataset = new DataSet(file->openDataSet(name));
                dataspace = new DataSpace(dataset->getSpace());
                int rank = dataspace->getSimpleExtentNdims();
                hsize_t current_dim[rank];
                dataspace->getSimpleExtentDims(current_dim,NULL);
                if (rank == 2 && current_dim[0] == nconv) {
                    std::cout << "Modifying " << name << " in place\n";
                    dataset->write(averages,PredType::NATIVE_FLOAT);
                }
                else {
                    /* Unlink if the dataset already exists */
                    try {
                        file->unlink(name);
                        std::cout << "Unlinking " << name << std::endl;
                        throw std::out_of_range("Unlinking");
                    }
                    catch( FileIException unlink_error ) {
                        std::cout << "\nFailed to unlink " << name << std::endl;
                        std::exit(1);
                    }
                }
            }
            catch(...)
            {
                std::cout << "Writing " << name << std::endl;
                hsize_t dim[2] = { nconv, 2 };
                dataspace = new DataSpace(2,dim);
                dataset = new DataSet(file->createDataSet(name,PredType::NATIVE_FLOAT,*dataspace));
                dataset->write(averages,PredType::NATIVE_FLOAT);

                /* Attribute */
                StrType str_type(PredType::C_S1,1024);
                hsize_t attrsize = 2;
                char unit_attr[attrsize][1024];
                hsize_t attr_dim[1] = { attrsize };
                sprintf(unit_attr[0],"average(%s)",options.datunits[i].c_str());
                sprintf(unit_attr[1],"stdev(%s)",options.datunits[i].c_str());
                DataSpace attr_dataspace(1,attr_dim);
                Attribute unit_attribute = dataset->createAttribute("Units",str_type,attr_dataspace);
                unit_attribute.write(str_type,unit_attr);

                dataset->close();
                delete dataset;
                delete dataspace;
            }
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}

void Interact_H5::h5_write_hist(const std::vector<std::vector<double> > &dat, int &index)
{
    try
    {
        Exception::dontPrint();
        DataSpace *dataspace;
        DataSet *dataset;
        Group *group;
        char name[1024];
        
        float adat[(int)dat.size()][2];
        for (int i=0;i<(int)dat.size();i++){
            for (int j=0;j<2;j++){
                adat[i][j] = dat[i][j];
            }
        }
        
        /* If the /Histogram group does not exist, we need to create it */
        try {
            group = new Group(file->openGroup("/Histogram"));
        }
        catch(...) {
            group = new Group(file->createGroup("/Histogram"));
        }
            
        /* Write data */
        sprintf(name,"/Histogram/%s",options.datnames[index].c_str());
        /*
          See if it exists already.
         If the dataspace already exists and is the same size as we need,
         just replace the data currently there.  Otherwise, unlink the
         dataset and make a new one (which will grow the file in size).        */
        try
        {
            dataset = new DataSet(file->openDataSet(name));
            dataspace = new DataSpace(dataset->getSpace());
            int rank = dataspace->getSimpleExtentNdims();
            hsize_t current_dim[rank];
            dataspace->getSimpleExtentDims(current_dim,NULL);
            if (rank == 2 && current_dim[0] == (int)dat.size()) {
                std::cout << "Modifying " << name << " in place\n";
                dataset->write(adat,PredType::NATIVE_FLOAT);
            }
            else {
                /* Unlink if the dataset already exists */
                try {
                    file->unlink(name);
                    std::cout << "Unlinking " << name << std::endl;
                    throw std::out_of_range("Unlinking");
                }
                catch( FileIException unlink_error ) {
                    std::cout << "\nFailed to unlink " << name << std::endl;
                    std::exit(1);
                }
            }
        }
        catch(...)
        {
            std::cout << "Writing " << name << std::endl;
            hsize_t dim[2] = { (int)dat.size(), 2 };
            dataspace = new DataSpace(2,dim);
            dataset = new DataSet(file->createDataSet(name,PredType::NATIVE_FLOAT,*dataspace));
            dataset->write(adat,PredType::NATIVE_FLOAT);
            
            /* Attribute */
            StrType str_type(PredType::C_S1,1024);
            hsize_t attrsize = 2;
            char unit_attr[attrsize][1024];
            hsize_t attr_dim[1] = { attrsize };
            sprintf(unit_attr[0],"Left Edge(%s)",options.datunits[index].c_str());
            sprintf(unit_attr[1],"Probability");
            DataSpace attr_dataspace(1,attr_dim);
            Attribute unit_attribute = dataset->createAttribute("Units",str_type,attr_dataspace);
            unit_attribute.write(str_type,unit_attr);
            dataset->close();
            delete dataset;
            delete dataspace;
        }
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        std::exit(1);
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        std::exit(1);
    }
    return;
}