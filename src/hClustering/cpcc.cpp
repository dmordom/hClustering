//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
//
// For more reference on the underlying algorithm and research they have been used for refer to:
// - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014).
//   A hierarchical method for whole-brain connectivity-based parcellation.
//   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528
// - Moreno-Dominguez, D. (2014).
//   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography.
//   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig.
//   ISBN 978-3-941504-45-5
//
// hClustering is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// http://creativecommons.org/licenses/by-nc/3.0
//
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------
//
//  cpcc
//
//  Compute the cophenetic correlation coefficient (Farris, 1969) of a hierarchical tree.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -t --tree:       File with the hierarchical tree to compute cpcc from.
//
//   -I --inputf:     Input data folder containing the blocks of the precomputed tract pairwise distance matrix.
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files].
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  * Usage example:
//
//   cpcc -t tree_lh.txt -I distBlocks/ -v
//
//
//  * Outputs:
//
//   - Introduces the cpcc value in the #cpcc field of the tree file defined at option -t.
//
//---------------------------------------------------------------------------


// std librabry
#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <stdexcept>
#include <fstream>

// boost library
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// classes
#include "fileManagerFactory.h"
#include "treeManager.h"


int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("cpcc");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string treeFilename;
        std::string distMatrixFolder;
        unsigned int threads(0);
        bool niftiMode( true );
        bool verbose(false);


        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "tree,t",  boost::program_options::value< std::string >(&treeFilename), "file with the hierarchical tree to compute the cpcc value from")
                ( "inputf,I", boost::program_options::value< std::string >(&distMatrixFolder), "Input data folder with the distance matrix files.")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "pthreads,p",  boost::program_options::value< unsigned int >(&threads), "[opt] number of processing cores to run the program in. Default: all available." )
                ;

        // Hidden options, will be allowed both on command line and in config file, but will not be shown to the user.
        boost::program_options::options_description hiddenOptions("Hidden options");
        //hiddenOptions.add_options() ;

        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions).add(configOptions).add(hiddenOptions);
        boost::program_options::options_description configFileOptions;
        configFileOptions.add(configOptions).add(hiddenOptions);
        boost::program_options::options_description visibleOptions("Allowed options");
        visibleOptions.add(genericOptions).add(configOptions);
        boost::program_options::positional_options_description posOpt; //this arguments do not need to specify the option descriptor when typed in
        posOpt.add("tree", -1);

        boost::program_options::variables_map variableMap;
        store(boost::program_options::command_line_parser(argc, argv).options(cmdlineOptions).positional(posOpt).run(), variableMap);

        std::ifstream ifs(configFilename.c_str());
        store(parse_config_file(ifs, configFileOptions), variableMap);
        notify(variableMap);

        if (variableMap.count("help"))
        {
            std::cout << "---------------------------------------------------------------------------" << std::endl;
            std::cout << std::endl;
            std::cout << " Project: hClustering" << std::endl;
            std::cout << std::endl;
            std::cout << " Whole-Brain Connectivity-Based Hierarchical Parcellation Project" << std::endl;
            std::cout << " David Moreno-Dominguez" << std::endl;
            std::cout << " d.mor.dom@gmail.com" << std::endl;
            std::cout << " moreno@cbs.mpg.de" << std::endl;
            std::cout << " www.cbs.mpg.de/~moreno" << std::endl;
            std::cout << std::endl;
            std::cout << " For more reference on the underlying algorithm and research they have been used for refer to:" << std::endl;
            std::cout << " - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014)." << std::endl;
            std::cout << "   A hierarchical method for whole-brain connectivity-based parcellation." << std::endl;
            std::cout << "   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528" << std::endl;
            std::cout << " - Moreno-Dominguez, D. (2014)." << std::endl;
            std::cout << "   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography." << std::endl;
            std::cout << "   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig." << std::endl;
            std::cout << "   ISBN 978-3-941504-45-5" << std::endl;
            std::cout << std::endl;
            std::cout << " hClustering is free software: you can redistribute it and/or modify" << std::endl;
            std::cout << " it under the terms of the GNU Lesser General Public License as published by" << std::endl;
            std::cout << " the Free Software Foundation, either version 3 of the License, or" << std::endl;
            std::cout << " (at your option) any later version." << std::endl;
            std::cout << " http://creativecommons.org/licenses/by-nc/3.0" << std::endl;
            std::cout << std::endl;
            std::cout << " hClustering is distributed in the hope that it will be useful," << std::endl;
            std::cout << " but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
            std::cout << " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
            std::cout << " GNU Lesser General Public License for more details." << std::endl;
            std::cout << std::endl;
            std::cout << "---------------------------------------------------------------------------" << std::endl << std::endl;
            std::cout << "cpcc" << std::endl << std::endl;
            std::cout << "Compute the cophenetic correlation coefficient (Farris, 1969) of a hierarchical tree." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -t --tree:       File with the hierarchical tree to compute cpcc from." << std::endl << std::endl;
            std::cout << " -I --inputf:     Input data folder containing the blocks of the precomputed tract pairwise distance matrix." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	     read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " cpcc -t tree_lh.txt -I distBlocks/ -v" << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs:" << std::endl << std::endl;
            std::cout << " - Introduces the cpcc value in the #cpcc field of the tree file defined at option -t." << std::endl;
            std::cout << std::endl;
            exit(0);
        }
        if (variableMap.count("version")) {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }
        if (variableMap.count("verbose")) {
            std::cout << "verbose output"<<std::endl;
            verbose=true;
        }
        if (variableMap.count("pthreads"))
        {
            if (threads==1)
            {
                std::cout <<"Using a single processor"<< std::endl;
            }
            else if(threads==0 || threads>=omp_get_max_threads())
            {
                threads = omp_get_max_threads();
                std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
            }
            else
            {
                std::cout <<"Using a maximum of "<< threads <<" processors "<< std::endl;
            }
            omp_set_num_threads( threads );
        }
        else
        {
            threads = omp_get_max_threads();
            omp_set_num_threads( threads );
            std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
        }

        if (variableMap.count("tree")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename))) {
                std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Input tree file: "<< treeFilename << std::endl;
        } else {
            std::cerr << "ERROR: no tree file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("inputf")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(distMatrixFolder))) {
                std::cerr << "ERROR: distance matrix folder \""<<distMatrixFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Distance matrix folder: "<< distMatrixFolder << std::endl;
        } else {
            std::cerr << "ERROR: no distance matrix folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if ( variableMap.count( "vista" ) )
        {
            if( verbose )
            {
                std::cout << "Using vista format" << std::endl;
            }
            fileManagerFactory fmf;
            fmf.setVista();
            niftiMode = false;
        }
        else
        {
            if( verbose )
            {
                std::cout << "Using nifti format" << std::endl;
            }
            fileManagerFactory fmf;
            fmf.setNifti();
            niftiMode = true;
        }

        /////////////////////////////////////////////////////////////////

        WHtree tree(treeFilename);
        if (verbose)
        {
            std::cout<<tree.getReport()<<std::endl;
        }
        if( !tree.isLoaded() )
        {
            std::cerr << "Error while loading tree "<< std::endl;
            exit(-1);
        }

        treeManager treeMngr(&tree, verbose);



        treeMngr.setDistMatrixFolder(distMatrixFolder);
        float cpcc(treeMngr.doCpcc());

        tree.writeTree(treeFilename,niftiMode);

        std::cout<<std::endl<<std::endl<<"CPCC: "<< cpcc <<std::endl<<std::endl;

        /////////////////////////////////////////////////////////////////

        // save and print total time
        time_t programEndTime(time(NULL));
        int totalTime( difftime(programEndTime,programStartTime) );
        std::cout << "Program Finished, total time: " << totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\""<< std::endl;



//    }
//    catch(std::exception& e)
//    {
//        std::cout << e.what() << std::endl;
//        return 1;
//    }
    return 0;
}

