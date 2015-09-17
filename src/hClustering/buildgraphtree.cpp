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
//  buildgraphtree
//
//  Build a graph linkage hierarchical tree from a distance matrix built with distBlocks.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -r --roi-file:   A text file with the seed voxel coordinates and the corresponding tractogram index (if tractogram naming is based on index rather than coordinates).
//
//   -g --graph-method: The graph linkage method to recalculate distances, use: 0=single, 1=complete, 2=average, 3=weighted, 4=ward.
//
//   -I --inputf:     Input data folder (containing the distance blocks).
//
//   -O --outputf:    Output folder where tree files will be written.
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) files].
//
//  [--debugout]:     write additional detailed outputs meant to be used for debugging.
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  * Usage example:
//
//   buildgraphtree -r roi_lh.txt -g 2 -I distblocks/ -O results/ -v
//
//
//  * Outputs (in output folder defined at option -O):
//
//   - 'LINKAGE.txt' - (where LINKAGE is a string defining the method chosen in option -g: single/complete/average/weighgted/ward) Contains the output hierarchical tree.
//   - 'buildgraphtree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program.
//
//   [extra outputs when using --debugout option)
//
//   - 'LINKAGE_debug.txt' - tree file with redundant information for debugging purposes.
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
#include "graphTreeBuilder.h"

size_t numComps(0);

bool verbose(false);

// A helper function to simplify the main part.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("buildgraphtree");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, inputFolder, outputFolder;
        unsigned int selector(0), threads(0);
        bool niftiMode( true ), debug( false );
        TG_GRAPHTYPE graphMethod;

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "roi-file,r", boost::program_options::value< std::string >(&roiFilename), "file with the seed voxels coordinates." )
                ( "graph-method,g",  boost::program_options::value< unsigned int >(&selector), "use N graph method (0=single, 1=complete, 2=average, 3=weighted, 4=ward)")
                ( "inputf,I",  boost::program_options::value< std::string >(&inputFolder), "input data folder (distance blocks)." )
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder" )
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "debugout", "[opt] write additional detailed outputs meant for debug." )
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
        //posOpt.add("roi-file", -1);

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
            std::cout << "buildgraphtree" << std::endl << std::endl;
            std::cout << "Build a graph linkage hierarchical tree from a distance matrix built with distBlocks." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -r --roi-file:   a text file with the seed voxel coordinates and the corresponding tractogram index (if tractogram naming is based on index rather than coordinates)." << std::endl << std::endl;
            std::cout << " -g --graph-method: The graph linkage method to recalculate distances, use: 0=single, 1=complete, 2=average, 3=weighted, 4=ward." << std::endl;
            std::cout << " -I --inputf:     input data folder (containing the distance blocks)." << std::endl << std::endl;
            std::cout << " -O --outputf:    output folder where tree files will be written." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	     read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << "[--debugout]:     write additional detailed outputs meant to be used for debugging." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " buildgraphtree -r roi_lh.txt -g 2 -I distblocks/ -O results/ -v" << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " - 'LINKAGE.txt' - (where LINKAGE is a string defining the method chosen in option -g: single/complete/average/weighgted/ward) Contains the output hierarchical tree." << std::endl;
            std::cout << " - 'buildgraphtree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
            std::cout << std::endl;
            std::cout << " [extra outputs when using --debugout option)" << std::endl << std::endl;
            std::cout << " - 'LINKAGE_debug.txt' - tree file with redundant information for debugging purposes." << std::endl;
            std::cout << std::endl;
            exit(0);
        }
        if ( variableMap.count( "verbose" ) ) {
            std::cout << "verbose output" << std::endl;
            verbose=true;
        }

        if ( variableMap.count( "pthreads" ) )
        {
            if ( threads == 1 )
            {
                std::cout << "Using a single processor" << std::endl;
            }
            else if( threads == 0 || threads >= omp_get_num_procs() )
            {
                threads = omp_get_num_procs();
                std::cout << "Using all available processors ( " << threads << " )." << std::endl;
            }
            else
            {
                std::cout << "Using a maximum of " << threads << " processors " << std::endl;
            }
            omp_set_num_threads( threads );
        }
        else
        {
            threads = omp_get_num_procs();
            omp_set_num_threads( threads );
            std::cout << "Using all available processors ( " << threads << " )." << std::endl;
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

        if ( variableMap.count( "debugout" ) )
        {
            if( verbose )
            {
                std::cout << "Debug output files activated" << std::endl;
            }
            debug = true;
        }

        if ( variableMap.count( "version" ) )
        {
            std::cout << progName << ", version 2.0" << std::endl;
            exit(0);
        }

        if ( variableMap.count( "roi-file" ) )
        {
            if( !boost::filesystem::is_regular_file( boost::filesystem::path( roiFilename ) ) )
            {
                std::cerr << "ERROR: roi file \"" <<roiFilename<< "\" is not a regular file" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            else if( verbose )
            {
                std::cout << "Seed voxels roi file: " << roiFilename << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no seed voxels roi file stated" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("inputf")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(inputFolder))) {
                std::cerr << "ERROR: input folder \""<<inputFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "input folder: "<< inputFolder << std::endl;
        } else {
            std::cerr << "ERROR: no input folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }

        if (variableMap.count("outputf")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(outputFolder))) {
                std::cerr << "ERROR: output folder \""<<outputFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "Output folder: "<< outputFolder << std::endl;
        } else {
            std::cerr << "ERROR: no output folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }


        if (variableMap.count("graph-method"))
        {
            if ( (selector<0)||(selector>4))
            {
                std::cerr << "ERROR: invalid graph method"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Graph method. "<< std::flush;
            if (selector==0)
            {
                graphMethod = TG_SINGLE;
                std::cout<<"Single linkage: Ds(k,i+j) = min[D(i,k),D(j,k)]"<<std::endl;
            }
            else if (selector==1)
            {
                graphMethod = TG_COMPLETE;
                std::cout<<"Complete linkage: Dc(k,i+j) = MAX[D(i,k),D(j,k)]"<<std::endl;
            }
            else if (selector==2)
            {
                graphMethod = TG_AVERAGE;
                std::cout<<"Average linkage: Da(k,i+j) = [D(i,k)*Size(i),D(j,k)*size(j)]/[size(i)+size(j)]"<<std::endl;
            }
            else if (selector==3)
            {
                graphMethod = TG_WEIGHTED;
                std::cout<<"Weighted linkage: Dwg(k,i+j) = [D(i,k)+D(i,k)]/2"<<std::endl;
            }
            else
            {
                graphMethod = TG_WARD;
                std::cout<<"Ward linkage: Dwd(k,i+j) = [(Si*Sj)/(Si+Sj)]*[Da(i,k)-Da(i,i)/2-Da(j,j)/2]"<<std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no graph method stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        std::string logFilename(outputFolder+"/"+progName+"_log.txt");
        std::ofstream logFile(logFilename.c_str());
        if(!logFile)
        {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }
        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile << "Verbose:\t" << verbose << std::endl;
        logFile << "Processors used:\t" << threads << std::endl;
        if( niftiMode )
        {
            logFile << "Using nifti file format" << std::endl;
        }
        else
        {
            logFile << "Using vista file format" << std::endl;
        }
        logFile <<"Roi file:\t"<< roiFilename <<std::endl;
        logFile <<"Input folder:\t"<< inputFolder <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;

        if (graphMethod==TG_SINGLE)
            logFile <<"Method used:\tSingle linkage: D(k,i+j) = min[D(i,k),D(j,k)]"<<std::endl;
        else if (graphMethod==TG_COMPLETE)
            logFile <<"Method used:\tComplete linkage: D(k,i+j) = MAX[D(i,k),D(j,k)]"<<std::endl;
        else if (graphMethod==TG_AVERAGE)
            logFile <<"Method used:\tAverage linkage: D(k,i+j) = [D(i,k)*Size(i),D(j,k)*size(j)]/[size(i)+size(j)]"<<std::endl;
        else if (graphMethod==TG_WEIGHTED)
            logFile <<"Method used:\tWeighted linkage: D(k,i+j) = [D(i,k)+D(i,k)]/2"<<std::endl;
        else if (graphMethod==TG_WARD)
            logFile <<"Method used:\tWard linkage: Dwd(k,i+j) = [(Si*Sj)/(Si+Sj)]*[Da(i,k)-Da(i,i)/2-Da(j,j)/2]"<<std::endl;
        else {
            std::cerr << "ERROR: unrecognized graph option"<<std::endl;
            exit(-1);
        }
        logFile << "Debug outputr:\t" << debug << std::endl;
        logFile <<"-------------"<<std::endl;

        /////////////////////////////////////////////////////////////////

        graphTreeBuilder builder(roiFilename, verbose);
        logFile <<"Roi size:\t"<< builder.roiSize() <<std::endl;

        builder.log(&logFile);
        builder.setInputFolder(inputFolder);
        builder.setOutputFolder(outputFolder);
        builder.setDebugOutput( debug );
        builder.buildGraph(graphMethod);

        /////////////////////////////////////////////////////////////////

        // save and print total time
        time_t programEndTime(time(NULL));
        int totalTime( difftime(programEndTime,programStartTime) );
        std::cout <<"Program Finished, total time: "<< totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\"   "<< std::endl;
        logFile <<"-------------"<<std::endl;
        logFile <<"Finish Time:\t"<< ctime(&programEndTime) <<std::endl;
        logFile <<"Elapsed time : "<< totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\""<< std::endl;


//    }
//    catch(std::exception& e)
//    {
//        std::cout << e.what() << std::endl;
//        return 1;
//    }
    return 0;
}

