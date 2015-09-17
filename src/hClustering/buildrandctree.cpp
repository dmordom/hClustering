//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
// This file is also part of OpenWalnut ( http://www.openwalnut.org ).
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
//  buildrandctree
//
//  Build a centroid hierarchical tree from a set of artificially pre-generated set of tractograms yoielding a uniformly random similarity matrix and a seed neighborhood information voxel list.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -r --roi-file:   A text file with the seed voxel coordinates and the corresponding tractogram index (if tractogram naming is based on index rather than coordinates).
//
//   -I --inputf:     Input data folder (containing the compact tractograms).
//
//   -O --outputf:    Output folder where tree files will be written.
//
//  [-d --maxnbdist]: Maximum dissimilarity a seed voxel tract must have to its most similar neighbor not be discarded.
//                     Valid values: (0,1] Use a value of 1 (default) if no discarding is desired.
//
//  [-c --cnbhood]:   Use centroid method with C neighborhood level. Valid values: 6, 18, 24(default), 32, 96, 124.
//
//  [-S --basesize]:  Merge homogeneous base nodes of size S. (mutually exclusive with -N option). Default: 0 (no homogeneous merging).
//
//  [-N --basenum]:   Merge N homogeneous base nodes. (mutually exclusive with -S option). Default: 0 (no homogeneous merging).
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files].
//
//  [-m --cache-mem]: Maximum amount of RAM memory (in GBytes) to use for temporal tractogram cache storing. Valid values [0.1,50]. Default: 0.5.
//
//  [-k --keep-disc]: Keep discarded voxel information in a specialiced section of the tree.
//
//  [--debugout]:     write additional detailed outputs meant to be used for debugging.
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  * Usage example:
//
//   buildrandctree -r roi_lh.txt -I tracograms/ -O results/ -c 26 -N 1000 -k -m 2 -v
//
//  * Outputs (in output folder defined at option -O):
//
//   - 'cX_bin_nmt.txt' - (where X is the neighborhood level defined at option -c) non-monotonic binary-branching hierarchical tree without tree processing (if desired use processtree command).
//   - 'baselist_nmt.txt' - meta-leaves (base nodes defined by the us of option -N or -S) list with IDs corresponding to the non-monotonic tree file.
//   - 'success.txt' - An empty file created when the program has sucessfully exited after completion (to help for automatic re-running scripting after failure).
//   - 'buildrandtree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program.
//
//   [extra outputs when using --debugout option)
//
//   - 'cX_bin_nmt_debug.txt' - version of the counterpart file without '_debug' suffix with redundant information for debugging purposes.
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
#include "randCnbTreeBuilder.h"
#include "common/fileManagerFactory.h"

size_t numComps(0);

bool verbose(false);



int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL) );
        boost::filesystem::path workingDir( boost::filesystem::current_path() );


        // ========== PROGRAM PARAMETERS ==========

        std::string progName( "buildrandctree" );
        std::string configFilename( "/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg" );

        // program parameters
        std::string roiFilename, inputFolder, outputFolder;
        float memory( 0.5 ), maxNbDist( 1 );
        unsigned int nbLevel( 26 ), threads( 0 );
        bool keepDiscarded( false ), niftiMode( true ), debug( false );
        TC_GROWTYPE growType( TC_GROWOFF );
        size_t baseSize( 0 );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions( "Generic options" );
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "roi-file,r", boost::program_options::value< std::string >(&roiFilename), "file with the seed voxels coordinates." )
                ( "inputf,I",  boost::program_options::value< std::string >(&inputFolder), "input data folder (seed tractograms)." )
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder" )
                ( "maxnbdist,d",  boost::program_options::value< float >(&maxNbDist)->implicit_value(1), "[opt] maximum dissimilarity a seed voxel tract must have to its most similar neighbor not be discarded. (0,1]." )
                ( "cnbhood,c",  boost::program_options::value< unsigned int >(&nbLevel)->implicit_value(26), "[opt] centroid method neighborhood level. Valid values: 6, 18, 26(default), 32, 96, 124." )
                ( "basesize,S",  boost::program_options::value< size_t >(&baseSize), "[opt] grow homogeneous base nodes (meta-leaves) of size S. (>=2)." )
                ( "basenum,N",  boost::program_options::value< size_t >(&baseSize), "[opt] grow N homogeneous base nodes (meta-leaves). (>=10)." )
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions( "Configuration" );
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "cache-mem,m",  boost::program_options::value< float >(&memory)->implicit_value(0.5), "[opt] maximum of memory (in GBytes) to use for tractogram cache memory. Default: 0.5." )
                ( "keep-disc,k", "[opt] keep discarded voxels data in a section of the tree file." )
                ( "debugout", "[opt] write additional detailed outputs meant for debug." )
                ( "pthreads,p",  boost::program_options::value< unsigned int >(&threads), "[opt] number of processing cores to run the program in. Default: all available." )
                ;

        // Hidden options, will be allowed both on command line and in config file, but will not be shown to the user.
        boost::program_options::options_description hiddenOptions( "Hidden options" );
        //hiddenOptions.add_options() ;

        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions).add(configOptions).add(hiddenOptions);
        boost::program_options::options_description configFileOptions;
        configFileOptions.add(configOptions).add(hiddenOptions);
        boost::program_options::options_description visibleOptions( "Allowed options" );
        visibleOptions.add(genericOptions).add(configOptions);
        boost::program_options::positional_options_description posOpt; //this arguments do not need to specify the option descriptor when typed in
        //posOpt.add( "roi-file", -1);

        boost::program_options::variables_map variableMap;
        store(boost::program_options::command_line_parser(argc, argv).options(cmdlineOptions).positional(posOpt).run(), variableMap);

        std::ifstream ifs(configFilename.c_str() );
        store(parse_config_file(ifs, configFileOptions), variableMap);
        notify( variableMap);



        if ( variableMap.count( "help" ) )
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
            std::cout << "buildrandctree" << std::endl << std::endl;
            std::cout << "Build a centroid hierarchical tree from a set of artificially pre-generated set of tractograms yoielding a uniformly random similarity matrix and a seed neighborhood information voxel list." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -r --roi-file:   a text file with the seed voxel coordinates and the corresponding tractogram index (if tractogram naming is based on index rather than coordinates)." << std::endl << std::endl;
            std::cout << " -I --inputf:     input data folder (containing the compact tractograms)." << std::endl << std::endl;
            std::cout << " -O --outputf:    output folder where tree files will be written." << std::endl << std::endl;
            std::cout << "[-d --maxnbdist]: maximum dissimilarity a seed voxel tract must have to its most similar neighbor not be discarded." << std::endl;
            std::cout << "                   Valid values: (0,1] Use a value of 1 (default) if no discarding is desired." << std::endl << std::endl;
            std::cout << "[-c --cnbhood]:   use centroid method with C neighborhood level. Valid values: 6, 18, 24(default), 32, 96, 124." << std::endl << std::endl;
            std::cout << "[-S --basesize]:  merge homogeneous base nodes of size S. (mutually exclusive with -N option). Default: 0 (no homogeneous merging)." << std::endl << std::endl;
            std::cout << "[-N --basenum]:   grow N homogeneous base nodes. (mutually exclusive with -S option). Default: 0 (no homogeneous merging)." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	     read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << "[-m --cache-mem]: maximum amount of RAM memory (in GBytes) to use for temporal tractogram cache storing. Valid values [0.1,50]. Default: 0.5." << std::endl << std::endl;
            std::cout << "[-k --keep-disc]: keep discarded voxel information in a specialiced section of the tree." << std::endl << std::endl;
            std::cout << "[--debugout]:     write additional detailed outputs meant to be used for debugging." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " buildrandctree -r roi_lh.txt -I tractograms/ -O results/ -c 26 -N 1000 -k -m 2 -v " << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " - 'cX_bin_nmt.txt' - (where X is the neighborhood level defined at option -c) non-monotonic binary-branching hierarchical tree without tree processing (if desired use processtree command)." << std::endl;
            std::cout << " - 'baselist_nmt.txt' - meta-leaves (base nodes defined by the us of option -N or -S) list with IDs corresponding to the non-monotonic tree file." << std::endl;
            std::cout << " - 'success.txt' - An empty file created when the program has sucessfully exited after completion (to help for automatic re-running scripting after failure)." << std::endl;
            std::cout << " - 'buildrandtree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
            std::cout << std::endl;
            std::cout << " [extra outputs when using --debugout option)" << std::endl << std::endl;
            std::cout << " - 'cX_bin_nmt_debug.txt' - version of the counterpart file without '_debug' suffix with redundant information for debugging purposes." << std::endl;
            std::cout << std::endl;
            exit(0);        }

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


        if( verbose )
        {
            std::cout << "Maximum distance to most similar neighbor: " << maxNbDist << std::endl;
        }

        if ( maxNbDist <= 0 || maxNbDist > 1 )
        {
            std::cerr << "ERROR: distance value used is out of bounds please use a value within (0,1]" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        else if ( maxNbDist == 1 && verbose )
        {
            std::cout << "No neighbor distance restrictions will be applied" << std::endl;
        }
        else if( verbose )
        {
            std::cout << "Seed voxels with no neighbors with tract dissimilarity lower than " << maxNbDist << " will be discarded as outliers" << std::endl;
        }

        if( verbose )
        {
            std::cout << "Centroid neighborhood level: " << nbLevel << std::endl;
        }

        if ( ( nbLevel != 6 ) && ( nbLevel != 18 ) && ( nbLevel != 26 ) && ( nbLevel != 32 ) && ( nbLevel != 92 ) && ( nbLevel != 124 ) )
        {
            std::cerr << "ERROR: invalid nbhood level, only (6,18,26,32,92,124) are accepted" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if ( ( variableMap.count( "basesize" ) && variableMap.count( "basenum" ) ) )
        {
            std::cerr << "ERROR: options --basesize (-S) and --basenum (-N) are mutually exclusive" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if ( variableMap.count( "basesize" ) )
        {
            if( baseSize <= 1 )
            {
                std::cerr << "ERROR: base node (meta-leaf) size must be greater than 1" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            else
            {
                if( verbose )
                {
                    std::cout << "Initial merging stage up to homogeneous base nodes of size: " << baseSize << std::endl;
                }
                growType = TC_GROWSIZE;
            }
        }
        if ( variableMap.count( "basenum" ) )
        {
            if( baseSize < 10 )
            {
                std::cerr << "ERROR: base node (meta-leaf) number must be equal or greater than 10" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            else
            {
                if( verbose )
                {
                    std::cout << "Initial merging stage up to " << baseSize << " homogeneous base nodes (meta-leaves)" << std::endl;
                }
                growType = TC_GROWNUM;
            }
        }

        if( growType == TC_GROWOFF && verbose )
        {
            std::cout << "No homogeneous merging stage" << std::endl;
        }

        if ( variableMap.count( "keep-disc" ) )
        {
            if( verbose )
            {
                std::cout << "Discarded voxel coordinates will be saved in an special section fo the tree file" << std::endl;
            }
            keepDiscarded = true;
        }
        else
        {
            if( verbose )
            {
                std::cout << "Discarded voxel coordinates will not be saved" << std::endl;
            }
            keepDiscarded = false;
        }

        if ( variableMap.count( "inputf" ) )
        {
            if( !boost::filesystem::is_directory( boost::filesystem::path( inputFolder ) ) )
            {
                std::cerr << "ERROR: input seed tractogram folder \"" <<inputFolder<< "\" is not a directory" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            else if( verbose )
            {
                std::cout << "Input seed tractogram folder: " << inputFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no input seed tractogram stated" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if ( variableMap.count( "outputf" ) )
        {
            if( !boost::filesystem::is_directory( boost::filesystem::path( outputFolder ) ) )
            {
                std::cerr << "ERROR: output folder \"" <<outputFolder<< "\" is not a directory" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            else if( verbose )
            {
                std::cout << "Output folder: " << outputFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no output folder stated" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if ( memory < 0.1 || memory > 50)
        {
            std::cerr << "ERROR: cache size must be a positive float between 0.1 and 50 (GB)" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        else if( verbose )
        {
            std::cout << "Tractogram cache memory: " << memory << " GBytes" << std::endl;
        }


        std::string logFilename(outputFolder+"/"+progName+"_log.txt" );
        std::ofstream logFile(logFilename.c_str() );
        if(!logFile)
        {
            std::cerr << "ERROR: unable to open log file: \"" <<logFilename<< "\"" << std::endl;
            exit(-1);
        }


        logFile << "Start Time:\t" << ctime(&programStartTime) << std::endl;
        logFile << "Working directory:\t" << workingDir.string() << std::endl;
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
        logFile << "Vista mode flag:\t" << verbose << std::endl;
        logFile << "Roi file:\t" << roiFilename << std::endl;
        logFile << "Max nb distance:\t" << maxNbDist << std::endl;
        logFile << "Nbhood restriction level:\t" <<nbLevel<< std::endl;
        switch(growType)
        {
        case TC_GROWOFF:
            logFile << "Region growing: None" << std::endl;
            break;
        case TC_GROWSIZE:
            logFile << "Region growing: Size: " << baseSize << std::endl;
            break;
        case TC_GROWNUM:
            logFile << "Region growing: Number: " << baseSize << std::endl;
            break;
        }
        logFile << "Input seed tract folder:\t" << inputFolder << std::endl;
        logFile << "Output folder:\t" << outputFolder << std::endl;
        logFile << "Memory cache size:\t" << memory << " GB" << std::endl;
        logFile << "Debug outputr:\t" << debug << std::endl;
        logFile << "-------------" << std::endl;


        /////////////////////////////////////////////////////////////////


        randCnbTreeBuilder builder( roiFilename, verbose );

        logFile << "Roi size:\t" << builder.roiSize() << std::endl;
        builder.log( &logFile );
        builder.setInputFolder( inputFolder );
        builder.setOutputFolder( outputFolder );
        builder.setDebugOutput( debug );
        builder.buildRandCentroid( nbLevel, memory, growType, baseSize, keepDiscarded );


        /////////////////////////////////////////////////////////////////


        // save and print total time
        time_t programEndTime(time(NULL) );
        int totalTime( difftime(programEndTime,programStartTime) );
        std::cout << "Program Finished, total time: " << totalTime/3600 << "h " <<  (totalTime%3600)/60 << "' " << ((totalTime%3600)%60) << "\"   " << std::endl;
        logFile << "-------------" << std::endl;
        logFile << "Finish Time:\t" << ctime(&programEndTime) << std::endl;
        logFile << "Elapsed time : " << totalTime/3600 << "h " <<  (totalTime%3600)/60 << "' " << ((totalTime%3600)%60) << "\"" << std::endl;


        // create file that indicates process was finished successfully
        std::string successFilename(outputFolder+"/success.txt" );
        std::ofstream successFile(successFilename.c_str() );
        if(!successFile)
        {
            std::cerr << "ERROR: unable to create success file: \"" <<successFile<< "\"" << std::endl;
            exit(-1);
        }
        successFile << "success";


//    }
//    catch(std::exception& e)
//    {
//        std::cout << e.what() << std::endl;
//        return 1;
//    }
    return 0;
}

