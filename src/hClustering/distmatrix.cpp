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
//  distmatrix
//
//  Compute a pairwise distance matrix between seed voxel compact tracts. Matrix will be divided in sub-blocks for easier & safer computing and storing.
//
//  * Notes:
//         - As matrix will be simmetrical only upper triangle is computed.
//         - Distance metric used is normalized dot product.
//         - Memory and CPU heavy.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -r --roi:        File with the seed voxel coordinates and corresponding tractogram IDs.
//
//   -I --inputf:     Input data folder (containing the seed voxel compact tractograms).
//
//   -O --outputf:    Output folder where distance matrix block files will be written.
//
//  [-t --threshold]: Number of streamlines relative to the total generated that must pass through a tract voxel to be considered for tract similarity
//                     (i.e.: minimum value of a normalized probabilistic tract in natural units to be considered above noise).
//                     Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired.
//
//  [-b --blocksize]: Desired size (in number of elements per row/column) of the blocks the distance matrix will be subdivided in. Choose 0 for maximum size according to available memory. Default: 5000.
//
//  [--start]:        A pair of row-column integers indicating the first block where to start the process. Previous blocks will not be computed.
//
//  [--finish]:       A pair of row-column integers indicating the last block where to finish the process. Posterior blocks will not be computed.
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [-V --vverbose]:  Very verbose output. Writes additional progress information in the standard output.
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files].
//
//  [-m --memory]:    Approximate RAM memory amount to be made available and used by the program (in GBytes). Valid values [0.1,50]. Default: 0.5.
//
//  [-z --zip]:       zip output files.
//
//  [--nolog]:        Treat input tractograms as being normalized in natural units rather than logarithmic.
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  * Usage example:
//
//   distmatrix -r roi_lh.txt -I tracograms/ -O results/ -t 0.001 -b 5000 -v -m 5 -z
//
//
//  * Outputs (in output folder defined at option -O):
//
//   - "roi_index.txt" - A file containing an index matching each seed coordinate to a block number and position within the block.
//   - "dist_block_X_Y.nii(.v)" - Files containing the distance values for the submatrix in pasition XY within the full distance matrix.
//   - "distmatrix_log.txt" - A text log file containing the parameter details and in-run and completion information of the program.
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
#include "distMatComputer.h"



int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("distmatrix");
        std::string configFilename("../../config/"+progName+".cfg");


        // program parameters
        std::string roiFilename;
        std::string tractFolder, outputFolder;
        std::vector<size_t> startBlock, finishBlock;
        std::pair<size_t, size_t> startPair, finishPair;
        size_t blocksize( 5000 );
        unsigned int threads(0);
        bool verbose( false ), veryVerbose( false ), niftiMode( true );
        bool doZip( false ), noLog( false );
        float memory(0.5), relativeThreshold( 0 );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "roi,r",  boost::program_options::value< std::string >(&roiFilename), "file with the seed voxel coordinates and corresponding tractogram IDs")
                ( "inputf,I", boost::program_options::value< std::string >(&tractFolder), "Input data folder (containing the seed voxel compact tractograms).")
                ( "outputf,O",      boost::program_options::value< std::string >(&outputFolder), "output folder where the distance matrix blocks will be written")
                ( "threshold,t", boost::program_options::value< float >(&relativeThreshold)->implicit_value(0), "[opt] noise threshold for the tractograms relative to number of streamlines per tract. [0,1). Default: 0 (no threshold)" )
                ( "blocksize,b", boost::program_options::value< size_t >(&blocksize)->implicit_value(5000), "[opt] size of the blocks in which the matrix will be divided. If 0 maximum size for the available memmory will be used. Default: 5000." )
                ( "start", boost::program_options::value< std::vector<size_t> >(&startBlock)->multitoken(), "[opt] A pair of row-column integers indicating the first block where to start the process. Previous blocks will not be computed." )
                ( "finish", boost::program_options::value< std::vector<size_t> >(&finishBlock)->multitoken(), "[opt] A pair of row-column integers indicating the last block where to finish the process. Posterior blocks will not be computed." )
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vverbose,V", "[opt] very verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "memory,m",  boost::program_options::value< float >(&memory)->implicit_value(0.5), "[opt] maximum of memory (in GBytes) to use for tractogram cache memory. Default: 0.5." )
                ( "zip,z", "[opt] zip output files.")
                ( "nolog", "[opt] treat input tractograms as being in natural units rather than logarithmic" )
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
            std::cout << "distmatrix" << std::endl << std::endl;
            std::cout << "Compute a pairwise distance matrix between seed voxel compact tracts. Matrix will be divided in sub-blocks for easier & safer computing and storing." << std::endl << std::endl;
            std::cout << "* Notes:" << std::endl;
            std::cout << "       - As matrix will be simmetrical only upper triangle is computed." << std::endl;
            std::cout << "       - Distance metric used is normalized dot product." << std::endl;
            std::cout << "       - Memory and CPU heavy." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       Produce extended program help message." << std::endl << std::endl;
            std::cout << " -r --roi:        File with the seed voxel coordinates and corresponding tractogram IDs." << std::endl << std::endl;
            std::cout << " -I --inputf:     Input data folder (containing the seed voxel compact tractograms)." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where distance matrix block files will be written." << std::endl << std::endl;
            std::cout << "[-t --threshold]: Number of streamlines relative to the total generated that must pass through a tract voxel to be considered for tract similarity." << std::endl;
            std::cout << "                    (i.e.: minimum value of a normalized probabilistic tract in natural units to be considered above noise)." << std::endl;
            std::cout << "                    Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired." << std::endl << std::endl;
            std::cout << "[-b --blocksize]: Desired size (in number of elements per row/column) of the blocks the distance matrix will be subdivided in. Choose 0 for maximum size according to available memory. Default: 5000." << std::endl << std::endl;
            std::cout << "[--start]:        A pair of row-column integers indicating the first block where to start the process. Previous blocks will not be computed." << std::endl << std::endl;
            std::cout << "[--finish]:       A pair of row-column integers indicating the last block where to finish the process. Posterior blocks will not be computed." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[-V --vverbose]:  Very verbose output. Writes additional progress information in the standard output." << std::endl << std::endl;
            std::cout << "[--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << "[-m --memory]:    Approximate RAM memory amount to be made available and used by the program (in GBytes). Valid values [0.1,50]. Default: 0.5." << std::endl << std::endl;
            std::cout << "[-z --zip]:       Zip output files." << std::endl << std::endl;
            std::cout << "[--nolog]:        Treat input tractograms as being normalized in natural units rather than logarithmic." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " distmatrix -r roi_lh.txt -I tracograms/ -O results/ -t 0.001 -b 5000 -v -m 5 -z" << std::endl << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " - 'roi_index.txt'' - A file containing an index matching each seed coordinate to a block number and position within the block." << std::endl;
            std::cout << " - 'dist_block_X_Y.nii(.v)'' - Files containing the distance values for the submatrix in pasition XY within the full distance matrix." << std::endl;
            std::cout << " - 'distmatrix_log.txt'' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
            std::cout << std::endl;
            exit(0);
        }


        if ( variableMap.count( "verbose" ) )
        {
            std::cout << "verbose output" << std::endl;
            verbose=true;
        }

        if ( variableMap.count( "vverbose" ) )
        {
            std::cout << "very verbose output" << std::endl;
            veryVerbose=true;
            verbose=true;
        }

        if ( variableMap.count( "version" ) )
        {
            std::cout << progName << ", version 2.0" << std::endl;
            exit(0);
        }

        if (variableMap.count("zip"))
        {
            if( verbose )
            {
                std::cout << "zipping output files"<<std::endl;
            }
            doZip=true;
        }
        else
        {
            doZip=false;
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


        if (variableMap.count("roi"))
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(roiFilename)))
            {
                std::cerr << "ERROR: roi file \""<<roiFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Input roi file: "<< roiFilename << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no input roi file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if (variableMap.count("inputf"))
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path(tractFolder)))
            {
                std::cerr << "ERROR: single tract folder \""<<tractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Single (leaf) tracts folder: "<< tractFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no single tract folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if (variableMap.count("outputf"))
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path(outputFolder)))
            {
                std::cerr << "ERROR: output folder \""<<outputFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Output folder: "<< outputFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no output folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if( verbose )
        {
            std::cout << "Tractogram relative threshold value: " << relativeThreshold << std::endl;
        }
        if ( relativeThreshold < 0 || relativeThreshold >= 1 )
        {
            std::cerr << "ERROR: Threshold value used is out of bounds please use a value within [0,1)" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        else if ( relativeThreshold == 0 && verbose )
        {
            std::cout << "No tractogram thresholding will be applied" << std::endl;
        }
        else if( verbose )
        {
            std::cout << "Tractogram voxels visited by less than " << relativeThreshold * 100 << " % of the streamlines generated will be set to 0 before dissimilarity computation" << std::endl;
        }


        if (memory<0.1 || memory>50) {
            std::cerr << "ERROR: memory size must be a positive float between 0.1 and 50"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if( verbose )
        {
            std::cout <<"Maximum memory available to the program: "<< memory <<" GBytes"<< std::endl;
        }

        if( verbose )
        {
            std::cout <<"Desired distance matrix block size: ";
            if( blocksize == 0 )
            {
                std::cout << "maximum for available memory."<< std::endl;
            }
            else
            {
                std::cout << blocksize << "x" << blocksize << " elements." << std::endl;
            }
        }

        if ( variableMap.count( "start" ) )
        {
            if( startBlock.size() != 2 )
            {
                std::cerr << "ERROR: start block input must be 2 integer indices. " << startBlock.size() << " numbers introduced." << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( startBlock[0] <= startBlock[1] )
            {
                startPair = std::make_pair< size_t, size_t >( startBlock[0], startBlock[1] );
            }
            else
            {
                startPair = std::make_pair< size_t, size_t >( startBlock[1], startBlock[0] );
                std::cerr << "WARNING: As only upper triangle matrix is computed, starting indices will be interpret in the reverse order." << std::endl;
            }
            if( verbose )
            {
                std::cout << "First block to be computed: " << startPair.first << "-" << startPair.second << " (previous blocks will be ignored)." <<  std::endl;
            }
        }
        else
        {
            startPair = std::make_pair< size_t, size_t >( 0, 0 );
        }

        if ( variableMap.count( "finish" ) )
        {
            if( finishBlock.size() != 2 )
            {
                std::cerr << "ERROR: finish block input must be 2 integer indices. " << finishBlock.size() << " numbers introduced." << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( finishBlock[0] <= finishBlock[1] )
            {
                finishPair = std::make_pair< size_t, size_t >( finishBlock[0], finishBlock[1] );
            }
            else
            {
                finishPair = std::make_pair< size_t, size_t >( finishBlock[1], finishBlock[0] );
                std::cerr << "WARNING: As only upper triangle matrix is computed, finishing indices will be interpret in the reverse order." << std::endl;
            }
            if( verbose )
            {
                std::cout << "Last block to be computed: " << finishPair.first << "-" << finishPair.second << " (later blocks will be ignored)." <<  std::endl;
            }
        }
        else
        {
            finishPair = std::make_pair< size_t, size_t >( 0, 0 );
        }


        if ( variableMap.count( "nolog" ) )
        {
            if( verbose )
            {
                std::cout << "Interpreting tracts as having only linear normalization:" <<  std::endl;
            }
            noLog = true;

        }

        ///////////////////////

        std::string logFilename(outputFolder+"/"+progName+"_log.txt" );
        std::ofstream logFile(logFilename.c_str() );
        if(!logFile)
        {
            std::cerr << "ERROR: unable to open log file: \"" <<logFilename<< "\"" << std::endl;
            exit(-1);
        }

        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile << "Verbose:\t" << verbose << std::endl;
        logFile << "Very verbose:\t" << veryVerbose << std::endl;

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
        logFile <<"Tract folder:\t"<< tractFolder <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Relative threshold:\t"<< relativeThreshold <<std::endl;
        logFile <<"Desired block size:\t"<< blocksize <<std::endl;
        logFile <<"Starting block:\t"<< startPair.first << "-" << startPair.second <<std::endl;
        logFile <<"Finishing block:\t";
        if( finishPair.first != 0 && finishPair.second != 0 )
        {
            logFile << finishPair.first << "-" << finishPair.second <<std::endl;
        }
        else
        {
            logFile << "---" << std::endl;
        }
        logFile <<"Tracts read as:\t";
        if( noLog )
        {
            logFile << "Linear" << std::endl;
        }
        else
        {
            logFile << "Logarithmic" << std::endl;
        }
        logFile <<"Zip flag:\t"<< doZip <<std::endl;
        logFile <<"Available memory:\t"<< memory <<" GB"<<std::endl;
        logFile <<"-------------"<<std::endl;


        // ==============================================================================


        distMatComputer distMat( roiFilename, relativeThreshold, verbose, noLog);
        distMat.setInputFolder( tractFolder );
        distMat.setOutputFolder( outputFolder );
        distMat.setBlockSize( memory, blocksize );
        if( startPair.first != 0 && startPair.second != 0 )
        {
            distMat.setStartingBlock( startPair.first, startPair.second );
        }
        if( finishPair.first != 0 && finishPair.second != 0 )
        {
            distMat.setFinishBlock( finishPair.first, finishPair.second );
        }
        if( doZip )
        {
            distMat.storeZipped();
        }
        else
        {
            distMat.storeUnzipped();
        }
        distMat.setVeryVerbose( veryVerbose );

        distMat.doDistBlocks();


        // ==============================================================================

        // save and print total time
        time_t programEndTime(time(NULL));
        int totalTime( difftime(programEndTime,programStartTime) );
        std::cout << "Program Finished, total time: " << totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\""<< std::endl;
        logFile <<"-------------"<<std::endl;
        logFile << "Program Finished, total time: " << totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\""<< std::endl;



//    }
//    catch(std::exception& e)
//    {
//        std::cout << e.what() << std::endl;
//        return 1;
//    }
    return 0;
}


