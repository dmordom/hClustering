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
//  fliptracts
//
//  Flip the desired precomputed compact tractograms in the X-axis. Requires the white matter mask used to compact them.
//   Options are: flip all leaf tracts defined in a roi file or flip the desired node tracts by ID or as defiend in a tree file.
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -m --mask:       White matter mask image that was used to compact the tracts.
//
//  [-r --roi]:      [xor with -t and -n] Roi file with leaf coordinates/trackIDs of tractograms to filp .
//
//  [-t --treebases]:[xor with -r] Hierarchical tree with base-node ids of precomputed tractograms to flip.
//
//  [-a --all]:      [use only alongside -t] When used, all precomputed tree-node corresponding tractograms will be flipped.
//
//  [-n --nodes]:    [xor with -r] A sequence of ids of precomputed node-tractograms to flip.
//
//   -I --inputf:     Folder with the input tractograms to be flipped. These should be leaf tracts for option -r and node tracts for options -t and -n.
//
//   -O --outputf:    Output folder where results be written.
//
//  [-o --no-ow]:     When used tracts will NOT be overwritten if the corresponding output file already exist.
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files.
//
//  [-z --zip]:       zip output files.
//
//  [-F --ufloat]:    use float32 representation to write output tracts (default is uint8).
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  example:
//
//  fliptracts -m wmask.nii -t tree.txt -n 20 234 -I meantracts/ -O results/ -v
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
#include "WHtree.h"
#include "roiLoader.h"

// A helper function to simplify the main part.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " " ) );
    return os;
}

int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("fliptracts");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, inputFolder, outputFolder, maskFilename, treeFilename;
        unsigned int threads(0);
        std::vector<size_t> inNodes;
        bool roiMode( false ), treeBaseMode ( false ), treeAllMode( false ), nodesMode(false), noOverwrite( false );
        bool verbose( false ), niftiMode( true ), doZip( false ), useFloat( false );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions( "Generic options" );
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "mask,m",  boost::program_options::value< std::string >( &maskFilename ), "White matter mask that was used to compact the tracts." )
                ( "roi,r",  boost::program_options::value< std::string >( &roiFilename ), "[opt | xor with -t and -n] Roi file. Will flip all leaf tracts from the roi file." )
                ( "treebases,t",  boost::program_options::value< std::string >( &treeFilename ), "[opt | xor with -r] Tree file. Will flip all precomputed tree base-nodes (meta-leaves) tracts." )
                ( "all,a", "[opt | use only with -t] Flip all precomputed tree nodes (not just base-nodes).")
                ( "nodes,n", boost::program_options::value< std::vector<size_t> >( &inNodes )->multitoken(), "[opt | xor with -r] Node tract IDs corresponding to mean tracts to flip." )
                ( "inputf,I",  boost::program_options::value< std::string >( &inputFolder ), "folder with the input tractograms. Should be leaf tracts for option -r and node tracts for options -t and -n." )
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder" )
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "no-ow,o", "[opt] Tracts will NOT be overwritten if they already exist.")
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "zip,z", "[opt] zip output files.")
                ( "ufloat,F", "[opt] use float32 representation to write tracts (default is uint8)")
                ( "pthreads,p",  boost::program_options::value< unsigned int >(&threads), "[opt] number of processing cores to run the program in. Default: all available." )
                ;

        // Hidden options, will be allowed both on command line and in config file, but will not be shown to the user.
        boost::program_options::options_description hiddenOptions("Hidden options");
        hiddenOptions.add_options()
                ;


        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions).add(configOptions).add(hiddenOptions);
        boost::program_options::options_description configFileOptions;
        configFileOptions.add(configOptions).add(hiddenOptions);
        boost::program_options::options_description visibleOptions("Allowed options");
        visibleOptions.add(genericOptions).add(configOptions);
        boost::program_options::positional_options_description posOpt; //this arguments do not need to specify the option descriptor when typed in
//        posOpt.add("mean-tracts", -1);



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
            std::cout << "fliptracts" << std::endl << std::endl;
            std::cout << "Compute the mean tractograms from a hierarchical tree nodes and the original leaf tracts and write them in compact form." << std::endl;
            std::cout << " Options are: flip all leaf tracts defined in a roi file or flip the desired node tracts by ID or as defiend in a tree file." << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -m --mask:       White matter mask image that was used to compact the tracts." << std::endl << std::endl;
            std::cout << "[-r --roi]:      [xor with -t and -n] Roi file with leaf coordinates/trackIDs of tractograms to filp ." << std::endl << std::endl;
            std::cout << "[-t --treebases]:[xor with -r] Hierarchical tree with base-node ids of precomputed tractograms to flip." << std::endl << std::endl;
            std::cout << "[-a --all]:      [use only alongside -t] When used, all precomputed tree-node corresponding tractograms will be flipped." << std::endl << std::endl;
            std::cout << "[-n --nodes]:    [xor with -r] A sequence of ids of precomputed node-tractograms to flip." << std::endl << std::endl;
            std::cout << " -I --inputf:     Folder with the input tractograms to be flipped. These should be leaf tracts for option -r and node tracts for options -t and -n." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where results be written." << std::endl << std::endl;
            std::cout << "[-o --no-ow]:     When used tracts will NOT be overwritten if the corresponding output file already exist." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   Verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files." << std::endl << std::endl;
            std::cout << "[-z --zip]:       Zip output files." << std::endl << std::endl;
            std::cout << "[-F --ufloat]:    Use float32 representation to write output tracts (default is uint8)." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "example:" << std::endl << std::endl;
            std::cout << "fliptracts -m wmask.nii -t tree.txt -n 20 234 -I meantracts/ -O results/ -v" << std::endl << std::endl;
            exit(0);
        }

        if ( variableMap.count( "version" ) )
        {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }
        if (variableMap.count( "verbose" ) )
        {
            std::cout << "verbose output"<<std::endl;
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

        if ( variableMap.count( "mask" ) )
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path( maskFilename ) ) )
            {
                std::cerr << "ERROR: mask file \""<<maskFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
            std::cout << "Mask file: "<< maskFilename << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no mask file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if( variableMap.count( "roi" ) && ( variableMap.count( "treebases" ) || variableMap.count( "nodes" ) ) )
        {
            std::cerr << "ERROR: option -r is mutually exclusive with -t and -n" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if( variableMap.count( "roi" ) )
        {
            if( !boost::filesystem::is_regular_file( boost::filesystem::path( roiFilename ) ) )
            {
                std::cerr << "ERROR: roi file \"" << roiFilename << "\" is not a regular file" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Flipping leaf tracts. Roi file: "<< roiFilename << std::endl;
            }
            roiMode = true;
        }

        if( variableMap.count( "treebases" ) )
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename))) {
                std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            treeBaseMode = true;
            if( variableMap.count( "all" ) )
            {
                if( verbose )
                {
                    std::cout << "Flipping all tree tratcts. Tree file: "<< treeFilename << std::endl;
                }
                treeAllMode = true;
            }
            else if( verbose )
            {
                std::cout << "Flipping base node tracts. Tree file: "<< treeFilename << std::endl;
            }
        }

        if( variableMap.count( "nodes" ) )
        {
            nodesMode = true;
            if( verbose )
            {
                std::cout << "Flipping tracts from node ids: "<< inNodes << std::endl;
            }
        }


        if ( variableMap.count( "inputf" ) )
        {
            if( !boost::filesystem::is_directory( boost::filesystem::path( inputFolder ) ) )
            {
                std::cerr << "ERROR: input tractogram folder \"" <<inputFolder<< "\" is not a directory" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            else if( verbose )
            {
                std::cout << "Input tractogram folder: " << inputFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no input tractogram stated" << std::endl;
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


        if (variableMap.count("ufloat"))
        {
            if( verbose )
            {
                std::cout << "writing in float"<<std::endl;
            }
            useFloat=true;
        }
        else
        {
            if( verbose )
            {
                std::cout << "writing in char"<<std::endl;
            }
            useFloat=false;
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

        if (variableMap.count( "no-ow") )
        {
            std::cout << "overwrite disabled"<<std::endl;
            noOverwrite=true;
        }

        std::string logFilename(outputFolder+"/"+progName+"_log.txt");
        std::ofstream logFile(logFilename.c_str());
        if(!logFile) {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }
        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Processors used:\t"<< threads <<std::endl;
        logFile <<"Mask file:\t"<< maskFilename <<std::endl;
        if( niftiMode )
        {
            logFile << "Using nifti file format" << std::endl;
        }
        else
        {
            logFile << "Using vista file format" << std::endl;
        }
        if( useFloat )
        {
            logFile << "Writing in float32" << std::endl;
        }
        else
        {
            logFile << "Writing in uint8" << std::endl;
        }
        if( doZip )
        {
            logFile << "Zipping output tracts" << std::endl;
        }
        logFile <<"No overwrite flag:\t"<< noOverwrite <<std::endl;
        logFile <<"Roi file:\t"<< roiFilename <<std::endl;
        logFile <<"Tree file:\t"<< treeFilename <<std::endl;
        logFile <<"All tree nodes flag:\t"<< treeAllMode <<std::endl;
        logFile <<"Node ids:\t"<< inNodes <<std::endl;
        logFile <<"Tracts folder:\t"<< inputFolder <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"-------------"<<std::endl;

        /////////////////////////////////////////////////////////////////


        time_t lastTime(time(NULL)), startTime(time(NULL));
        size_t tractProg(0), realDone(0);

        fileManagerFactory tractReaderFileMF( inputFolder );
        fileManager& tractReader( tractReaderFileMF.getFM() );
        tractReader.readAsLog();
        tractReader.readAsUnThres();
        tractReader.loadMaskImage(maskFilename);

        fileManagerFactory tractWriterFileMF( outputFolder );
        fileManager& tractWriter( tractWriterFileMF.getFM() );
        if ( useFloat )
        {
            tractWriter.writeInFloat();
        }
        else
        {
            tractWriter.writeInChar();
        }

        if ( doZip )
        {
            tractWriter.storeZipped();
        }
        else
        {
            tractWriter.storeUnzipped();
        }

        if (roiMode)
        {
            std::vector<WHcoord> roi, flippedRoi;
            std::vector<size_t> trackids;
            WHcoord datasetSize;
            HC_GRID datasetGrid( HC_SURF );
            size_t numStreamlines( 0 );

            RoiLoader roiLoader( niftiMode );
            roiLoader.readRoi( roiFilename, &datasetGrid, &datasetSize, &numStreamlines, &roi, &trackids );
            flippedRoi = roi;
            if(verbose)
            {
                std::cout<<"Roi loaded, "<< roi.size() <<" seed voxels"<<std::endl;
                std::cout<<"Dataset size is: "+datasetSize.getNameString()+" in "+getGridString(datasetGrid)+" format";
            }

            if (datasetGrid != HC_VISTA && datasetGrid != HC_NIFTI )
            {
                    std::cerr<<"Roi coordinates are neither in nifti nor vista format. Only nifti or vista is supported for this operation"<<std::endl;
                    return 1;
            }

            if (verbose)
            {
                std::cout<<"Flipping and saving leaf tracts to folder: "<< outputFolder <<"..." << std::endl;
            }
            #pragma omp parallel for schedule(guided)
            for (size_t i=0; i<roi.size(); ++i)
            {
                #pragma omp atomic
                    ++tractProg;

                if (roiMode)
                {
                    WHcoord flippedCoord(roi[i]);
                    flippedCoord.m_x= datasetSize.m_x-1-(flippedCoord.m_x);
                    flippedRoi[i] = flippedCoord;
                    std::string tractFilename( tractWriter.getLeafTractFilename(i, trackids, flippedRoi) );

                    if( noOverwrite )
                    {
                        // if tract already exists dont overwrite1
                        if(boost::filesystem::is_regular_file(boost::filesystem::path(tractFilename)))
                            continue;
                    }
                    compactTract tract;
                    tractReader.readLeafTract( i, trackids, roi, &tract );
                    tractReader.flipXtract( &tract );
                    tractWriter.writeLeafTract (i, trackids, flippedRoi, tract );
                }

                #pragma omp atomic
                    ++realDone;

                #pragma omp single nowait // only one thread executes output
                if (verbose) {
                    time_t currentTime(time(NULL));
                    if(currentTime-lastTime>1) {
                        // output progress every 2 seconds
                        lastTime=currentTime;
                        size_t currentCount(tractProg);
                        float progress = (currentCount)*100./(roi.size());
                        int expected_remain(difftime(currentTime,startTime)*((100.-progress)/progress));
                        std::cout<<"\r"<<(int)progress<<" % Completed ("<< currentCount <<" flipped tracts)" <<". Expected remaining time: ";
                        std::cout<<expected_remain/3600 <<"h "<<  (expected_remain%3600)/60 <<"' "<< ((expected_remain%3600)%60) <<"\"  "<< std::flush;
                    }
                }

            } // end for

        }
        else
        {
            std::vector< size_t >nodeVector;
            if( treeBaseMode )
            {
                WHtree thisTree( treeFilename );
                if( verbose )
                {
                    std::cout<< "Tree loaded: " << thisTree.getReport() << std::endl;
                }

                if( treeAllMode )
                {
                    nodeVector.reserve( thisTree.getNumNodes() );
                    for( size_t i=0; i< thisTree.getNumNodes(); ++i )
                    {
                        nodeVector.push_back( i );
                    }
                }
                else
                {

                    if( !thisTree.testRootBaseNodes() )
                    {
                        std::cerr<<"Tree base nodes are not valid meta-leaves (have both leaf and node children)"<<std::endl;
                        return 1;
                    }
                    nodeVector = thisTree.getRootBaseNodes();
                }
            }

            if( nodesMode )
            {
                nodeVector.insert( nodeVector.end(), inNodes.begin(), inNodes.end() );
            }

            std::sort( nodeVector.begin(), nodeVector.end() );
            std::unique( nodeVector.begin(), nodeVector.end() );


            if (verbose)
            {
                std::cout<<"Flipping and saving leaf tracts to folder: "<< outputFolder <<"..." << std::endl;
            }



            #pragma omp parallel for schedule(guided)
            for (size_t i=0; i<nodeVector.size(); ++i)
            {

                #pragma omp atomic
                ++tractProg;

                std::string tractFilename( tractWriter.getNodeTractFilename( nodeVector[i] ) );

                if( noOverwrite )
                {
                    // if tract already exists dont overwrite1
                    if(boost::filesystem::is_regular_file(boost::filesystem::path(tractFilename)))
                        continue;
                }
                compactTract tract;
                tractReader.readNodeTract(nodeVector[i],&tract);
                tractReader.flipXtract( &tract );
                tractWriter.writeNodeTract( nodeVector[i],tract );


                #pragma omp atomic
                ++realDone;

                #pragma omp single nowait // only one thread executes output
                if (verbose)
                {
                    time_t currentTime(time(NULL));
                    if(currentTime-lastTime>1)
                    {
                        // output progress every 2 seconds
                        lastTime=currentTime;
                        size_t currentCount(tractProg);
                        float progress = (currentCount)*100./(nodeVector.size());
                        int expected_remain(difftime(currentTime,startTime)*((100.-progress)/progress));
                        std::cout<<"\r"<<(int)progress<<" % Completed ("<< currentCount <<" flipped tracts)" <<". Expected remaining time: ";
                        std::cout<<expected_remain/3600 <<"h "<<  (expected_remain%3600)/60 <<"' "<< ((expected_remain%3600)%60) <<"\"  "<< std::flush;
                    }
                }

            }
        }

        if (verbose)
        {
            std::cout<<"\r100 % Completed ("<< tractProg <<" flipped tracts)" << std::endl;
            if( noOverwrite )
            {
                std::cout<<"( Only "<< realDone <<" computed on this run, the rest were already present)" << std::endl;
            }

        }

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

