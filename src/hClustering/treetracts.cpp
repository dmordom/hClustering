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
//  treetracts
//
//  Compute the mean tractograms from a hierarchical tree nodes and the original leaf tracts and write them in compact form.
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -t --tree:       File with the hierarchical tree to compute node tractograms from.
//
//   -I --inputf:     Input data folder (containing the seed voxel compact tractograms).
//
//   -O --outputf:    Output folder where tree files will be written.
//
//  [-n --nodes]:     Write tracts for the following node ids (separated with whitespaces).
//
//  [-b --bases]:     Write tracts corresponding to the base-nodes (meta-leaves).
//
//  [-a --all]:       Write tracts for all the tree nodes.
//
//  [-f --full]:      Write full 3D image tracts instead of compact tracts, indicate location of wm mask file here.
//
//  [-c --clustmsk]:  Write for each tract the corresponding 3D mask of all the seed voxels contained in the corresponding cluster, indicate location of wm mask file here.
//
//  [--notracts]:     [use only with -c] Do not write tracts to file (only cluster masks).
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files].
//
//  [-m --cache-mem]: Maximum amount of RAM memory (in GBytes) to use for temporal tractogram cache storing. Valid values [0.1,50]. Default: 0.5.
//
//  [-z --zip]:       zip output files.
//
//  [-F --ufloat]:    use float32 representation to write output tracts (default is uint8).
//
//  [--debugout]:     write additional detailed outputs meant to be used for debugging.
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  example:
//
//  treetracts -t tree_lh.txt -I tracograms/ -O results/ -n 40 65 -b -m 2 -v -c
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
#include "treeManager.h"



int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("treetracts");
        std::string configFilename("../../config/"+progName+".cfg");


        // program parameters
        std::string treeFilename, maskFilename, maskFilename2;
        std::string tractFolder, outputFolder;
        std::vector<size_t> inputNodes;
        unsigned int threads(0);
        bool verbose( false ), debug( false ), niftiMode( true );
        bool bases(false), allNodes( false ), fullTracts( false );
        bool clustMasks( false ), onlyClustMasks( false );
        bool useFloat( false ), doZip( false );
        float memory(0.5);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "tree,t",  boost::program_options::value< std::string >(&treeFilename), "file with the hierarchical tree")
                ( "inputf,I", boost::program_options::value< std::string >(&tractFolder), "Input data folder (containing the seed voxel compact tractograms).")
                ( "outputf,O",      boost::program_options::value< std::string >(&outputFolder), "output folder where tractograms will be written")
                ( "nodes,n", boost::program_options::value< std::vector<size_t> >(&inputNodes)->multitoken(), "[opt] input nodes to compute mean tracts of")
                ( "bases,b", "[opt] write tracts for all base nodes")
                ( "all,a", "[opt] write tracts for all tree nodes")
                ( "full,f", boost::program_options::value< std::string >(&maskFilename), "[opt] write full 3D image tracts instead of compact tracts, must be followed my wm mask filepath")
                ( "clustmsk,c", boost::program_options::value< std::string >(&maskFilename2),  "[opt] write for each mean tract the corresponding 3D mask of contained seed voxels,, must be followed my wm mask filepath")
                ( "notracts", "[opt, use only with -c] do not write tracts to file (write only cluster masks)")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "cache-mem,m",  boost::program_options::value< float >(&memory)->implicit_value(0.5), "[opt] maximum of memory (in GBytes) to use for tractogram cache memory. Default: 0.5." )
                ( "zip,z", "[opt] zip output files.")
                ( "ufloat,F", "[opt] use float32 representation to write tracts (default is uint8)")
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
        posOpt.add("tree", 1);

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
            std::cout << "treetracts" << std::endl << std::endl;
            std::cout << "Compute the mean tractograms from a hierarchical tree nodes and the original leaf tracts and write them in compact form." << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -t --tree:       File with the hierarchical tree to compute node tractograms from." << std::endl << std::endl;
            std::cout << " -I --inputf:     Input data folder (containing the seed voxel compact tractograms)." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where tractogram files will be written." << std::endl << std::endl;
            std::cout << "[-n --nodes]:     Write tracts for the following node ids (separated with whitespaces)." << std::endl << std::endl;
            std::cout << "[-b --bases]:     Write only the tracts corresponding to the base-nodes (meta-leaves)." << std::endl << std::endl;
            std::cout << "[-a --all]:       Write tracts for all the tree nodes." << std::endl << std::endl;
            std::cout << "[-f --full]:      Write full 3D image tracts instead of compact tracts, indicate location of wm mask file here." << std::endl << std::endl;
            std::cout << "[-c --clustmsk]:  Write for each tract the corresponding 3D mask of all the seed voxels contained in the corresponding cluster, indicate location of wm mask file here." << std::endl << std::endl;
            std::cout << "[--notracts]:     [use only with -c] Do not write tracts to file (only cluster masks)." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	     read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << "[-m --cache-mem]: maximum amount of RAM memory (in GBytes) to use for temporal tractogram cache storing. Valid values [0.1,50]. Default: 0.5." << std::endl << std::endl;
            std::cout << "[-z --zip]:       zip output files." << std::endl << std::endl;
            std::cout << "[-F --ufloat]:    use float32 representation to write output tracts (default is uint8)." << std::endl << std::endl;
            std::cout << "[--debugout]:     write additional detailed outputs meant to be used for debugging." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "example:" << std::endl << std::endl;
            std::cout << "treetracts -t tree_lh.txt -I tracograms/ -O results/ -n 40 65 -b -m 2 -v -c" << std::endl << std::endl;
            exit(0);
        }


        if ( variableMap.count( "verbose" ) ) {
            std::cout << "verbose output" << std::endl;
            verbose=true;
        }

        if (variableMap.count("version")) {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }

        if (variableMap.count("all"))
        {
            std::cout << "writing tracts for all tree nodes"<<std::endl;
            allNodes=true;
        }
        else
        {
            if (variableMap.count("bases"))
            {
                std::cout << "writing tracts for all base nodes"<<std::endl;
                bases=true;
            }
            if (variableMap.count("nodes"))
            {
                std::cout << "writing tracts for " << inputNodes.size()<< " tree nodes"<<std::endl;
            }
        }

        if ( variableMap.count("full") )
        {
            if( !boost::filesystem::is_regular_file( boost::filesystem::path( maskFilename ) ) )
            {
                std::cerr << "ERROR: white matter mask file \"" <<maskFilename<< "\" is not a regular file" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            else
            {
                std::cout << "writing tracts in full 3D image form. White mattter mask filepath: " << maskFilename << std::endl;
                fullTracts=true;
            }
        }

        if ( variableMap.count( "clustmsk" ) )
        {
            if( !boost::filesystem::is_regular_file( boost::filesystem::path( maskFilename2 ) ) )
            {
                std::cerr << "ERROR: white matter mask file \"" <<maskFilename2<< "\" is not a regular file" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            else
            {
                clustMasks = true;
                if ( variableMap.count( "notracts" ) )
                {
                    std::cout << "writing only cluster seed voxel masks (no tracts will be written). " << std::flush;
                    onlyClustMasks = true;
                }
                else
                {
                    std::cout << "writing the corresponding cluster seed voxel mask for each tractogram. " << std::flush;
                    verbose=true;
                }

                std::cout << "White mattter mask filepath: " << maskFilename2 << std::endl;

                if ( fullTracts )
                {
                    if( maskFilename2 != maskFilename )
                    {
                        std::cerr << "ERROR: white matter mask files from -f and -c options do not match" << std::endl;
                        exit(-1);
                    }
                }
                maskFilename = maskFilename2;
            }
        }

        if (variableMap.count("ufloat"))
        {
            std::cout << "writing in float"<<std::endl;
            useFloat=true;
        }
        else
        {
            std::cout << "writing in char"<<std::endl;
            useFloat=false;
        }

        if (variableMap.count("zip"))
        {
            std::cout << "zipping output files"<<std::endl;
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
        if ( variableMap.count( "debugout" ) )
        {
            if( verbose )
            {
                std::cout << "Debug output files activated" << std::endl;
            }
            debug = true;
        }

        if (variableMap.count("tree")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename)))
            {
                std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Input tree file: "<< treeFilename << std::endl;
        } else {
            std::cerr << "ERROR: no input tree file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if (variableMap.count("inputf"))
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path(tractFolder))) {
                std::cerr << "ERROR: single tract folder \""<<tractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Single (leaf) tracts folder: "<< tractFolder << std::endl;
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
            std::cout << "Output folder: "<< outputFolder << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no output folder stated"<<std::endl;
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

        if (memory<0.1 || memory>50) {
            std::cerr << "ERROR: cache size must be a positive float between 0.1 and 50"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        std::cout <<"Tractogram cache memory: "<< memory <<" GBytes"<< std::endl;


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
        logFile <<"Tree file:\t"<< treeFilename <<std::endl;
        logFile <<"Tract folder:\t"<< tractFolder <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Processors used:\t"<< threads <<std::endl;
        if(!inputNodes.empty())
        {
            logFile <<"writing tracts for nodes: ";
            for(size_t i=0; i<inputNodes.size(); ++i)
            {
                logFile << inputNodes[i] << " ";
            }
            logFile << std::endl;
        }
        logFile <<"Writing tracts for all bases:\t"<< bases <<std::endl;
        logFile <<"Writing tracts for all nodes:\t"<< allNodes <<std::endl;
        logFile <<"Writing node cluster mask images:\t"<< clustMasks <<std::endl;
        logFile <<"Writing ONLY cluster mask images:\t"<< onlyClustMasks <<std::endl;

        logFile <<"Float32 flag:\t"<< useFloat <<std::endl;
        logFile <<"Zip flag:\t"<< doZip <<std::endl;
        logFile <<"Cache size:\t"<< memory <<" GB"<<std::endl;
        logFile <<"-------------"<<std::endl;


        // ==============================================================================

        WHtree tree(treeFilename);
        if (verbose)
        {
            std::cout<<tree.getReport()<<std::endl;
        }
        logFile <<tree.getReport()<<std::endl;

        logFile <<"-------------"<<std::endl;

        treeManager treeMngr(&tree, verbose );
        treeMngr.log(&logFile);
        treeMngr.setFullTractFolder(outputFolder);
        treeMngr.setSingleTractFolder(tractFolder);
        if( !maskFilename.empty() )
        {
            treeMngr.setMaskFilename( maskFilename );
        }

        if( useFloat )
        {
            treeMngr.writeInFloat();
        }
        else
        {
            treeMngr.writeInChar();
        }
        if( doZip )
        {
            treeMngr.storeZipped();
        }
        else
        {
            treeMngr.storeUnzipped();
        }
        treeMngr.setDebugOutput( debug );

        std::vector<size_t> nodes2write;
        if(allNodes)
        {
            nodes2write.reserve(tree.getNumNodes());
            for(size_t i=0; i<tree.getNumNodes(); ++i)
            {
                nodes2write.push_back(i);
            }
        }
        else
        {
            if( bases )
            {
                if( tree.testRootBaseNodes() )
                {
                    nodes2write = tree.getRootBaseNodes();
                }
                else
                {
                    std::cerr<< "WARNING: tree is not base node tree. Ignoring base nodes..."<<std::endl;
                }
            }

            if( !inputNodes.empty() )
            {
                nodes2write.insert(nodes2write.end(),inputNodes.begin(),inputNodes.end());
                std::sort(nodes2write.begin(), nodes2write.end());
                std::unique(nodes2write.begin(), nodes2write.end());
            }
        }

        if( !onlyClustMasks )
        {
            if( fullTracts )
            {
                treeMngr.writeFullTract( nodes2write );
            }
            else
            {
                treeMngr.setMeanTractFolder(outputFolder);
                if(allNodes)
                {
                    treeMngr.writeAllNodeTracts( memory );
                }
                else
                {
                    treeMngr.writeMeanTracts(nodes2write);
                }
            }
        }
        if(clustMasks)
        {
            treeMngr.writeClusterMasks(nodes2write);
            for( size_t i=0; i < nodes2write.size() ; ++i  )
            {
                WHcoord centre( tree.getMeanCoordinate4node(nodes2write[i]) );
                std::cout << "node " << nodes2write[i] << " -> " << centre.getNameString() << std::endl;
            }
        }

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



