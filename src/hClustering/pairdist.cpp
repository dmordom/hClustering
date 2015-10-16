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
//  pairdist
//
//  Retrieve the distance (dissimilarity) value between two tree leaves or nodes as encoded in the corresponding hierarchical tree.
//   Additionally distance between leaves can be retrieved from a distance matrix, and those from leaves/nodes computed directly from leaf/node tractograms.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -t --tree:       File with the hierarchical tree.
//
//   -i --IDs:        Input node IDs to compute the distance from, insert a pair of values, one for each node ID.
//
//  [-l --leaves]:    Interpret input node IDs as leaf IDs.
//
//  [-T --threshold]: Threshold to apply directly to the tractogram values before computing the dissimilarity (in order to avoid tractography noise affect the result).
//                     Unlike in other hClustering commands, this threshold value is an absolute value to apply to the tractogram data as is, not a relative threshold.
//                     Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired.
//
//   -L --leaftractf: Folder with the leaf seed voxel probabilistic tracts. Will trigger direct computation of tractogram distance (and prior computation of mean tractograms in case of node IDS).
//                     Tracts must be normalized.
//
//   -N --nodetractf: Folder with the node mean tracts. Tracts must be normalized. Do not use together with -leaves option.
//
//   -M --matrixf:    Folder with the dissimilarity matrix files. Use only together with -leaves option.
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files.
//
//
//  * Usage example:
//
//   pairdist -t tree.txt -i 234 368 -T 0.4 -L leaftracts/ -N nodetracts/ -M matrix/
//
//  * Outputs:
//
//   Results are displayed on standard output (screen).
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
#include "WHtree.h"
#include "compactTract.h"
#include "compactTractChar.h"
#include "fileManagerFactory.h"
#include "distBlock.h"
#include "treeManager.h"
#include "WStringUtils.h"



int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        //boost::filesystem::path workingDir( boost::filesystem::current_path());

        // ========== PROGRAM PARAMETERS ==========

        std::string progName("pairdist");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string treeFilename;
        std::string leafTractFolder,nodeTractFolder, distMatrixFolder;
        std::vector<size_t> inNodes;
        bool byLeafTract(false), byNodeTract(false), byMatrix(false), areLeaves( false );
        size_t idA( 0 ), idB( 0 );
        float thresValue(0);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "tree,t",  boost::program_options::value< std::string >(&treeFilename), "tree file")
                ( "IDs,i", boost::program_options::value< std::vector< size_t > >(&inNodes)->multitoken(), "input node IDs to compute the distance from, insert a pair of values, one for each node ID")
                ( "leaves,l", "[opt] interpret input nodes as leaf IDs")
                ( "threshold,T", boost::program_options::value<float> (&thresValue)->implicit_value(0), "[opt] threshold to apply before dissimilarity computation. Default 0 (no threshold). Use only for options -L and -N")
                ( "leaftractf,L", boost::program_options::value< std::string >(&leafTractFolder), "[opt] folder with the leaf seed voxel probabilistic tracts. Tracts must be normalized")
                ( "nodetractf,N", boost::program_options::value< std::string >(&nodeTractFolder), "[opt] folder with the node mean tracts. Tracts must be normalized")
                ( "matrixf,M",  boost::program_options::value< std::string >(&distMatrixFolder), "[opt] folder with the distance matrix files")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "vista", "[opt] Write output tree in vista coordinates (default is nifti)." )
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
        //posOpt.add("rest", -1);


        boost::program_options::variables_map variableMap;
        store(boost::program_options::command_line_parser(argc, argv).options(cmdlineOptions).positional(posOpt).run(), variableMap);

        std::ifstream ifs(configFilename.c_str());
        store(parse_config_file(ifs, configFileOptions), variableMap);
        notify(variableMap);






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
            std::cout << "pairdist" << std::endl << std::endl;
            std::cout << "Retrieve the distance (dissimilarity) value between two tree leaves or nodes as encoded in the corresponding hierarchical tree." << std::endl;
            std::cout << " Additionally distance between leaves can be retrieved from a distance matrix, and those from leaves/nodes computed direclty from leaf/node tractograms." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -t --tree:       File with the hierarchical tree." << std::endl << std::endl;
            std::cout << " -i --IDs:        Input node IDs to compute the distance from, insert a pair of values, one for each node ID." << std::endl << std::endl;
            std::cout << "[-l --leaves]:    Interpret input node IDs as leaf IDs." << std::endl << std::endl;
            std::cout << "[-T --threshold]: Threshold to apply directly to the tractogram values before computing the dissimilarity (in order to avoid tractography noise affect the result)." << std::endl;
            std::cout << "                   Unlike in other hClustering commands, this threshold value is an absolute value to apply to the tractogram data as is, not a relative threshold." << std::endl;
            std::cout << "                   Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired." << std::endl << std::endl;
            std::cout << " -L --leaftractf: Folder with the leaf seed voxel probabilistic tracts. Will trigger direct computation of tractogram distance (and prior computation of mean tractograms in case of node IDS)." << std::endl;
            std::cout << "                   Tracts must be normalized." << std::endl << std::endl;
            std::cout << " -N --nodetractf: Folder with the node mean tracts. Tracts must be normalized. Do not use together with -leaves option." << std::endl << std::endl;
            std::cout << " -M --matrixf:    Folder with the dissimilarity matrix files. Use only together with -leaves option." << std::endl << std::endl;
            std::cout << "[--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " pairdist -t tree.txt -i 234 368 -T 0.4 -L leaftracts/ -N nodetracts/ -M matrix/" << std::endl << std::endl;
            std::cout << "* Outputs:" << std::endl << std::endl;
            std::cout << " Results are displayed on standard output (screen)." << std::endl;
            std::cout << std::endl;
            exit(0);
        }

        if (variableMap.count( "version" ) )
        {
            std::cout << progName << ", version 2.0" <<std::endl;
            exit(0);
        }


        if ( variableMap.count( "vista" ) )
        {
            std::cout << "Using vista format" << std::endl;
            fileManagerFactory fmf;
            fmf.setVista();
        }
        else
        {
            std::cout << "Using nifti format" << std::endl;
            fileManagerFactory fmf;
            fmf.setNifti();
        }

        if (variableMap.count( "tree" ) )
        {
            if( !boost::filesystem::is_regular_file(boost::filesystem::path( treeFilename ) ) )
            {
                std::cerr << "ERROR: tree file \"" << treeFilename << "\" is not a regular file" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tree file: "<< treeFilename << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tree file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if ( variableMap.count("IDs") )
        {
            if(inNodes.size()!=2) {
                std::cerr << "ERROR: 2 node/leaf IDs must be entered" <<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            idA = inNodes[0];
            idB = inNodes[1];
            std::cout<< "Input IDs: " << idA << " and " << idB << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no input node/leaf IDs stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if ( variableMap.count("leaves") )
        {
            std::cout<< "Interpreting inputs as leaf IDs" << std::endl;
            areLeaves = true;
        }
        else
        {
            std::cout<< "Interpreting inputs as node IDs" << std::endl;
            areLeaves = false;
        }

        if ( variableMap.count( "matrixf" ) )
        {
            if( !areLeaves )
            {
                std::cerr << "WARNING: using node IDs, distance matrix input will be ignored"<<std::endl;
                byMatrix=false;
            }
            else
            {
                if( !boost::filesystem::is_directory( boost::filesystem::path( distMatrixFolder ) ) )
                {
                    std::cerr << "ERROR: distance matrix folder \""<<distMatrixFolder<<"\" is not a directory"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }
                std::cout << "Distance matrix folder: "<< distMatrixFolder << std::endl;
                byMatrix=true;
            }
        }

        if ( variableMap.count( "leaftractf" ) )
        {
            if( !boost::filesystem::is_directory(boost::filesystem::path( leafTractFolder ) ) )
            {
                std::cerr << "ERROR: leaf seed tract folder \""<<leafTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Leaf seed tract folder: "<< leafTractFolder << std::endl;
            byLeafTract=true;
        }

        if ( variableMap.count( "nodetractf" ) )
        {
            if( areLeaves )
            {
                std::cerr << "WARNING: using leaf IDs, mean node tract folder input will be ignored"<<std::endl;
                byNodeTract=false;
            }
            else
            {
                if( !boost::filesystem::is_directory( boost::filesystem::path( nodeTractFolder ) ) )
                {
                    std::cerr << "ERROR: node mean tract folder \""<<nodeTractFolder<<"\" is not a directory"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }
                std::cout << "Node mean tract folder: "<< nodeTractFolder << std::endl;
                byNodeTract=true;
            }
        }

        if (variableMap.count( "threshold" ) )
        {
            if( byNodeTract || byLeafTract )
            {
                if(thresValue < 0 || thresValue >= 1)
                {
                    std::cerr << "ERROR: tthreshold must be [0,1)"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }
                std::cout << "Tractogram threshold: "<< thresValue << std::endl;
            }
            else
            {
                std::cerr << "WARNING: Not using tractogram sources (option -L or -N), threshold input will be ignored"<<std::endl;
            }
        }
        else
        {
            if( byNodeTract || byLeafTract )
            {
                std::cout << "No tractogram threshold will be applied" << std::endl;
            }
        }



        // ========== OBTAIN DISTANCES ==========

        WHtree tree( treeFilename );
        std::cout << tree.getReport() << std::endl;

        size_t id1( inNodes[0] ), id2( inNodes[1] );
        nodeID_t fullID1, fullID2;
        if( areLeaves )
        {
            fullID1 = std::make_pair< bool, size_t >( false, id1 );
            fullID2 = std::make_pair< bool, size_t >( false, id2 );
            std::cout << "Seed A. ID: " << id1 << ". Coords: " << tree.getCoordinate4leaf( id1 ).getNameString() << ". Trackid: " << tree.getTrackID( id1 ) << std::endl;
            std::cout << "Seed B. ID: " << id2 << ". Coords: " << tree.getCoordinate4leaf( id2 ).getNameString() << ". Trackid: " << tree.getTrackID( id2 ) << std::endl;

        }
        else
        {
            fullID1 = std::make_pair< bool, size_t >( true, id1 );
            fullID2 = std::make_pair< bool, size_t >( true, id2 );
        }

        std::cout << std::endl;

        dist_t copheneticDist( tree.getDistance( fullID1, fullID2 ) );
        std::cout << "Cophenetic distance:\t" << string_utils::toString( copheneticDist ) << std::endl;

        if(byMatrix)
        {
            distBlock distanceBlock( distMatrixFolder );
            WHcoord coord1( tree.getCoordinate4leaf( id1 ) ), coord2( tree.getCoordinate4leaf( id2 ) );
            distanceBlock.loadBlock( coord1, coord2 );
            dist_t matrixDist( distanceBlock.getDistance( coord1, coord2 ) );
            std::cout<< "Matrix distance:\t" << string_utils::toString( matrixDist ) << std::endl;
        }

        if( byNodeTract )
        {
            fileManagerFactory nodeFileMF( nodeTractFolder );
            fileManager& nodeFM( nodeFileMF.getFM() );
            compactTract tract1, tract2;
            nodeFM.readAsLog();
            nodeFM.readAsUnThres();
            nodeFM.readNodeTract( id1, &tract1 );
            nodeFM.readNodeTract( id2, &tract2 );
            tract1.threshold( thresValue );
            tract2.threshold( thresValue );
            tract1.computeNorm();
            tract2.computeNorm();
            dist_t nodeDist( tract1.tractDistance( tract2 ) );
            std::cout<< "Distance by node tracts:\t" << string_utils::toString( nodeDist ) << std::endl;
        }

        if( byLeafTract )
        {
            if( areLeaves )
            {
                fileManagerFactory leafFileMF( leafTractFolder );
                fileManager& leafFM( leafFileMF.getFM() );
                compactTract tract1, tract2;
                leafFM.readAsLog();
                leafFM.readAsUnThres();
                leafFM.readLeafTract( id1, tree.getTrackids(), tree.getRoi(), &tract1 );
                leafFM.readLeafTract( id2, tree.getTrackids(), tree.getRoi(), &tract2 );
                tract1.threshold( thresValue );
                tract2.threshold( thresValue );
                tract1.computeNorm();
                tract2.computeNorm();
                dist_t leafDist( tract1.tractDistance( tract2 ) );
                std::cout<< "Distance by leaf tracts:\t" << string_utils::toString( leafDist ) << std::endl;
            }
            else
            {
                treeManager treeMgr(&tree, true);
                treeMgr.setSingleTractFolder( leafTractFolder );
                compactTract meanTract1( treeMgr.getMeanTract( id1 ) );
                compactTract meanTract2( treeMgr.getMeanTract( id2 ) );
                meanTract1.threshold( thresValue );
                meanTract2.threshold( thresValue );
                meanTract1.computeNorm();
                meanTract2.computeNorm();
                dist_t meanDist( meanTract1.tractDistance( meanTract2 ) );
                std::cout<< "Distance by averaged leaf tracts:\t" << string_utils::toString( meanDist ) << std::endl;
            }
        }

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

