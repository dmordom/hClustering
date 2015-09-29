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
//  comparetrees
//
//  Matches leaves or meta-leaves (base-nodes) across trees and computes tree comparison values (tcpcc and triples).
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   --cl:            [xor with --cg and --cr] direct leaf-wise correspondence. Use for matching trees built over the same seed voxel tractograms.
//                     Due to expected high number of leaves and to reduce computing time triples will be subsampled by 1/LEAF_TRIPLES_FREQ (to change this value modify at source code).
//
//   --cg:            [xor with --cl and --cr] greedy-match base-node-wise correspondence, indicate file where to write/load base-node dissimilarity matrix.
//                     Matches with a dissimilarity higher than DISSIM_THRESHOLD will not be considered a match (to change this value modify at source code).
//
//   --cr:            [xor with --cl and --cg] random base-node-wise correspondence. Used to obtain a random chance baseline for tcpcc and triples value to compare to.
//                     RAND_REPEAT repetitions will be computed and triples will be subsampled by 1/RAND_TRIPLES_FREQ (to change this value modify at source code).
//
//   --t1:            File with first tree to be matched and compared.
//
//   --t2:            File with second tree to be matched and compared.
//
//   --f1:            Folder with the tracts for the first tree. If --cl is chosen the folder should contain leaf tracts.
//                     If --cg or --cr options are chosen, it should contain base-node tracts and cluster masks warped to a common space.
//
//   --f2:            Folder with the tracts for the second tree. If --cl is chosen the folder should contain leaf tracts.
//                     If --cg or --cr options are chosen, it should contain base-node tracts and cluster masks warped to a common space.
//
//   -O --outputf:    Output folder where result files will be written.
//
//  [-t --threshold]: Number of streamlines relative to the total generated that must pass through a tract voxel to be considered for tract similarity
//                     (i.e.: minimum value of a normalized probabilistic tract in natural units to be considered above noise).
//                     Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired. Default: 0 (no threshold).
//
//  [-d --eucdist]:   Maximum euclidean distance (in number of isotorpic voxel distance units) between matched base-node cluster center coordinates to be accepted as a valid match.
//                     Base-nodes considered for match with a higher euclidean distance (in common space) will be considered without match if no better matching possibilities exist.
//                     [use only with --cg or --cr] Default: 20 voxel distance units.
//
//  [-n --noise]:     [use only with --cg] matching-noise correction. insert alpha value (0,1]. Matching noise will not take into account for comparison any tree structure below the noise level.
//                     The noise level for a given node in the tree is computed as the average matching distance of the contained base nodes multiplied by a linear alpha coefficient to control noise weighting.
//                     An alpha value of 0 will compute results at the full [0,1] alpha value range at 0.05 intervals. Refer to (Moreno-Dominguez, 2014) for more information on the matching-noise scheme.
//
//  [--nocomp]:       Only obtain tree correspondence, not the trree comparison values (tcpcc nor triples). Ignored if in --cr mode.
//
//  [--notriples]:    Only obtain the correspondence and tcpcc value, not the triples (the latter is significantly more time-consuming).
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files].
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  example:
//
//  comparetrees -cg distMatrix.nii --t1 tree1.txt --t2 tree2.txt --f1 tracts1/ --f2 tracts2/ -O results/ -t 0.001 -d 20 -v
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
#include "treeComparer.h"
#include "WStringUtils.h"


#define RAND_REPEAT         100
#define RAND_TRIPLES_FREQ   3
#define LEAF_TRIPLES_FREQ   10
#define DISSIM_THRESHOLD    0.9


int main( int argc, char *argv[] )
{



//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("comparetrees");
        std::string configFilename("../../config/"+progName+".cfg");

        typedef enum
        {
            CRSP_DIRECT,
            CRSP_GREEDY,
            CRSP_RAND
        }
        CRSP_MODE;


        // program parameters
        std::string treeFilename1, treeFilename2, tractFolder1, tractFolder2, outputFolder, matrixFilename;
        unsigned int threads(0);
        bool verbose(false), noComp( false ), noTriples( false ), niftiMode( true );
        CRSP_MODE compMode(CRSP_DIRECT);
        bool matchNoise( false ), noiseLoop( false );
        float noiseAlpha( 0 ), maxPhysDist( 20 ), relativeThreshold( 0 );

        // Declare a group of options that will be allowed only on command line
        std::string clmessage( "[xor with --cg and --cr] direct leaf-wise correspondence. Triples will be subsampled by 1/" + string_utils::toString(LEAF_TRIPLES_FREQ) );
        std::string cgmessage( "[xor with --cl and --cr] greedy-match base-node-wise correspondence, indicate file where to write/load base-node dissimilarity matrix. Maximum dissimilarity for a valid match: " + string_utils::toString(DISSIM_THRESHOLD) );
        std::string crmessage( "[xor with --cl and --cg] random base-node-wise correspondence. " + string_utils::toString(RAND_REPEAT) +" repetitions will be computed and triples will be subsampled by 1/" + string_utils::toString(RAND_TRIPLES_FREQ) );


        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "cl", clmessage.c_str())
                ( "cg", boost::program_options::value< std::string >(&matrixFilename), cgmessage.c_str() )
                ( "cr", crmessage.c_str() )
                ( "t1",  boost::program_options::value< std::string >(&treeFilename1), "file with first tree")
                ( "t2",  boost::program_options::value< std::string >(&treeFilename2), "file with second tree")
                ( "f1",  boost::program_options::value< std::string >(&tractFolder1), "folder with the tracts for the first tree")
                ( "f2",  boost::program_options::value< std::string >(&tractFolder2), "folder with the tracts for the second tree")
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder where results will be written")
                ( "threshold,t", boost::program_options::value< float >(&relativeThreshold)->implicit_value(0), "[opt] noise threshold for the tractograms relative to number of streamlines per tract. [0,1)." )
                ( "eucdist,d",  boost::program_options::value< float >(&maxPhysDist)->implicit_value(20), "[opt | use only with --cg or --cr] maximum euclidean distance between cluster centers for a valid match (in number of isotropic voxel distance units ). Default: 20")
                ( "noise,n", boost::program_options::value< float >(&noiseAlpha), "[opt | use only with --cg] matching-noise correction. insert alpha value (0,1]. A value of 0 will compute results at the full [0,1] range at 0.05 intervals ")
                ( "nocomp", "[opt] only obtain tree correspondence, not the trree comparison values (tcpcc nor triples)")
                ( "notriples", "[opt] only obtain the correspondence and tcpcc value, not the triple (the latter is significantly more time-consuming")
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
        //posOpt.add("roi-file", -1);

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
            std::cout << "buildctree" << std::endl << std::endl;
            std::cout << "Build a centroid hierarchical tree from a set of seed voxels and their corresponding tractograms." << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " --cl:            [xor with --cg and --cr] direct leaf-wise correspondence. Use for matching trees built over the same seed voxel tractograms." << std::endl;
            std::cout << "                   Due to expected high number of leaves and to reduce computing time triples will be subsampled by 1/LEAF_TRIPLES_FREQ(1/" << LEAF_TRIPLES_FREQ << ") (to change this value modify at source code)." << std::endl << std::endl;
            std::cout << " --cg:            [xor with --cl and --cr] greedy-match base-node-wise correspondence, indicate file where to write/load base-node dissimilarity matrix." << std::endl;
            std::cout << "                   Matches with a dissimilarity higher than DISSIM_THRESHOLD(" << DISSIM_THRESHOLD << ") will not be considered a match (to change this value modify at source code)." << std::endl << std::endl;
            std::cout << " --cr:            [xor with --cl and --cg] random base-node-wise correspondence. Used to obtain a random chance baseline for tcpcc and triples value to compare to." << std::endl;
            std::cout << "                   RAND_REPEAT(" << RAND_REPEAT << ") repetitions will be computed and triples will be subsampled by 1/RAND_TRIPLES_FREQ(1/" << RAND_TRIPLES_FREQ << ") (to change this value modify at source code)." << std::endl << std::endl;
            std::cout << " --t1:            File with first tree to be matched and compared." << std::endl << std::endl;
            std::cout << " --t2:            File with second tree to be matched and compared." << std::endl << std::endl;
            std::cout << " --f1:            Folder with the tracts for the first tree. If --cl is chosen the folder should contain leaf tracts." << std::endl;
            std::cout << "                   If --cg or --cr options are chosen, it should contain base-node tracts and cluster masks warped to a common space." << std::endl << std::endl;
            std::cout << " --f2:            Folder with the tracts for the second tree. If --cl is chosen the folder should contain leaf tracts." << std::endl;
            std::cout << "                   If --cg or --cr options are chosen, it should contain base-node tracts and cluster masks warped to a common space." << std::endl << std::endl;
            std::cout << " -O --outputf:    output folder where result files will be written." << std::endl << std::endl;
            std::cout << "[-t --threshold]: number of streamlines relative to the total generated that must pass through a tract voxel to be considered for tract similarity" << std::endl;
            std::cout << "                   (i.e.: minimum value of a normalized probabilistic tract in natural units to be considered above noise)." << std::endl;
            std::cout << "                   Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired." << std::endl << std::endl;

            std::cout << "[-d --eucdist]:   Maximum euclidean distance (in number of isotorpic voxel distance units) between matched base-node cluster center coordinates to be accepted as a valid match." << std::endl;
            std::cout << "                   Base-nodes considered for match with a higher euclidean distance (in common space) will be considered without match if no better matching possibilities exist." << std::endl;
            std::cout << "                   [use only with --cg or --cr] Default: 20 voxel distance units." << std::endl << std::endl;
            std::cout << "[-n --noise]:     [use only with --cg] matching-noise correction. insert alpha value (0,1]. Matching noise will not take into account for comparison any tree structure below the noise level." << std::endl;
            std::cout << "                   The noise level for a given node in the tree is computed as the average matching distance of the contained base nodes multiplied by a linear alpha coefficient to control noise weighting." << std::endl;
            std::cout << "                   An alpha value of 0 will compute results at the full [0,1] alpha value range at 0.05 intervals. Refer to (Moreno-Dominguez, 2014) for more information on the matching-noise scheme." << std::endl << std::endl;
            std::cout << "[--nocomp]:       Only obtain tree correspondence, not the trree comparison values (tcpcc nor triples). Ignored if in --cr mode." << std::endl << std::endl;
            std::cout << "[--notriples]:    Only obtain the correspondence and tcpcc value, not the triples (the latter is significantly more time-consuming)." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   Verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	    Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "example:" << std::endl << std::endl;
            std::cout << "comparetrees -cg distMatrix.nii --t1 tree1.txt --t2 tree2.txt --f1 tracts1/ --f2 tracts2/ -O results/ -t 0.001 -d 20 -v" << std::endl << std::endl;
            exit(0);
        }

        if ( variableMap.count( "verbose" ) )
        {
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

        if ( variableMap.count( "version" ) )
        {
            std::cout << progName << ", version 2.0" << std::endl;
            exit(0);
        }



        if ( variableMap.count( "t1" ) )
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename1))) {
                std::cerr << "ERROR: tree1 file \""<<treeFilename1<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tree1 file: "<< treeFilename1 << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tree1 file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if ( variableMap.count( "t2" ) )
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename2))) {
                std::cerr << "ERROR: tree2 file \""<<treeFilename2<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tree2 file: "<< treeFilename2 << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tree1 file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if ( variableMap.count( "f1" ) )
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path(tractFolder1))) {
                std::cerr << "ERROR: Tract folder for tree1 \""<<tractFolder1<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "Tract folder for tree1: "<< tractFolder1 << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tract folder for tree1 stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }

        if ( variableMap.count( "f2") )
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path(tractFolder2))) {
                std::cerr << "ERROR: Tract folder for tree2 \""<<tractFolder2<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "Tract folder for tree2: "<< tractFolder2 << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tract folder for tree2 stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }



        if (  variableMap.count("outputf" ) )
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path(outputFolder))) {
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



        size_t countComp(0);
        if (variableMap.count("cl"))
        {
            std::cout << "Direct leaf-wise correspondence"<<std::endl;
            compMode=CRSP_DIRECT;
            ++countComp;
        }
        if (variableMap.count("cg"))
        {
            std::cout << "Greedy matching baseNode-wise correspondence"<<std::endl;
            compMode=CRSP_GREEDY;
            ++countComp;
        }
        if (variableMap.count("cr"))
        {
            std::cout << "Random baseNode-wise correspondence"<<std::endl;
            compMode=CRSP_RAND;
            ++countComp;
        }


        if (countComp==0)
        {
            std::cerr << "ERROR: no comparison mode stated. Please choose either --cl, --cg or --cr"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if (countComp!=1)
        {
            std::cerr << "ERROR: More than one comparison mode stated, only one mode allowed. Please choose either --cl, --cg or --cr"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("nocomp"))
        {
            std::cout << "Obtaining only tree correspondence"<<std::endl;
            noComp = true;
        }
        else if (variableMap.count("notriples"))
        {
            std::cout << "Obtaining only cpcc and not triples"<<std::endl;
            noTriples = true;
        }

        if (variableMap.count("noise") && variableMap.count("cg"))
        {
            if( noiseAlpha < 0)
            {
                std::cout << "Inserted matching-noise correction alpha value is invalid (negative), not using matching noise correction" <<std::endl;
                matchNoise = false;
                noiseAlpha = 0;
            }
            else if (noiseAlpha > 1)
            {
                std::cout << "Inserted matching-noise correction alpha value is too high ( >1 ), using an alpha of 1" <<std::endl;
                matchNoise = true;
                noiseLoop = false;
                noiseAlpha = 1;
            }
            else if (noiseAlpha == 0)
            {
                std::cout << "Using matching-noise correction. Looping through alpha values from 0 to 1 with 0.05 increases" <<std::endl;
                matchNoise = true;
                noiseLoop = true;
            }
            else
            {
                std::cout << "Using matching-noise correction. Matching noise alpha value: "<<noiseAlpha <<std::endl;
                matchNoise = true;
                noiseLoop = false;
            }
        }
        else if( variableMap.count("noise") && verbose )
        {
            std::cout << "WARNING: matching-noise correction parameter will be ignored when not in greedy base-node correspondence mode" <<std::endl;
        }
        else if( variableMap.count("cg") && verbose )
        {
            std::cout << "NOT Using matching noise correction. ";
            matchNoise = false;
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
        else if ( relativeThreshold == 0 )
        {
            std::cout << "No tractogram thresholding will be applied" << std::endl;
        }
        else if( verbose )
        {
            std::cout << "Tractogram voxels visited by less than " << relativeThreshold * 100 << " % of the streamlines generated will be set to 0 before dissimilarity computation" << std::endl;
        }

        if( verbose )
        {
            std::cout << "Maximum node-base cluster-centers distance to qualify for matching: " << maxPhysDist << " voxels" << std::endl;
        }

        if ( maxPhysDist < 0 )
        {
            std::cerr << "ERROR: negative distance value used, please use a positive value" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        else if ( maxPhysDist == 0 && verbose )
        {
            std::cout << "No Maximum distance restrictions will be applied" << std::endl;
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
        if( niftiMode )
        {
            logFile << "Using nifti file format" << std::endl;
        }
        else
        {
            logFile << "Using vista file format" << std::endl;
        }
        logFile <<"Tree1 file:\t"<< treeFilename1 <<std::endl;
        logFile <<"Tree2 file:\t"<< treeFilename2 <<std::endl;
        logFile <<"Tracts folder for tree1:\t"<< tractFolder1 <<std::endl;
        logFile <<"Tracts folder for tree2:\t"<< tractFolder2 <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile << "Relative threshold:\t" << relativeThreshold << std::endl;
        logFile << "Max matching euclidean distance:\t" << maxPhysDist << " voxels" << std::endl;
        logFile <<"-------------"<<std::endl;


        /////////////////////////////////////////////////////////////////


        WHtree tree1(treeFilename1);
        WHtree tree2(treeFilename2);


        if (!tree1.isLoaded() || !tree2.isLoaded()) {
            throw std::runtime_error ("ERROR @ compareTrees(): trees are not loaded");
        }

        logFile <<"Tree1: "<< tree1.getReport(false) <<std::endl;
        logFile <<"Tree2: "<< tree2.getReport(false) <<std::endl;

        if (tree1.getDataSize() != tree2.getDataSize()) {
            std::cerr <<"Tree1: "<< tree1.getReport() <<std::endl;
            std::cerr <<"Tree2: "<< tree2.getReport() <<std::endl;
            throw std::runtime_error ("ERROR @ compareTrees() datasets to compare have different dimensions");
        }

        if (verbose)
        {
            std::cout <<"Tree1: "<< tree1.getReport(false) <<std::endl;
            std::cout <<"Tree2: "<< tree2.getReport(false) <<std::endl;
        }

        treeComparer comparer(&tree1,&tree2, verbose);
        comparer.log(&logFile);
        comparer.setMaxPhysDist( maxPhysDist );
        comparer.setRelativeThreshold( relativeThreshold );


        if (compMode==CRSP_DIRECT)
        {
            comparer.setCoordsFromFile( false );
            comparer.setMeanTractsFromFile( false );
            comparer.setSingleTractFolder1(tractFolder1);
            comparer.setSingleTractFolder2(tractFolder2);
        }
        else
        {
            std::cout <<comparer.reportBaseNodes()<<std::endl;
            if (!comparer.areRealBaseNodes())
            {
                throw std::runtime_error ("ERROR @ compareTrees(): base nodes are of mixed type (contain both leaves and nodes), cannot compute base-node-wise matching");
            }
            comparer.setCoordsFromFile( true );
            comparer.setMeanTractsFromFile( true );
            comparer.setMeanTractFolder1(tractFolder1);
            comparer.setMeanTractFolder2(tractFolder2);
        }



        ///////////////////////////////////////////////////


        if (compMode == CRSP_RAND)
        {
            logFile <<"Correspondance mode:\t Random"<<std::endl;

            const size_t randRepetitions( RAND_REPEAT );
            const size_t randTriplesFreq( RAND_TRIPLES_FREQ );

            std::cout<<std::endl<< randRepetitions<<" rep loop: "<<std::endl;
            std::string outCpctFilename(outputFolder+"/randCpct.txt");
            std::string outStriplesFilename(outputFolder+"/randStriples.txt");
            std::ofstream outCpctFile(outCpctFilename.c_str());
            if(!outCpctFile) {
                std::cerr << "ERROR: unable to open out file: \""<<outCpctFilename<<"\""<<std::endl;
                exit(-1);
            }
            std::ofstream outStriplestFile(outStriplesFilename.c_str());
            if(!outStriplestFile) {
                std::cerr << "ERROR: unable to open out file: \""<<outStriplesFilename<<"\""<<std::endl;
                exit(-1);
            }

            std::cout<<"Fetching nodes and coordinates"<<std::endl;
            comparer.fetchBaseNodes( true );

            std::cout<<"Running random matchings"<<std::endl;


            outCpctFile <<  "weighted-tCPCC    simple-sCPCC   %usedPairs   effectiveGranularity" << std::endl;

            for(int i=0; i<randRepetitions; ++i) {

                WHtree treeRand1(tree1);
                WHtree treeRand2(tree2);
                treeComparer randComparer(&treeRand1,&treeRand2,comparer);

                randComparer.randomCorrespondence();
                std::pair< std::pair< float, float >, std::pair< float, float > >  tCPCC(randComparer.doTcpcc());
                outCpctFile <<  tCPCC.first.first << " " <<  tCPCC.first.second << " " << std::flush;
                outCpctFile <<  tCPCC.second.first*100 << " " <<  tCPCC.second.second << std::endl;


                if( !noTriples )
                {
                    std::pair< float, float >sTriplets(randComparer.simpleTriplets(randTriplesFreq));
                    outStriplestFile << sTriplets.first << " " << sTriplets.second << std::endl;
                }
            }

        }
        else
        {

            size_t sampleFreq(1);

            bool redoCoords( true );
            std::string outputFileName(outputFolder+"/compValues.txt");
            std::ofstream outFile(outputFileName.c_str());
            if(!outFile) {
                std::cerr << "ERROR: unable to open out file: \""<<outputFileName<<"\""<<std::endl;
                exit(-1);
            }
            std::string outputCompactFileName(outputFolder+"/compactValues.txt");
            std::ofstream outCompactFile(outputCompactFileName.c_str());
            if(!outCompactFile) {
                std::cerr << "ERROR: unable to open out file: \""<<outputCompactFileName<<"\""<<std::endl;
                exit(-1);
            }

            std::vector<float> correspRates;

            switch (compMode)
            {
            case CRSP_DIRECT:
                logFile <<"Correspondance mode:\t Direct"<<std::endl;
                if (verbose)
                {
                    std::cout<<"Equalizing leaves.."<<std::endl;
                }
                if(!comparer.leafCorrespondence())
                {
                    std::cout<<"Tree coordinates match, no changes made."<<std::endl;
                }

                sampleFreq = LEAF_TRIPLES_FREQ;
                break;

            case CRSP_GREEDY:
                logFile <<"Correspondance mode:\t Greedy"<<std::endl;
                if( boost::filesystem::is_regular_file( boost::filesystem::path( matrixFilename ) ) )
                {
                    if (verbose)
                    {
                        std::cout <<"Getting zipped similarity matrix from file: "<< matrixFilename << " ..." << std::flush;
                    }
                    comparer.readBaseDistMatrix( matrixFilename );
                    if (verbose)
                    {
                        std::cout <<"Done"<< std::endl;
                    }
                    logFile <<"Similarity matrix read from file: " << matrixFilename << std::endl;

                }
                else if ( boost::filesystem::is_regular_file( boost::filesystem::path( matrixFilename + ".gz" ) ) )
                {
                    if (verbose)
                    {
                        std::cout <<"Getting zipped similarity matrix from file: "<<  matrixFilename + ".gz" << " ..." << std::flush;
                    }
                    comparer.readBaseDistMatrix( matrixFilename + ".gz" );
                    if (verbose)
                    {
                        std::cout <<"Done"<< std::endl;
                    }
                    logFile <<"Similarity matrix read from file: " << matrixFilename << ".gz" << std::endl;
                }
                else
                {
                    if (verbose)
                    {
                        std::cout <<"Similarity matrix file does not exist, it will be computed and saved in file"<< std::endl;
                    }
                    comparer.getBaseDistMatrix();
                    comparer.writeBaseDistMatrix( matrixFilename );
                    if (verbose)
                    {
                        std::cout << "Similarity matrix saved in file: " << matrixFilename << std::endl;
                    }
                    logFile <<"Similarity matrix saved in file: " << matrixFilename << std::endl;

                }
                if (noComp)
                {
                    redoCoords = false;
                }


                comparer.greedyCorrespondence( DISSIM_THRESHOLD, redoCoords );
                comparer.writeFullCorrespondence( outputFolder + "/fullCorresp.txt" );
                comparer.writeCorrespondence( outputFolder + "/correspTable.txt" );

                correspRates = ( comparer.rateCorrespondence() );
                if( correspRates.size() != 6 )
                {
                    std::cerr << "ERROR: correspondance ratings vector has wrong size" << std::endl;
                }
                else
                {

                    outFile << "Size-Match_Correlation: " << correspRates[0] << std::endl;
                    logFile << "Size-Match_Correlation: " << correspRates[0] << std::endl;

                    outFile << "Mean_Match_Distance: " << correspRates[1] << std::endl;
                    logFile << "Mean_Match_Distance: " << correspRates[1] << std::endl;

                    outFile << "Size-Weighted_Match_Distance: " << correspRates[2] << std::endl;
                    logFile << "Size-Weighted_Match_Distance: " << correspRates[2] << std::endl;
                    outCompactFile << "Match_Distance: " << correspRates[2] << std::endl;

                    outFile << "%_of_matches: " << 100.0 * correspRates[3] << std::endl;
                    logFile << "%_of_matches: " << 100.0 * correspRates[3] << std::endl;

                    outFile << "Mean_Euclidean_Distance: " << correspRates[4] << std::endl;
                    logFile << "Mean_Euclidean_Distance: " << correspRates[4] << std::endl;

                    outFile << "Size-Weighted_Euclidean_Distance: " << correspRates[5] << std::endl;
                    logFile << "Size-Weighted_Euclidean_Distance: " << correspRates[5] << std::endl;
                    outCompactFile << "Euclidean_Distance: " << correspRates[5] << std::endl;

                }

                break;

            default:
                std::cerr << "ERROR: unrecognized correspondance mode"<<std::endl;
                exit(-1);
                break;
            } // end switch



            if( matchNoise )
            {

                std::string outpuAlphatFileName(outputFolder+"/alphaTest.txt");
                std::ofstream outAlpha(outpuAlphatFileName.c_str());
                if(!outAlpha) {
                    std::cerr << "ERROR: unable to open out file: \""<<outpuAlphatFileName<<"\""<<std::endl;
                    exit(-1);
                }

                outAlpha << "Alpha gran1 gran2 meanGran tCPCC %pairs effGran: "<<std::endl;
                float alphaStart( 0 ), alphaEnd( 0 );
                if( noiseLoop )
                {
                    if( verbose )
                    {
                        std::cout<< "Starting alpha loop "<< std::endl<< std::endl;
                    }
                    alphaStart = 0;
                    alphaEnd = 1.01;
                }
                else
                {
                    alphaStart = noiseAlpha;
                    alphaEnd = noiseAlpha;
                }


                for( float thisNoiseAlpha = alphaStart; thisNoiseAlpha <= alphaEnd; thisNoiseAlpha+= 0.05 )
                {

                    outAlpha << thisNoiseAlpha << " " << std::flush;
                    WHtree treeloop1(tree1);
                    WHtree treeloop2(tree2);
                    treeComparer comparerLoop(&treeloop1, &treeloop2, comparer);

                    outFile << "matching_noise_alpha_value: "<<thisNoiseAlpha<<std::endl;
                    logFile << "matching_noise_alpha_value: "<<thisNoiseAlpha<<std::endl;


                    std::pair< float, float > matchGran( comparerLoop.applyNoiseBaseline( thisNoiseAlpha ) );

                    outFile << "Size_of_maxgran_part_1: " << matchGran.first << std::endl;
                    outFile << "Size_of_maxgran_part_2: " << matchGran.second << std::endl;
                    outFile << "Average_maxgran_size: "<< ( ( matchGran.first + matchGran.second )  / 2.0 ) << std::endl;

                    logFile << "Size_of_maxgran_part_1: " << matchGran.first << std::endl;
                    logFile << "Size_of_maxgran part_2: " << matchGran.second << std::endl;
                    logFile << "Average_maxgran_size: "<< ( ( matchGran.first + matchGran.second )  / 2.0 ) << std::endl;

                    outAlpha << matchGran.first << " " << matchGran.second << " " << ( matchGran.first + matchGran.second )  / 2.0  << " " << std::flush;


                    if( !noComp )
                    {
                        std::pair< std::pair< float, float >, std::pair< float, float > > tCPCC(comparerLoop.doTcpcc());

                        outFile << "weighted_tCPCT: " << tCPCC.first.first << std::endl;
                        logFile << "weighted_tCPCC: " << tCPCC.first.first << std::endl;

                        outFile << "simple_CPCC: " << tCPCC.first.second << std::endl;
                        logFile << "simple_CPCC: " << tCPCC.first.second << std::endl;
                        outFile << "%_used_pairs: " << 100*tCPCC.second.first << std::endl;
                        logFile << "%_used_pairs: " << 100*tCPCC.second.first << std::endl;
                        outFile << "effect_Granularity: " << tCPCC.second.second << std::endl;
                        logFile << "effect_Granularity: " << tCPCC.second.second << std::endl;

                        outAlpha << tCPCC.first.first << " " << 100*tCPCC.second.first << " " << tCPCC.second.second  << " " << std::flush;


                        if( !noTriples )
                        {
                            std::pair< float, float > sTriplets(comparerLoop.simpleTriplets(sampleFreq));

                            outFile << "Simple_Triplets_Unweighted: " << sTriplets.first << std::endl;
                            logFile << "Simple_Triplets_Unweighted: " << sTriplets.first << std::endl;

                            outFile << "Simple_Triplets_Size-Weighted: " << sTriplets.second << std::endl;
                            logFile << "Simple_Triplets_Size-Weighted: " << sTriplets.second << std::endl;
                        }
                    }


                    logFile << "Final Tree 1: " << treeloop1.getReport( false ) << std::endl;
                    logFile << "Final Tree 2: " << treeloop2.getReport( false ) << std::endl;
                    logFile << comparerLoop.reportBaseNodes() << std::endl;



                    if (verbose) {
                        std::cout << "Final Tree 1: " << treeloop1.getReport( false ) << std::endl;
                        std::cout << "Final Tree 2: " << treeloop2.getReport( false ) << std::endl;
                        std::cout << comparerLoop.reportBaseNodes() << std::endl;
                    }

                    if(verbose)
                    {
                        std::cout<< std::endl;
                    }

                    outAlpha << std::endl;

                } // end alpha for
            } // end matchnoise if
            else
            {
                if( !noComp )
                {
                    std::pair< std::pair< float, float >, std::pair< float, float > > tCPCC(comparer.doTcpcc());

                    outFile << "weighted_tCPCT: " << tCPCC.first.first << std::endl;
                    logFile << "weighted_tCPCC: " << tCPCC.first.first << std::endl;
                    outCompactFile << "tCPCT: " << tCPCC.first.first << std::endl;

                    outFile << "simple_CPCC: " << tCPCC.first.second << std::endl;
                    logFile << "simple_CPCC: " << tCPCC.first.second << std::endl;
                    outFile << "%_used_pairs: " << 100*tCPCC.second.first << std::endl;
                    logFile << "%_used_pairs: " << 100*tCPCC.second.first << std::endl;
                    outFile << "effect_Granularity: " << tCPCC.second.second << std::endl;
                    logFile << "effect_Granularity: " << tCPCC.second.second << std::endl;


                    if( !noTriples )
                    {
                        std::pair< float, float > sTriplets(comparer.simpleTriplets(sampleFreq));

                        outFile << "Simple_Triplets_Unweighted: " << sTriplets.first << std::endl;
                        logFile << "Simple_Triplets_Unweighted: " << sTriplets.first << std::endl;
                        outFile << "Simple_Triplets_Size-Weighted: " << sTriplets.second << std::endl;
                        logFile << "Simple_Triplets_Size-Weighted: " << sTriplets.second << std::endl;
                        outCompactFile << "wTriples: " << sTriplets.second << std::endl;
                    }
                }// end matchnoise else


                logFile << "Final Tree 1: " << tree1.getReport( false ) << std::endl;
                logFile << "Final Tree 2: " << tree2.getReport( false ) << std::endl;
                logFile << comparer.reportBaseNodes() << std::endl;



                if (verbose) {
                    std::cout << "Final Tree 1: " << tree1.getReport( false ) << std::endl;
                    std::cout << "Final Tree 2: " << tree2.getReport( false ) << std::endl;
                    std::cout << comparer.reportBaseNodes() << std::endl;
                }

                if(verbose)
                {
                    std::cout<< std::endl;
                }
            }

            tree1.writeTree( outputFolder + "/treeCompared1.txt" );
            tree2.writeTree( outputFolder + "/treeCompared2.txt" );

        } // end no-rand if

        /////////////////////////////////////////////////////////////////

        // save and print total time
        time_t programEndTime( time( NULL ) );
        int totalTime( difftime( programEndTime, programStartTime ) );
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

