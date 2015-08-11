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
//  matchpartition
//
// Finds the best matching corresponding partitions in a target tree to those present in an unrelated reference tree (meta-leaf matching across these two trees must have been precomputed using comparetrees).
//  Two partition matching algorithms are available: signature matching and overlap matching. Found target partitions will be color-matched as best as possible.
//  There is also the possibility of only color-matching predefined partitions of the target tree to predefined partitions of the reference tree.
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -r --reference:  The tree file with the reference partitioned tree.
//
//   -t --target:     The tree file with the target tree to find matching partitions in (or with partitions to be color-matched).
//
//   -m --leafmatch   File with the meta-leaf matching information across both trees (output of comparetrees command).
//
//   -O --outputf:    Output folder where partition/color matched tree files will be written.
//
//  [-s --signature]: Signature-based partition matching, instert lambda coefficient value. [xor with -o and -c].
//                     In this method a pair signature matrices are computed for each reference-target partitions to find the quality of the match.
//                     Each signature matrix defines a value for each pair of base-nodes of the tree it belongs to: 1 the base nodes are found in the same cluster, 0 if otherwise.
//                     The higher the correlation between the reference and target-derived matrices, the best match is the target tree partition to the reference tree one.
//                     A smart hierarchical search through possible partritions is conducted to find the one with best signature matching.
//                    The lambda coefficient determines if and how a similar number of clusters in both partitions affects the matching quality value,
//                     Lambda=0 -> cluster number does not affect the quality value. Lambda=1 -> cluster value similarity has as much weight as singature correlation.
//
//  [-o --overlap]:   Overlap-based partition matching. [xor with -o and -c].
//                     A match between two partititionsis found by iteratively matching clusters with higher base-node overlap and resolving possible ambiguities.
//                     The matching quality between parittions is defined as the number of base-nodes pairs that are classified in the same way in both partitions
//                     (both in the smae cluster r both in different clusters) against the total number of pair combinations.
//                     A smart hierarchical search through possible partritions is conducted to find the one with best signature matching.
//
//  [-d --depth]:     Partition search depth (for signature and overlap matching. A higher value will mean a more exhaustive search of the possible partitions,
//                     but also a higher computation time, specially if the partition to be matched has a high number of clusters (>100).
//                     The default value (0, recommended) will adaptively give high search depth to low-cluster partitions and lower search depth to high-cluster partittions.
//
//  [-c --justcolor]: Perform only color matching across reference and target tree parttitions (both trees need to have the same number of precompouted partitions).
//                     In multiple-to-one matching cases clusters from the reference tree might also be recolored to better identify matching relationships across partitions.
//
//  [-x --excl]:      Color exclusively clusters that have a match, clusters without match will be recolored white (on both reference and target trees)
//
//  [-v --verbose]:   Verbose output (recommended).
//
//  [--vista]:        Write output tree files in vista coordinates (default is nifti).
//
//
//  example:
//
//  matchpartition -r refTree.txt -t targetTree.txt -m matching.txt -O results/ -s 0.5 -v
//    -or-
//  matchpartition -r refTree.txt -t targetTree.txt -m matching.txt -O results/ -o -v
//    -or-
//  matchpartition -r refTree.txt -t targetTree.txt -m matching.txt -O results/ -c -v
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
#include "partitionMatcher.h"
#include "WStringUtils.h"
#include "fileManagerFactory.h"

int main( int argc, char *argv[] )
{
//    try {
        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());

        // ========== PROGRAM PARAMETERS ==========

        std::string progName("matchpartition");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string refTreeFilename, targetTreeFilename, matchTableFilename, outputFolder;
        unsigned int searchDepth(1);
        float lambda(0);
        bool signaturePart(false), colorMatching(false), overlapPart(false), exclusive(false);
        bool verbose(false), niftiMode( true );


        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "reference,r",  boost::program_options::value< std::string >(&refTreeFilename), "file with reference partitioned tree" )
                ( "target,t",  boost::program_options::value< std::string >(&targetTreeFilename), "file with target tree to be partitioned-matched" )
                ( "leafmatch,m",  boost::program_options::value< std::string >(&matchTableFilename), "file with meta-leaves (base-nodes) matching table" )
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder where partition-matched trees will be written" )
                ( "signature,s",  boost::program_options::value< float >(&lambda), "[xor with -o and -c] Signature-based partition matching, inster lambda coefficient value" )
                ( "overlap,o", "[xor with -s and -c] Meta-leaf overlap-based partition matching")
                ( "depth,d", boost::program_options::value< unsigned int >(&searchDepth)->implicit_value(0), "[opt] partition search depth. Default: 0 (automatic partition-size based adaptive depth, recommended)")
                ( "justcolor,c",  "[xor with -s and -o] Perform only olor matching (requires pre-computed partitions in both trees)")
                ( "excl,x",  "[opt] color exclusively clusters that have a match, clusters without match will be white")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] Write output tree in vista coordinates (default is nifti)." )
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

        if (variableMap.count( "help" ) )
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
            std::cout << "matchpartition" << std::endl << std::endl;
            std::cout << "Finds the best matching corresponding partitions in a target tree to those present in an unrelated reference tree (meta-leaf matching across these two trees must have been precomputed using comparetrees)." << std::endl;
            std::cout << " Two partition matching algorithms are available: signature matching and overlap matching. Found target partitions will be color-matched as best as possible." << std::endl;
            std::cout << " There is also the possibility of only color-matching predefined partitions of the target tree to predefined partitions of the reference tree." << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -r --reference:  The tree file with the reference partitioned tree." << std::endl << std::endl;
            std::cout << " -t --target:     The tree file with the target tree to find matching partitions in (or with partitions to be color-matched)." << std::endl << std::endl;
            std::cout << " -m --leafmatch   File with the meta-leaf matching information across both trees (output of comparetrees command)." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where partitioned/color matched tree files will be written." << std::endl << std::endl;
            std::cout << "[-s --signature]: Signature-based partition matching, instert lambda coefficient value. [xor with -o and -c]." << std::endl;
            std::cout << "                   In this method a pair signature matrices are computed for each reference-target partitions to find the quality of the match." << std::endl;
            std::cout << "                   Each signature matrix defines a value for each pair of base-nodes of the tree it belongs to: 1 the base nodes are found in the same cluster, 0 if otherwise." << std::endl;
            std::cout << "                   The higher the correlation between the reference and target-derived matrices, the best match is the target tree partition to the reference tree one." << std::endl;
            std::cout << "                   A smart hierarchical search through possible partritions is conducted to find the one with best signature matching." << std::endl;
            std::cout << "                  The lambda coefficient determines if and how a similar number of clusters in both partitions affects the matching quality value," << std::endl;
            std::cout << "                   Lambda=0 -> cluster number does not affect the quality value. Lambda=1 -> cluster value similarity has as much weight as singature correlation." << std::endl << std::endl;
            std::cout << "[-o --overlap]:   Overlap-based partition matching. [xor with -o and -c]." << std::endl;
            std::cout << "                   A match between two partititionsis found by iteratively matching clusters with higher base-node overlap and resolving possible ambiguities." << std::endl;
            std::cout << "                   The matching quality between parittions is defined as the number of base-nodes pairs that are classified in the same way in both partitions" << std::endl;
            std::cout << "                   (both in the smae cluster r both in different clusters) against the total number of pair combinations." << std::endl;
            std::cout << "                   A smart hierarchical search through possible partritions is conducted to find the one with best signature matching." << std::endl << std::endl;
            std::cout << "[-d --depth]:     Partition search depth (for signature and overlap matching. A higher value will mean a more exhaustive search of the possible partitions," << std::endl;
            std::cout << "                   but also a higher computation time, specially if the partition to be matched has a high number of clusters (>100)." << std::endl;
            std::cout << "                   The default value (0, recommended) will adaptively give high search depth to low-cluster partitions and lower search depth to high-cluster partittions." << std::endl << std::endl;
            std::cout << "[-c --justcolor]: Perform only color matching across reference and target tree parttitions (both trees need to have the same number of precompouted partitions)." << std::endl;
            std::cout << "                   In multiple-to-one matching cases clusters from the reference tree might also be recolored to better identify matching relationships across partitions." << std::endl << std::endl;
            std::cout << "[-x --excl]:      Color exclusively clusters that have a match, clusters without match will be recolored white (on both reference and target trees)" << std::endl << std::endl;
            std::cout << "[-v --verbose]:   Verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]:        Write output tree files in vista coordinates (default is nifti)." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "example:" << std::endl << std::endl;
            std::cout << "matchpartition -r refTree.txt -t targetTree.txt -m matching.txt -O results/ -s 0.5 -v" << std::endl << std::endl;
            exit(0);
        }
        if (variableMap.count( "version" ) )
        {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }
        if ( variableMap.count( "verbose" ) )
        {
            std::cout << "verbose output" << std::endl;
            verbose=true;
        }

        if ( variableMap.count( "vista" ) )
        {
            if( verbose )
            {
                std::cout << "Using vista coordinates" << std::endl;
            }
            fileManagerFactory fmf;
            fmf.setVista();
            niftiMode = false;
        }
        else
        {
            if( verbose )
            {
                std::cout << "Using nifti coordinates" << std::endl;
            }
            fileManagerFactory fmf;
            fmf.setNifti();
            niftiMode = true;
        }


        if (variableMap.count("reference"))
        {
            if(!boost::filesystem::is_regular_file( boost::filesystem::path( refTreeFilename ) ) )
            {
                std::cerr << "ERROR: reference tree file \""<<refTreeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Reference tree file: "<< refTreeFilename << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no reference tree file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count( "target" ) )
         {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path( targetTreeFilename ) ) )
            {
                std::cerr << "ERROR: target tree file \""<<targetTreeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Target tree file: "<< targetTreeFilename << std::endl;
        } else {
            std::cerr << "ERROR: no target tree file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count( "leafmatch" ) )
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path( matchTableFilename ) ) )
            {
                std::cerr << "ERROR: match table file \""<<matchTableFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Match table file: "<< matchTableFilename << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no match Table file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if( variableMap.count( "signature" ) + variableMap.count( "overlap" ) + variableMap.count( "justcolor" ) > 1 )
        {
            std::cerr << "ERROR: multiple matching types selected, please use only one from  -s, -o, -c."<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count( "signature" ) )
        {
            std::cout << "Performing Signature partition matching (and color matching)" << std::endl;
            std::cout <<" Using a lambda factor of "<< lambda << std::endl;
            signaturePart = true;
            colorMatching = true;
        }
        else if (variableMap.count( "overlap" ) )
        {
            std::cout << "Performing Overlap partition matching (and color matching): " << std::endl;
            overlapPart = true;
            colorMatching = true;
        }
        else if (variableMap.count( "justcolor" ) )
        {
            std::cout << "Performing only color matching: " << std::endl;
            colorMatching = true;
        }
        else
        {
            std::cerr << "ERROR: no matching type selected, select signature, overlap or color matching"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count( "excl" ) )
        {
            std::cout << "Color exclusively matched clusters (unmatched clusters will be white) " << std::endl;
            exclusive = true;
        }

        if( signaturePart || overlapPart )
        {
            if( searchDepth > 5 )
            {
                std::cout << "Level depth indicated: " << searchDepth << " is too high, setting to a maximum of 5" << std::endl;
                searchDepth = 5;
            }
            else if ( searchDepth = 0 )
            {
                std::cout << "Using automatic parttion-size based adaptive search depth " << std::endl;
            }
            else
            {
                std::cout << "Using a search depth of: " << searchDepth << std::endl;
            }
        }


        if (variableMap.count( "outputf" ) )
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path( outputFolder ) ) )
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
        if(!logFile) {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }
        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Reference tree:\t"<< refTreeFilename <<std::endl;
        logFile <<"Target tree:\t"<< targetTreeFilename <<std::endl;
        logFile <<"Matching table:\t"<< matchTableFilename <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        if( niftiMode )
        {
            logFile << "Using nifti coordinates" << std::endl;
        }
        else
        {
            logFile << "Using vista coordinates" << std::endl;
        }
        logFile <<"-------------"<<std::endl;


        /////////////////////////////////////////////////////////////////


        WHtree refTree( refTreeFilename );
        WHtree targetTree( targetTreeFilename );


        if (!refTree.isLoaded() || !targetTree.isLoaded() )
        {
            throw std::runtime_error ("ERROR @ compareTrees(): trees are not loaded");
        }

        logFile <<"Reference Tree: "<< refTree.getReport(false) <<std::endl;
        logFile <<"Target Tree: "<< targetTree.getReport(false) <<std::endl;

        if (refTree.getDataSize() != targetTree.getDataSize() )
        {
            std::cerr <<"Reference Tree: "<< refTree.getReport() <<std::endl;
            std::cerr <<"Target Tree: "<< targetTree.getReport() <<std::endl;
            throw std::runtime_error ("ERROR @ compareTrees() datasets have different dimensions");
        }

        if (verbose) {
            std::cout <<"Reference Tree: "<< refTree.getReport(false) <<std::endl;
            std::cout <<"Target Tree: "<< targetTree.getReport(false) <<std::endl;
        }


        partitionMatcher matcher(&refTree,&targetTree,matchTableFilename, verbose);
        std::string depthString;
        if( searchDepth > 0 )
        {
            depthString = "_d"+string_utils::toString< unsigned int >( searchDepth ) ;
        }

        std::cout <<matcher.reportBaseNodes()<<std::endl;
        std::string suffixPart("_pm_Signature_l" + string_utils::toString< float >( lambda ) + depthString + ".txt");
        std::string suffixNew("_pm_Overlap" + depthString + ".txt");
        std::string suffixColor("_colorMatch.txt");
        bool refTreeColorsChanged( false );


        if( signaturePart )
        {
            logFile <<  "Signature Matching" << std::endl;
            logFile <<  "Lambda:\t" << lambda <<std::endl;
            logFile <<  "Search depth:\t" << searchDepth <<std::endl;
            matcher.findMatchingPartitions( lambda );
        }
        else if ( overlapPart )
        {
            logFile <<  "Overlap Matching" <<std::endl;
            logFile <<  "Search depth:\t" << searchDepth <<std::endl;
            matcher.findMatchingPartitions( -1 );
        }

        if( colorMatching )
        {
            logFile <<  "Color Matching" <<std::endl;
            refTreeColorsChanged = matcher.matchColors( exclusive );

        }

        std::string refOutput;
        std::string targetOutput;

        if ( signaturePart )
        {
            targetOutput = outputFolder + "/" + targetTree.getName() + suffixPart;
        }
        else if ( overlapPart )
        {
            targetOutput = outputFolder + "/" + targetTree.getName() + suffixNew ;
        }
        else
        {
            targetOutput = outputFolder + "/" + targetTree.getName() + suffixColor;
        }

        if( refTreeColorsChanged )
        {
            refOutput = outputFolder + "/" + refTree.getName() + suffixColor;
        }
        else
        {
            refOutput = outputFolder + "/" + refTree.getName() + ".txt";
        }

        if( verbose )
        {
            std::cout <<  "Writing output target tree file to " << targetOutput << std::endl;
            std::cout <<  "Writing output reference tree file to " << refOutput << std::endl;
        }

        targetTree.writeTree( targetOutput, niftiMode );
        refTree.writeTree( refOutput, niftiMode );

        logFile <<  "Written output target tree file to " << targetOutput << std::endl;
        logFile <<  "Written output reference tree file to " << refOutput << std::endl;

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

