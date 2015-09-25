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
//  processtree
//
//  Do tree processing, either full raw tree preprocessing (monotonicity correction, base-node flattening, debinarization...) or/and linear node collapse.
//   For more information on the preprocessing steps refer to (Moreno-Dominguez, 2014).
//   For an interactive tree processing management with more options please use the Hierarchcial Clustering module developed in OpenWalnut (www.openwalnut.org).
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -t --tree:       File with the hierarchical tree to preprocess.
//
//   -O --outputf:    Output folder where processed tree files will be written.
//
//  [-n --name]:      Prefix for the output tree filename.
//
//  [-c --collapse]:  Perform linear node collapse in order to de-binarize non-binary structures.
//                     Recommended collapse factorvalue: 0.05.
//
//  [--ignorebases]:  (use only alongiside option -c) ignore node-base status when performing node collapse.
//
//  [-r --raw]:       Do full processing of binary raw input tree.
//
//  [-b --bases]:     (use only alongiside option -r) do base-nodes (meta-leaves) flattening.
//                     Requires file with base-nodes indentifiers.
//
//  [-m --monmult]:   Monotonicity error multiplier. Increase if monotonicity correction enters an infinite loop.
//                     Default value: 1 (no multiplier).
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        write output tree in vista coordinates (default is nifti).
//
//  [--debugout]:     write additional detailed outputs meant to be used for debugging.
//
//
//  * Usage example:
//
//   processtree -t tree_lh.txt -O results/ -n processedtree -raw -c -v
//
//
//  * Outputs (in output folder defined at option -O):
//
//   - The processed tree file with either the same original name as the one defined by option -t, or the name defined by option -n when used.
//   - If both option -r and -c are used, the previous statement refers to the file with the processed raw-tree and the furthermore collapes output will be written with the "_collapsed" suffix.
//   - 'processtree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program.
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
#include "WHtreeProcesser.h"
#include "fileManagerFactory.h"

size_t numComps(0);



int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("processtree");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string treeFilename, outputFolder, treeName, basesFilename;
        float collapseFactor( 0.5 ), errorMult( 1 );
        bool raw( false ), collapse( false ), bases( false );
        bool ignoreBases( false ), verbose( false ), debug( false ), niftiMode( true );


        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "tree,t",  boost::program_options::value< std::string >(&treeFilename), "path to the tree file")
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder where processed tree(s) will be written")
                ( "name,n",  boost::program_options::value< std::string >(&treeName), "[opt] name prefix for output tree")
                ( "collapse,c", boost::program_options::value< float >(&collapseFactor), "[opt] perform linear collapse, enter collapse factor (recommended value: 0.05)")
                ( "ignorebases", "[use only with -c] allow level 1 nodes (with leaf children) to be eliminated in the collapse process (use only on trees without defined base-nodes)")
                ( "raw,r", "[opt] full tree processing from raw binary input tree")
                ( "bases,b",  boost::program_options::value< std::string >(&basesFilename), "[use only with -r] do base-nodes (metaleaves) flattening.")
                ( "monmult,m",boost::program_options::value< float >(&errorMult),  "[use only with -r] monotonicity error multiplier.Default: 1 (no multiplier)")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "debugout", "[opt] write additional detailed outputs meant for debug." )
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

        if (variableMap.count("verbose"))
        {
            std::cout << "verbose output"<<std::endl;
            verbose=true;
        }


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
            std::cout << "processtree" << std::endl << std::endl;
            std::cout << "Do tree processing, either full raw tree preprocessing (monotonicity correction, base-node flattening, debinarization...) or/and linear node collapse." << std::endl;
            std::cout << " For more information on the preprocessing steps refer to (Moreno-Dominguez, 2014)." << std::endl;
            std::cout << " For an interactive tree processing management with more options please use the Hierarchcial Clustering module developed in OpenWalnut (www.openwalnut.org)." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -t --tree:       File with the hierarchical tree to preprocess." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where processed tree files will be written." << std::endl << std::endl;
            std::cout << "[-n --name]:      Prefix for the output tree filename." << std::endl << std::endl;
            std::cout << "[-c --collapse]:  Perform linear node collapse in order to de-binarize non-binary structures." << std::endl;
            std::cout << "                   Recommended collapse factorvalue: 0.05." << std::endl << std::endl;
            std::cout << "[--ignorebases]:  (use only alongiside option -c) ignore node-base status when performing node collapse." << std::endl << std::endl;
            std::cout << "[-r --raw]:       Do full processing of binary raw input tree." << std::endl << std::endl;
            std::cout << "[-b --bases]:     (use only alongiside option -r) do base-nodes (meta-leaves) flattening." << std::endl;
            std::cout << "                   Requires file with base-nodes indentifiers." << std::endl << std::endl;
            std::cout << "[-m --monmult]:   Monotonicity error multiplier. Increase if monotonicity correction enters an infinite loop." << std::endl;
            std::cout << "                   Default value: 1 (no multiplier)." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	    write output tree in vista coordinates (default is nifti)." << std::endl << std::endl;
            std::cout << "[--debugout]:     write additional detailed outputs meant to be used for debugging." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " processtree -t tree_lh.txt -O results/ -n processedtree -raw -c -v" << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " - The processed tree file with either the same original name as the one defined by option -t, or the name defined by option -n when used." << std::endl;
            std::cout << " - If both option -r and -c are used, the previous statement refers to the file with the processed raw-tree and the furthermore collapes output will be written with the '_collapsed'' suffix." << std::endl;
            std::cout << " - 'processtree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
            std::cout << std::endl;
            exit(0);
        }
        if (variableMap.count("version"))
        {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }

        if ( variableMap.count( "debugout" ) )
        {
            if( verbose )
            {
                std::cout << "Debug output files activated" << std::endl;
            }
            debug = true;
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

        if (variableMap.count("tree"))
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename)))
            {
                std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Tree file: "<< treeFilename << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no tree file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("name"))
        {
            if( verbose )
            {
                std::cout << "Output tree name prefix: "<< treeName << std::endl;
            }
        }
        else
        {
            treeName = boost::filesystem::path(treeFilename).stem().string();
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

        if (variableMap.count("collapse"))
        {
            if( verbose )
            {
                std::cout << "Linear collapse of nodes with collapse factor: " << collapseFactor << std::endl;
            }
            collapse = true;
            if (variableMap.count("ignorebases"))
            {
                if( verbose )
                {
                    std::cout << "Ignoring base-node status during node collapse"<<std::endl;
                }
                ignoreBases = true;
            }
        }


        if (variableMap.count("raw"))
        {
            if( verbose )
            {
                std::cout << "Full tree processing (from binary raw tree)"<<std::endl;
            }
            raw = true;


            if (variableMap.count("bases"))
            {
                if(!boost::filesystem::is_regular_file(boost::filesystem::path(basesFilename)))
                {
                    std::cerr << "ERROR: bases file \""<<basesFilename<<"\" is not a regular file"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }
                if( verbose )
                {
                    std::cout << "Including base-node processing. Bases file: "<< basesFilename <<std::endl;
                }
                bases = true;
            }

            if (variableMap.count("monmult"))
            {
                if(errorMult < 1 )
                {
                    std::cerr << "WARNING: Invalid monotonicity error multipier ("<< errorMult << ") values must be >= 1. Setting it to 1 (no multiplier)" << std::endl;
                    errorMult =1;
                }
                else if( verbose )
                {
                    std::cout << "monotonicity error multipier: "<< errorMult << std::endl;
                }
            }
        }

        /////////////////////////////////////////////////////////////////

        std::string logFilename(outputFolder+"/"+progName+"_log.txt");
        std::ofstream logFile(logFilename.c_str());
        if(!logFile)
        {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }
        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Tree file:\t"<< treeFilename <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Output tree prefix:\t"<< treeName <<std::endl;
        logFile <<"ODebug outputs flag:\t"<< debug <<std::endl;


        WHtree tree(treeFilename);
        WHtreeProcesser treePrcssr(&tree);
        std::vector< size_t > baseVector;
        std::vector< size_t > prunedVector;

        logFile << tree.getReport( false ) <<std::endl;
        std::cout<<tree.getReport( false )<<std::endl;



        if( raw )
        {
            if( bases )
            {
                WFileParser parser( basesFilename );

                if( !parser.readFile() )
                {
                    std::cerr << "ERROR: Parser error when reading bases" << std::endl;
                    exit(-1);
                }
                std::vector<std::string> lines = parser.getRawLines();
                if( lines.size() == 0 )
                {
                    std::cerr << "ERROR: bases file is empty" << std::endl;
                    exit(-1);
                }
                {
                    std::vector< std::vector< std::string> >baseStrings = parser.getLinesForTagSeparated( "bases" );
                    baseVector.reserve( baseStrings.size() );
                    for( size_t i = 0; i < baseStrings.size(); ++i )
                    {
                        if( baseStrings[i].size() != 1 )
                        {
                            std::cerr << "ERROR: multiple base IDs in the same line, check format" << std::endl;
                            exit(-1);
                        }
                        baseVector.push_back( boost::lexical_cast<size_t>( baseStrings[i][0] ) );
                    }

                    std::vector< std::vector< std::string> >prunedStrings = parser.getLinesForTagSeparated( "pruned" );
                    prunedVector.reserve( prunedStrings.size() );
                    for( size_t i = 0; i < prunedStrings.size(); ++i )
                    {
                        if( prunedStrings[i].size() != 1 )
                        {
                            std::cerr << "ERROR: multiple pruned IDs in the same line, check format" << std::endl;
                            exit(-1);
                        }
                        prunedVector.push_back( boost::lexical_cast<size_t>( prunedStrings[i][0] ) );
                    }
                }
            }
            treePrcssr.flagLeaves( prunedVector );

            WHtree treeUp(tree);
            WHtreeProcesser treePrcssrUp(&treeUp);
            WHtree treeDown(tree);
            WHtreeProcesser treePrcssrDown(&treeDown);

            if( verbose )
            {
                std::cout<< "Starting full tree preprocessing..." << std::endl;
            }


            if( verbose )
            {
                std::cout<< "forcing monotonicity." << std::endl;
            }
            logFile<< "forcing monotonicity." << std::endl;
            treePrcssr.forceMonotonicity(errorMult);
            logFile << tree.getReport( false ) <<std::endl;
            if( verbose )
            {
                std::cout<<tree.getReport( false )<<std::endl;
            }

            if( debug )
            {
                std::string outFilename = outputFolder + "/" + treeName + "_bin.txt";
                tree.writeTree( outFilename, niftiMode );
                if( verbose )
                {
                    std::cout << "written to: "<< outFilename <<std::endl;
                }
                logFile << "written to: "<< outFilename <<std::endl;

                if( verbose )
                {
                    std::cout<< "forcing monotonicity up." << std::endl;
                }
                logFile<< "forcing monotonicity up." << std::endl;
                treePrcssrUp.forceMonotonicityUp();
                if( verbose )
                {
                    std::cout<<treeUp.getReport( false )<<std::endl;
                }
                logFile << treeUp.getReport( false ) <<std::endl;
                outFilename = outputFolder + "/" + treeName + "_bin_Up.txt";
                treeUp.writeTree( outFilename, niftiMode );
                if( verbose )
                {
                    std::cout << "written to: "<< outFilename <<std::endl;
                }
                logFile << "written to: "<< outFilename <<std::endl;

                if( verbose )
                {
                    std::cout<< "forcing monotonicity down." << std::endl;
                }
                logFile<< "forcing monotonicity down." << std::endl;
                treePrcssrDown.forceMonotonicityDown();
                if( verbose )
                {
                    std::cout<<treeDown.getReport( false )<<std::endl;
                }
                logFile << treeDown.getReport( false ) <<std::endl;
                outFilename = outputFolder + "/" + treeName + "_bin_Down.txt";
                treeDown.writeTree( outFilename, niftiMode );
                if( verbose )
                {
                    std::cout << "written to: "<< outFilename <<std::endl;
                }
                logFile << "written to: "<< outFilename <<std::endl;
            }

            if( bases )
            {
                if( verbose )
                {
                    std::cout<< "Flattening base nodes and pruning out unconnected voxels." << std::endl;
                }
                treePrcssr.flattenSelection( baseVector, false );
                logFile << tree.getReport( false ) <<std::endl;
                if( verbose )
                {
                    std::cout << tree.getReport( false )<<std::endl;
                }

                if( debug )
                {
                    if( verbose )
                    {
                        std::cout<< "Flattening base nodes and pruning out unconnected voxels (UP-tree)." << std::endl;
                    }
                    treePrcssrUp.flattenSelection( baseVector, false );
                    logFile << treeUp.getReport( false ) <<std::endl;
                    if( verbose )
                    {
                        std::cout << treeUp.getReport( false )<<std::endl;
                    }

                    if( verbose )
                    {
                        std::cout<< "Flattening base nodes and pruning out unconnected voxels (DOWN-tree)." << std::endl;
                    }
                    treePrcssrDown.flattenSelection( baseVector, false );
                    logFile << treeDown.getReport( false ) <<std::endl;
                    if( verbose )
                    {
                        std::cout << treeDown.getReport( false )<<std::endl;
                    }
                }
            }

            if( verbose )
            {
                std::cout<< "Debinarizing." << std::endl;
            }
            treePrcssr.debinarize( false );
            logFile << tree.getReport( false ) <<std::endl;
            if( verbose )
            {
                std::cout << tree.getReport( false )<<std::endl;
            }

            std::string outFilename = outputFolder + "/" + treeName + ".txt";
            tree.writeTree( outFilename, niftiMode );
            if( verbose )
            {
                std::cout << "written to: "<< outFilename <<std::endl;
            }
            logFile << "written to: "<< outFilename <<std::endl;
            if( tree.testRootBaseNodes() )
            {
                std::vector<size_t> baseVector = tree.getRootBaseNodes();
                std::sort( baseVector.begin(), baseVector.end() );
                treePrcssr.writeBases( baseVector, outputFolder + "/baselist.txt" );

                if( verbose )
                {
                    std::cout << "Final base list written in: "<< outputFolder << "/baselist.txt" << std::endl;
                }
                logFile << "Final base list written in: "<< outputFolder << "/baselist.txt" << std::endl;
            }
            else
            {
                if( verbose )
                {
                    std::cout << "Final tree is not a pure basenode tree" << std::endl;
                }
                logFile << "Final tree is not a pure basenode tree" << std::endl;
            }

            if( debug )
            {

                if( verbose )
                {
                    std::cout<< "Debinarizing (UP-tree)." << std::endl;
                }
                treePrcssrUp.debinarize( false );
                logFile << treeUp.getReport( false ) <<std::endl;
                if( verbose )
                {
                    std::cout << treeUp.getReport( false )<<std::endl;
                }

                outFilename = outputFolder + "/" + treeName + "_Up.txt";
                treeUp.writeTree( outFilename, niftiMode );
                if( verbose )
                {
                    std::cout << "written to: "<< outFilename <<std::endl;
                }
                logFile << "written to: "<< outFilename <<std::endl;
                if( treeUp.testRootBaseNodes() )
                {
                    std::vector<size_t> baseVector = treeUp.getRootBaseNodes();
                    std::sort( baseVector.begin(), baseVector.end() );
                    treePrcssrUp.writeBases( baseVector, outputFolder + "/baselist_Up.txt" );

                    if( verbose )
                    {
                        std::cout << "Final base list written in: "<< outputFolder << "/baselist_Up.txt" << std::endl;
                    }
                    logFile << "Final base list written in: "<< outputFolder << "/baselist_up.txt" << std::endl;
                }
                else
                {
                    if( verbose )
                    {
                        std::cout << "Final UP-tree is not a pure basenode tree" << std::endl;
                    }
                    logFile << "Final UP-tree is not a pure basenode tree" << std::endl;
                }


                if( verbose )
                {
                    std::cout<< "Debinarizing (DOWN-tree)." << std::endl;
                }
                treePrcssrDown.debinarize( false );
                logFile << treeDown.getReport( false ) <<std::endl;
                if( verbose )
                {
                    std::cout << treeDown.getReport( false )<<std::endl;
                }

                outFilename = outputFolder + "/" + treeName + "_Down.txt";
                treeDown.writeTree( outFilename, niftiMode );
                if( verbose )
                {
                    std::cout << "written to: "<< outFilename <<std::endl;
                }
                logFile << "written to: "<< outFilename <<std::endl;
                if( treeDown.testRootBaseNodes() )
                {
                    std::vector<size_t> baseVector = treeDown.getRootBaseNodes();
                    std::sort( baseVector.begin(), baseVector.end() );
                    treePrcssrDown.writeBases( baseVector, outputFolder + "/baselist_Down.txt" );

                    if( verbose )
                    {
                        std::cout << "Final base list written in: "<< outputFolder << "/baselist_Down.txt" << std::endl;
                    }
                    logFile << "Final base list written in: "<< outputFolder << "/baselist_Down.txt" << std::endl;
                }
                else
                {
                    if( verbose )
                    {
                        std::cout << "Final DOWN-tree is not a pure basenode tree" << std::endl;
                    }
                    logFile << "Final tree is not a pure basenode tree" << std::endl;
                }
            }

        } // end raw

        if( collapse )
        {
            if( verbose )
            {
                std::cout<< "Performing linear node collapse, collapse factor: "<< collapseFactor << std::endl;
            }
            logFile<< "Performing linear node collapse, collapse factor: "<< collapseFactor << std::endl;
            treePrcssr.collapseTreeLinear(collapseFactor , !ignoreBases);
            logFile << tree.getReport( false ) <<std::endl;
            if( verbose )
            {
                std::cout << tree.getReport( false )<<std::endl;
            }
            std::string outFilename;
            if( raw )
            {
                outFilename = outputFolder + "/" + treeName + "_collapsed.txt";
            }
            else
            {
                outFilename = outputFolder + "/" + treeName + ".txt";
            }
            tree.writeTree( outFilename, niftiMode );
            if( verbose )
            {
            std::cout << "written to: "<< outFilename <<std::endl;
            }
            logFile << "written to: "<< outFilename <<std::endl;

            if( raw && debug )
            {
                WHtree treeUp(  outputFolder + "/" + treeName + "_Up.txt" );
                WHtree treeDown(  outputFolder + "/" + treeName + "_Down.txt" );
                WHtreeProcesser treeProcesserUp( &treeUp );
                WHtreeProcesser treeProcesserDown( &treeDown );

                if( verbose )
                {
                    std::cout<< "Performing linear node collapse on UP-tree "<< std::endl;
                }
                logFile<< "Performing linear node collapse on UP-tree "<< std::endl;
                treeProcesserUp.collapseTreeLinear(collapseFactor , !ignoreBases);
                logFile << treeUp.getReport( false ) <<std::endl;
                if( verbose )
                {
                    std::cout << treeUp.getReport( false )<<std::endl;
                }
                outFilename = outputFolder + "/" + treeName + "_collapsed_Up.txt";
                treeUp.writeTree( outFilename, niftiMode );
                if( verbose )
                {
                    std::cout << "written to: "<< outFilename <<std::endl;
                }
                logFile << "written to: "<< outFilename <<std::endl;

                if( verbose )
                {
                    std::cout<< "Performing linear node collapse on DOWN-tree "<< std::endl;
                }
                logFile<< "Performing linear node collapse on DOWN-tree "<< std::endl;
                treeProcesserDown.collapseTreeLinear(collapseFactor , !ignoreBases);
                logFile << treeDown.getReport( false ) <<std::endl;
                if( verbose )
                {
                    std::cout << treeDown.getReport( false )<<std::endl;
                }
                outFilename = outputFolder + "/" + treeName + "_collapsed_Down.txt";
                treeDown.writeTree( outFilename, niftiMode );
                if( verbose )
                {
                    std::cout << "written to: "<< outFilename <<std::endl;
                }
                logFile << "written to: "<< outFilename <<std::endl;
            }
        }// end collapse


        /////////////////////////////////////////////////////////////////

        // save and print total time
        time_t programEndTime(time(NULL));
        int totalTime( difftime(programEndTime,programStartTime) );
        std::cout <<"Program Finished, total time: "<< totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\"   "<< std::endl;
        logFile <<"-------------"<<std::endl;
        logFile <<"Finish Time:\t"<< ctime(&programEndTime) <<std::endl;
        logFile <<"Elapsed time : "<< totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\""<< std::endl;

        std::cout <<"Total correlations: "<< numComps << std::endl;
        logFile <<"Total correlations: "<< numComps << std::endl;


//    }
//    catch(std::exception& e)
//    {
//        std::cout << e.what() << std::endl;
//        return 1;
//    }
    return 0;
}


#if 0

// Declare a group of options that will be allowed only on command line
boost::program_options::options_description genericOptions("Generic options");
genericOptions.add_options()
        ("version,V", "print version string")
        ("help,h", "produce help message")
        ("monotonic,m", "force tree monotonicity")
        ("tree-file,t",  boost::program_options::value< std::string >(&treeFilename), "file with the tree file")
        ("prune-size,s", boost::program_options::value< size_t >(&joinSize), "prune branches smaller than safe-size that join a node equal or bigger than N")
        ("prune-ratio,r", boost::program_options::value< float >(&sizePruneRatio), "prune branches smaller than safe-size whose size ratio to its parent node is smaller than N")
        ("prune-level,l", boost::program_options::value< float >(&distPruneLevel), "prune branches smaller than safe-size whose parent its at a distance level grater than N")
        ("prune-random,r", boost::program_options::value< size_t >(&randPruneNumber), "prune N randomly selected leaves")
        ("flatten,f", boost::program_options::value< float >(&flatGap), "flatten all those nodes whose distance to their respective parents/children are smaller than N")
        ("coarse,c", boost::program_options::value< unsigned int >(&coarseRatio), "decimate tree results to fit an N-times coarser anatomical grid")
        ("smoothen,h", "performs tree smoothing, requires prune-size, prune-level and safe-size parameters, other options will be ignored")
        ("base-tree,b", "converts base nodes into leaves and writes the simplififed tree")
        ("verbose,v", "verbose option")
        ;


if (variableMap.count("help")) {
    std::cout << visibleOptions << std::endl;
    exit(0);
}
if (variableMap.count("version")) {
    std::cout << progName <<", version 1.0"<<std::endl;
    exit(0);
}
if (variableMap.count("verbose")) {
    std::cout << "verbose output"<<std::endl;
    verbose=true;
}

if (variableMap.count("tree-file")) {
    if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename))) {
        std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
        std::cerr << visibleOptions << std::endl;
        exit(-1);
    }
    std::cout << "Roi voxels file: "<< treeFilename << std::endl;
} else {
    std::cerr << "ERROR: no tree file stated"<<std::endl;
    std::cerr << visibleOptions << std::endl;
    exit(-1);
}


if (variableMap.count("output")) {
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

/////////////////////////////////////////////////////////////////


std::string outFile;

std::string logFilename(outputFolder+"/"+progName+"_log.txt");
std::ofstream logFile(logFilename.c_str());
if(!logFile) {
    std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
    exit(-1);
}
logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
logFile <<"Verbose:\t"<< verbose <<std::endl;
logFile <<"Tree file:\t"<< treeFilename <<std::endl;
logFile <<"Output folder:\t"<< outputFolder <<std::endl;

WHtree tree(treeFilename);

logFile << tree.getReport( false ) <<std::endl;
std::cout<<tree.getReport( false )<<std::endl;


WHtreeProcesser treePrcssr(&tree);

if (variableMap.count("monotonic") )
{
    WHtree treeUp(tree);
    WHtreeProcesser treePrcssrUp(&treeUp);
    std::cout<< "forcing monotonicity up." << std::endl;
    logFile<< "forcing monotonicity up." << std::endl;
    treePrcssrUp.forceMonotonicityUp();
    logFile << treeUp.getReport( false ) <<std::endl;
    std::cout<<treeUp.getReport( false )<<std::endl;
    outFile = outputFolder + "/" + treeUp.getName() + "_monoUp.txt";
    treeUp.writeTree( outFilename, niftiMode );
    std::cout << "written to: "<< outFile <<std::endl;
    logFile << "written to: "<< outFile <<std::endl;

    WHtree treeDown(tree);
    WHtreeProcesser treePrcssrDown(&treeDown);
    std::cout<< "forcing monotonicity down." << std::endl;
    logFile<< "forcing monotonicity down." << std::endl;
    treePrcssrDown.forceMonotonicityDown();
    logFile << treeDown.getReport( false ) <<std::endl;
    std::cout<<treeDown.getReport( false )<<std::endl;
    outFile = outputFolder + "/" + treeDown.getName() + "_monoDown.txt";
    treeDown.writeTree( outFilename, niftiMode );
    std::cout << "written to: "<< outFile <<std::endl;
    logFile << "written to: "<< outFile <<std::endl;

    std::cout<< "forcing monotonicity weighted." << std::endl;
    logFile<< "forcing monotonicity weighted." << std::endl;
    treePrcssr.forceMonotonicity(errorMult);
    logFile << tree.getReport( false ) <<std::endl;
    std::cout<<tree.getReport( false )<<std::endl;
    outFile = outputFolder + "/" + tree.getName() + "_monoW.txt";
    tree.writeTree( outFilename, niftiMode );
    std::cout << "written to: "<< outFile <<std::endl;
    logFile << "written to: "<< outFile <<std::endl;
}


if (variableMap.count("smoothen") )
{
    if (variableMap.count("safe-size") && variableMap.count("prune-level") && variableMap.count("prune-size")) {
        std::cout << "Smoothing... safe-size: "<<safeSize<<". smooth-size: "<<joinSize<<". smooth-gap: "<<distPruneLevel<<std::endl;
        logFile <<"Smoothing => safe-size: "<<safeSize<<". smooth-size: "<<joinSize<<". smooth-gap: "<<distPruneLevel<<std::endl;

        if (variableMap.count("prune-ratio") || variableMap.count("prune-random") || variableMap.count("flatten")) {
            std::cerr << "WARNING: prune-ratio, prune-random and flatten options will be ignored"<<std::endl;
        }

        treePrcssr.smoothTree(safeSize,joinSize,distPruneLevel);

        std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
        logFile <<tree.getReport(false) <<std::endl;
        outFile = outputFolder + "/" + tree.getName() + "_smooth.txt";
        tree.writeTree( outFilename, niftiMode );
        outFile = outputFolder + "/" + tree.getName() + "_smooth_debug.txt";
        tree.writeTreeDebug(outFile);


    }
    else
    {
        std::cerr << "ERROR: smoothing requires the additional parameters: safe-size, prune-level and prune-size"<<std::endl;
        exit(-1);
    }
}
else
{


    if (variableMap.count("prune-ratio") || variableMap.count("prune-level") || variableMap.count("prune-size"))
    {
        if (variableMap.count("prune-random"))
        {
            std::cerr << "ERROR: pruned random may not be used simultaneously with other pruning methods"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        std::string safeString;
        if(variableMap.count("safe-size"))
        {
            safeString = ("smaller than "+boost::lexical_cast<std::string>(safeSize));
            logFile <<"Pruning-safe size:\t"<< safeSize <<std::endl;

            if (variableMap.count("prune-size"))
            {
                std::cout << "Pruning by minimum size if a cluster smaller than "<< safeSize <<" joins a node bigger than more than "<< joinSize <<". it will be pruned"<< std::endl;
                logFile <<"Min-size pruning:\t"<< safeSize << ":"<< joinSize <<std::endl;

                treePrcssr.pruneTree(joinSize,safeSize,HTPR_JOINSIZE);

                std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
                logFile <<tree.getReport(false) <<std::endl;
                outFile = outputFolder + "/" + tree.getName() + "_prunedS.txt";
                tree.writeTree( outFilename, niftiMode );
                outFile = outputFolder + "/" + tree.getName() + "_prunedS_debug.txt";
                tree.writeTreeDebug(outFile);

            }
        }


        if (variableMap.count("prune-ratio"))
        {
            std::cout << "Pruning by size ratio: if a cluster "<< safeString <<" joins a node more than "<< sizePruneRatio <<" times its size, it will be pruned"<< std::endl;
            logFile <<"Size-pruning ratio:\t"<< sizePruneRatio <<std::endl;

            treePrcssr.pruneTree(sizePruneRatio,safeSize,HTPR_SIZERATIO);

            std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
            logFile <<tree.getReport(false) <<std::endl;
            outFile = outputFolder + "/" + tree.getName() + "_prunedR.txt";
            tree.writeTree( outFilename, niftiMode );
            outFile = outputFolder + "/" + tree.getName() + "_prunedR_debug.txt";
            tree.writeTreeDebug( outFilename, niftiMode );

        }
        else if (variableMap.count("prune-level"))
        {
            if (safeSize==0)
            {
                std::cerr << "ERROR: when pruning by level a positive safe size must be stated"<<std::endl;
                exit(-1);
            }
            std::cout << "Pruning by distance level: if a cluster "<< safeString <<" joins a node at a level higher than "<< distPruneLevel <<" times its size, it will be pruned"<< std::endl;
            logFile <<"Level-pruning distance:\t"<< distPruneLevel <<std::endl;

            treePrcssr.pruneTree(distPruneLevel,safeSize,HTPR_JOINLEVEL);

            std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
            logFile <<tree.getReport(false) <<std::endl;
            outFile = outputFolder + "/" + tree.getName() + "_prunedL.txt";
            tree.writeTree( outFilename, niftiMode );
            outFile = outputFolder + "/" + tree.getName() + "_prunedL_debug.txt";
            tree.writeTreeDebug( outFilename, niftiMode );
        }

    }

    if (variableMap.count("prune-random"))
    {
        std::cout << "Randomly pruning "<< randPruneNumber <<" leaves from the tree"<< std::endl;
        logFile <<"Randomly pruned leaves:\t"<< randPruneNumber <<std::endl;


        std::cout << tree.getReport(false) <<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        treePrcssr.pruneRandom(randPruneNumber);

        std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_prunedRand.txt";
        tree.writeTree( outFilename, niftiMode );

        std::cout << "written to: "<< outFile <<std::endl;
        logFile << "written to: "<< outFile <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_prunedRand_debug.txt";
        tree.writeTreeDebug( outFilename, niftiMode );




    }

    if (variableMap.count("flatten"))
    {
        std::cout << "Flattening branches with intervals smaller than "<< flatGap << std::endl;
        logFile <<"Flattening gap:\t"<< flatGap <<std::endl;

        std::cout << tree.getReport(false) <<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        treePrcssr.collapseTree(flatGap,1);

        std::cout << tree.getReport(false) <<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_flat.txt";
        tree.writeTree( outFilename, niftiMode );

        std::cout << "written to: "<< outFile <<std::endl;
        logFile << "written to: "<< outFile <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_flat_debug.txt";
        tree.writeTreeDebug( outFilename, niftiMode );


    }

}

if (variableMap.count("coarse"))
{
    std::cout << "Decimating tree to a grid "<< coarseRatio <<" times coarser" << std::endl;
    logFile <<"Coarsing factor:\t"<< coarseRatio <<std::endl;

    treePrcssr.coarseTree(coarseRatio);

    std::cout<<"Done. "<<tree.getReport()<<std::endl;
    logFile <<tree.getReport() <<std::endl;
    outFile = outputFolder + "/" + tree.getName() + "_coarse.txt";
    tree.writeTree( outFilename, niftiMode );
    outFile = outputFolder + "/" + tree.getName() + "_coarse_debug.txt";
    tree.writeTreeDebug(outFile);

}

if (variableMap.count("base-tree"))
{
    logFile <<tree.getReport() <<std::endl;
    std::cout << tree.getReport() <<std::endl;

    std::cout << "Collapsing base nodes into leaves..." << std::endl;
    logFile << "Collapsing base nodes into leaves" << std::endl;
    treePrcssr.baseNodes2Leaves();

    std::cout<<"Done. "<<tree.getReport()<<std::endl;

    logFile <<tree.getReport() <<std::endl;
    std::cout << tree.getReport() <<std::endl;

    outFile = outputFolder + "/" + tree.getName() + "_baset.txt";
    tree.writeTreeSimple(outFile);
}


/////////////////////////////////////////////////////////////////


// save and print total time
time_t programEndTime(time(NULL));
int totalTime( difftime(programEndTime,programStartTime) );
std::cout <<"Program Finished, total time: "<< totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\"   "<< std::endl;
logFile <<"-------------"<<std::endl;
logFile <<"Finish Time:\t"<< ctime(&programEndTime) <<std::endl;
logFile <<"Elapsed time : "<< totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\""<< std::endl;

std::cout <<"Total correlations: "<< numComps << std::endl;
logFile <<"Total correlations: "<< numComps << std::endl;



#endif

