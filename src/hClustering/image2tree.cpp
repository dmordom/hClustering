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
//  image2tree
//
//  Creates a tree with base nodes matching those of an input tree file and structure matching that of an input single partition 3D image.
//   It uses the partition and roi coordinates information in the 3D image to assign each base-node (meta-leaf) of a hierarchical tree
//   to one of the partition clusters, then create a new tree with the same base nodes as the orginal but with only one partition in the hierarchical structure:
//   the most similar to the one defined in the 3D partition label image. This single-partition tree can then be used to perform tree-comparison statistics between the
//   original tree and the single-partition tree.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -t --tree:       File with the hierarchical tree to be used as base node template.
//
//  [-b --bases:]     File with the tree meta-leaves (base-nodes) identifiers. If omitted base nodes will be calculated from the tree.
//
//   -i --image:      File with the 3D partition label image that wishes to be projected into a tree with matching base nodes to the input tree.
//
//   -O --outputf:    Folder where the output partition tree file will be written.
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Write output tree in vista coordinates (default is nifti).
//
//
//  * Usage example:
//
//  image2tree -t tree.txt -i partition.nii -O results/ -v
//
//
//  * Outputs (in output folder defined at option -O):
//
//   - 'partitionTree.txt' - A copy of the original tree file with the best-matched partitions to the 3D label file included in the relevant fields.
//   - 'success.txt' - An empty file created when the program has sucessfully exited after completion (to help for automatic re-running scripting after failure).
//   - 'image2tree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program.
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
#include "common/image2treeBuilder.h"

size_t numComps(0);



int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("image2tree");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string treeFilename, baseFilename, imageFilename, outputFolder;
        bool verbose(false), niftiMode( true );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "tree,t",  boost::program_options::value< std::string >(&treeFilename), "file with the tree to be used as base-node template")
                ( "bases,b",  boost::program_options::value< std::string >(&baseFilename), "[opt] file with the tree base-nodes identifiers")
                ( "image,i",  boost::program_options::value< std::string >(&imageFilename), "3D partition image file")
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "Output folder where partition files will be written")
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
            std::cout << "image2tree" << std::endl << std::endl;
            std::cout << "Creates a tree with base nodes matching those of an input tree file and structure matching that of an input single partition 3D image." << std::endl;
            std::cout << " It uses the partition and roi coordinates information in the 3D image to assign each base-node (meta-leaf) of a hierarchical tree" << std::endl;
            std::cout << " to one of the partition clusters, then create a new tree with the same base nodes as the orginal but with only one partition in the hierarchical structure:" << std::endl;
            std::cout << " the most similar to the one defined in the 3D partition label image. This single-partition tree can then be used to perform tree-comparison statistics between the" << std::endl;
            std::cout << " original tree and the single-partition tree." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -t --tree:       File with the hierarchical tree to be used as base node template." << std::endl << std::endl;
            std::cout << "[-b --bases:]     File with the tree meta-leaves (base-nodes) identifiers. If omitted base nodes will be calculated from the tree." << std::endl << std::endl;
            std::cout << " -i --image:      File with the 3D partition label image that wishes to be projected into a tree with matching base nodes to the input tree." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where partition files will be written." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	    write output tree in vista coordinates (default is nifti)." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " image2tree -t tree.txt -i partition.nii -O results/ -v" << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " - 'partitionTree.txt' - A copy of the original tree file with the best-matched partitions to the 3D label file included in the relevant fields." << std::endl;
            std::cout << " - 'success.txt' - An empty file created when the program has sucessfully exited after completion (to help for automatic re-running scripting after failure)." << std::endl;
            std::cout << " - 'image2tree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
            std::cout << std::endl;
            exit(0);
        }

        if (variableMap.count("version"))
        {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }
        if (variableMap.count("verbose"))
        {
            std::cout << "verbose output"<<std::endl;
            verbose=true;
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


        if ( variableMap.count( "tree" ) )
        {
            if( !boost::filesystem::is_regular_file( boost::filesystem::path( treeFilename ) ) )
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

        if (  variableMap.count( "bases" ) )
        {
            if( !boost::filesystem::is_regular_file( boost::filesystem::path( baseFilename ) ) )
            {
                std::cerr << "ERROR: bases file \""<<baseFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Bases file: "<< baseFilename << std::endl;
            }
        }
        else if( verbose )
        {
            std::cout << "No bases file stated, getting bases from tree structure"<<std::endl;

        }

        if ( variableMap.count( "image" ) )
        {
            if( !boost::filesystem::is_regular_file(boost::filesystem::path( imageFilename ) ) )
            {
                std::cerr << "ERROR: partition image file \""<<imageFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Partition image file: "<< imageFilename << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no partition image file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if ( variableMap.count( "outputf" ) )
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path( outputFolder ) ) )
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


        std::string logFilename(outputFolder+"/"+progName+"_log.txt");
        std::ofstream logFile(logFilename.c_str());
        if(!logFile)
        {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }


        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile <<"Tree file:\t"<< treeFilename <<std::endl;
        logFile <<"Bases file:\t"<< baseFilename <<std::endl;
        logFile <<"Partition Image file:\t"<< imageFilename <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        if( niftiMode )
        {
            logFile << "Using nifti file format" << std::endl;
        }
        else
        {
            logFile << "Using vista file format" << std::endl;
        }
        logFile <<"-------------"<<std::endl;

        /////////////////////////////////////////////////////////////////


        image2treeBuilder builder( imageFilename, treeFilename, verbose, baseFilename );
        logFile <<" Roi size:\t" << builder.roiSize() << std::endl;
        builder.log( &logFile );
        builder.setOutputFolder( outputFolder );
        builder.importImagePart();

        /////////////////////////////////////////////////////////////////

        // save and print total time
        time_t programEndTime( time( NULL ) );
        int totalTime( difftime(programEndTime, programStartTime ) );
        std::cout <<" Program Finished, total time: " << totalTime/3600 << "h " <<  ( totalTime % 3600 ) / 60 << "' " << ( ( totalTime % 3600 ) % 60) <<"\"   "<< std::endl;
        logFile << "-------------" << std::endl;
        logFile << "Finish Time:\t" << ctime( &programEndTime ) << std::endl;
        logFile << "Elapsed time : " << totalTime / 3600 << "h " <<  ( totalTime % 3600 ) / 60 << "' " << ( ( totalTime % 3600 ) % 60 ) << "\"" << std::endl;


        // create file that indicates process was finished successfully
        std::string successFilename(outputFolder+"/success.txt");
        std::ofstream successFile(successFilename.c_str());
        if(!successFile)
        {
            std::cerr << "ERROR: unable to create success file: \""<<successFile<<"\""<<std::endl;
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

