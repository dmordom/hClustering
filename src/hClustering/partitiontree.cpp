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
//  partitiontree
//
//  Obtain tree partitions at all granularity levels using the Spread-Separation method (finding the the partition with highest SS index at each granularity).
//   Optimal SS value for each partition is searched within a defined search-depth hierarchical levels. Final partitions can be filtered with a defined kernel size.
//   to keep local SS maxima within that kernel. For SS index refer to (Moreno-Dominguez, 2014).
//   For an interactive 3D partition management with more options please use the Hierarchcial Clustering module developed in OpenWalnut (www.openwalnut.org).
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -t --tree:       File with the hierarchical tree to extract partitions from.
//
//   -O --outputf:    Output folder where partition files will be written.
//
//  [-d --search-depth]:  Search optimal partition for each granularity within d hierarchical levels.
//                         A higher value will produce more optimized partition but will increase computing time.
//                         Default: 3. Recommendened values: 3 for good quality and fast computation, 4 for enhanced quality.
//
//  [-r --filter-radius]: Filter output partitions to keep only local SS (partition quality) maxima
//                         within a r-sized kernel across the granularity dimension.
//
//  [-h --hoz]:       Write horizontal cut partitions instead of SS ones (optimal partition search is still based on SS index).
//
//  [-m --maxgran]:   Compute and write only the maximum granularity (meta-leaves) partition.
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        write output tree in vista coordinates (default is nifti).
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  * Usage example:
//
//   partitiontree -t tree_lh.txt -O results/ -d 3 -r 50 -v
//
//
//  * Outputs (in output folder defined at option -O):
//
//   (default outputs)
//   - 'allSSparts_dX.txt' - (where X is the search depth level defined at parameter -d) Contains a summary of the partition information (cut value and size) for all granularities.
//   - 'TREE_SSparts_dX.txt' - (where TREE is the filename of the input tree defined at parameter -t) contains a copy of the original tree file with the partitions at all granularities included in the relevant fields.
//   - 'partitiontree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program.
//
//   (additional if using option -r)
//   - 'filtSSparts_dX_rY.txt' - (where Y is the filter radius defined at parameter -r) Contains a summary of the resulting filtered partitions.
//   - 'TREE_SSparts_dX_rY.txt' - contains a copy of the original tree file with the resulting filtered partitions included in the relevant fields.
//
//   (when using --hoz option, the prefix 'SS' will be replaced by 'Hoz')
//
//   (alternative outputs when using option --maxgran)
//   - 'fmaxgranPart.txt' - Contains the size information of the resulting maximal granularity partition for that tree.
//   - 'TREE_maxgranPart.txt' - contains a copy of the original tree file with the resulting max granularity partition included in the relevant fields.
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
#include "WHtreePartition.h"
#include "fileManagerFactory.h"

// parallel execution
#include <omp.h>




int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("partitiontree");
        std::string configFilename("../../config/"+progName+".cfg");
        unsigned int threads(0), levelDepth(3), filterRadius(0);
        bool verbose(false), niftiMode( true );

        // program parameters
        std::string treeFilename, outputFolder;

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "tree,t",  boost::program_options::value< std::string >(&treeFilename), "file with the tree to compute partitions from")
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder where partition files will be written")
                ( "search-depth,d", boost::program_options::value< unsigned int >(&levelDepth)->implicit_value(3), "[opt] optimal partition search depth (default = 3)")
                ( "filter-radius,r", boost::program_options::value< unsigned int >(&filterRadius)->implicit_value(0), "[opt] output partition filter kernel radius (default = 0 | no filtering)")
                ( "hoz", "[opt] obtain horizontal cut partitions (instead of Spread-Separation ones)")
                ( "maxgran,m", "[opt] obtain only the maximum granularity partition")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "pthreads,p",  boost::program_options::value< unsigned int >(&threads), "[opt] number of processing threads to run the program in parallel, default: all available")
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
            std::cout << "partitiontree" << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << "Obtain tree partitions at all granularity levels using the Spread-Separation method (finding the the partition with highest SS index at each granularity)." << std::endl;
            std::cout << " Optimal SS value for each partition is searched within a defined search-depth hierarchical levels. Final partitions can be filtered with a defined kernel size." << std::endl;
            std::cout << " to keep local SS maxima within that kernel. For SS index refer to (Moreno-Dominguez, 2014)" << std::endl;
            std::cout << " For an interactive 3D partition management with more options please use the Hierarchcial Clustering module developed in OpenWalnut (www.openwalnut.org)." << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -t --tree:       File with the hierarchical tree to extract partitions from." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where partition files will be written." << std::endl << std::endl;
            std::cout << "[-d --search-depth]:  Search optimal partition for each granularity within d hierarchical levels." << std::endl;
            std::cout << "                       A higher value will produce more optimized partition but will increase computing time." << std::endl;
            std::cout << "                       Default: 3. Recommendened values: 3 for good quality and fast computation, 4 for enhanced quality." << std::endl << std::endl;
            std::cout << "[-r --filter-radius]: Filter output partitions to keep only local SS (partition quality) maxima" << std::endl;
            std::cout << "                       within a r-sized kernel across the granularity dimension." << std::endl << std::endl;
            std::cout << "[-h --hoz]:       Write horizontal cut partitions instead of SS ones (optimal partition search is still based on SS index)." << std::endl << std::endl;
            std::cout << "[-m --maxgran]:   Compute and write only the maximum granularity (meta-leaves) partition." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]: 	    write output tree in vista coordinates (default is nifti)." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " partitiontree -t tree_lh.txt -O results/ -d 3 -r 50 -v" << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " (default outputs)" << std::endl;
            std::cout << " - 'allSSparts_dX.txt' - (where X is the search depth level defined at parameter -d) Contains a summary of the partition information (cut value and size) for all granularities." << std::endl;
            std::cout << " - 'TREE_SSparts_dX.txt' - (where TREE is the filename of the input tree defined at parameter -t) contains a copy of the original tree file with the partitions at all granularities included in the relevant fields." << std::endl;
            std::cout << " - 'partitiontree_log.txt' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
            std::cout << std::endl;
            std::cout << " (additional if using option -r)" << std::endl;
            std::cout << " - 'filtSSparts_dX_rY.txt' - (where Y is the filter radius defined at parameter -r) Contains a summary of the resulting filtered partitions." << std::endl;
            std::cout << " - 'TREE_SSparts_dX_rY.txt' - contains a copy of the original tree file with the resulting filtered partitions included in the relevant fields." << std::endl;
            std::cout << std::endl;
            std::cout << " (when using --hoz option, the prefix 'SS' will be replaced by 'Hoz'')" << std::endl;
            std::cout << std::endl;
            std::cout << " (alternative outputs when using option --maxgran)" << std::endl;
            std::cout << " - 'fmaxgranPart.txt' - Contains the size information of the resulting maximal granularity partition for that tree." << std::endl;
            std::cout << " - 'TREE_maxgranPart.txt' - contains a copy of the original tree file with the resulting max granularity partition included in the relevant fields." << std::endl;
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

        if (variableMap.count("pthreads"))
        {
            if (threads==1)
            {
                std::cout <<"Using a single processor"<< std::endl;
            }
            else if(threads==0 || threads>=omp_get_num_procs())
            {
                threads = omp_get_num_procs();
                std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
            }
            else
            {
                std::cout <<"Using a maximum of "<< threads <<" processors "<< std::endl;
            }
            omp_set_num_threads( threads );
        }
        else
        {
            threads = omp_get_num_procs();
            omp_set_num_threads( threads );
            std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
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
            std::cout << "Roi voxels file: "<< treeFilename << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tree file stated"<<std::endl;
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



        if (variableMap.count("maxgran"))
        {
            std::cout<<"Obtaining only max. granularity partition..."<<std::endl;

            WHtree tree(treeFilename);
            std::cout<<tree.getReport( false )<<std::endl;
            if( tree.testRootBaseNodes() )
            {
                std::vector<size_t > maxpart( tree.getRootBaseNodes() );
                std::vector<std::vector<size_t > > partitionVector( 1, maxpart);
                std::vector<float > partitionValues(1,0);
                std::cout<<"maxgranpart size: "<<std::endl<<maxpart.size()<<std::endl;
                WHtreePartition partitioner(&tree);
                std::string outPartFilename( outputFolder + "/maxgranPart.txt" );
                partitioner.writePartitionSet( outPartFilename, partitionValues,partitionVector);
                tree.insertPartitions( partitionVector, partitionValues );
                std::string outTreeFilename( outputFolder + "/" + tree.getName() + "_maxgranPart" );
                outTreeFilename += ( ".txt" );
                tree.writeTree( outTreeFilename, niftiMode );
                return 0;
            }
            else
            {
                std::cout<<"ERROR: tree  does not have a maximum granularity meta-leaf partition"<<std::endl;
                return(-1);
            }
        }

        if( levelDepth > 5 )
        {
            std::cout << "Level depth indicated: " << levelDepth << " is too high, setting to a maximum of 5" << std::endl;
            levelDepth = 5;
        }
        std::cout << "Using a search depth of: " << levelDepth << std::endl;

        if( filterRadius > 1000 )
        {
            std::cout << "filter radius indicated: " << filterRadius << " is too high (max is 1000), setting to 100" << std::endl;
            filterRadius = 10;
        }
        if( filterRadius == 0 )
        {
            std::cout << "using no filtering (radius 0)" << std::endl;
        }
        else if( filterRadius < 0 )
        {
            std::cout << "filter radius indicated: " << filterRadius << " must be positive. using no filtering (radius 0)" << std::endl;
            filterRadius = 0;
        }
        else
        {
            std::cout << "Using a filter radius of: " << filterRadius << std::endl;
        }

        /////////////////////////////////////////////////////////////////



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
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        if( niftiMode )
        {
            logFile << "Using nifti file format" << std::endl;
        }
        else
        {
            logFile << "Using vista file format" << std::endl;
        }

        WHtree tree(treeFilename);

        logFile << tree.getReport( false ) <<std::endl;
        std::cout<<tree.getReport( false )<<std::endl;

        std::vector< float > partitionValues;
        std::vector< std::vector< size_t> > partitionVector;

        WHtreePartition treePartition(&tree);

        std::string prefix;

        if (variableMap.count("hoz"))
        {
            prefix = "Hoz";
            std::cout <<"getting hoz partitions at all levels..." <<std::endl;
            treePartition.scanHozPartitions( &partitionValues, &partitionVector );

            std::cout << partitionValues.size() << " Partitions obtained, writing to file..." <<std::endl;
            logFile <<"Initial partitions:\t"<< partitionValues.size() <<std::endl;
            std::string outPartFilename( outputFolder + "/all" + prefix + "parts.txt" );
            treePartition.writePartitionSet( outPartFilename, partitionValues, partitionVector);

            tree.insertPartitions( partitionVector, partitionValues );
            std::string outTreeFilename( outputFolder + "/" + tree.getName() + "_" + prefix + "parts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outTreeFilename += ( ".txt" );
            tree.writeTree( outTreeFilename, niftiMode );
        }
        else
        {

            prefix = "SS";
            std::cout <<"getting SS partitions at all levels..." <<std::endl;
            treePartition.scanOptimalPartitions( levelDepth, &partitionValues, &partitionVector );

            std::cout << partitionValues.size() << " Partitions obtained, writing to file..." <<std::endl;
            logFile <<"Initial partitions:\t"<< partitionValues.size() <<std::endl;
            std::string outPartFilename( outputFolder + "/all" + prefix + "parts_d" + boost::lexical_cast<std::string>(levelDepth) + ".txt" );
            treePartition.writePartitionSet( outPartFilename, partitionValues, partitionVector);

            tree.insertPartitions( partitionVector, partitionValues );
            std::string outTreeFilename( outputFolder + "/" + tree.getName() + "_" + prefix + "parts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outTreeFilename += ( ".txt" );
            tree.writeTree( outTreeFilename, niftiMode );

        }


        std::vector < unsigned int > filterRadii;
        //filterRadii.reserve( 6 );
        //        filterRadii.push_back( 1 );
        //        filterRadii.push_back( 2 );
        //        filterRadii.push_back( 5 );
        //        filterRadii.push_back( 10 );
        //        filterRadii.push_back( 15 );
        //        filterRadii.push_back( 20 );
        filterRadii.push_back( filterRadius );



        for(size_t i=0; i< filterRadii.size(); ++i)
        {
            if( filterRadii[i] <= 0 )
            {
                continue;
            }
            std::vector< float > filtPartValues( partitionValues );
            std::vector< std::vector< size_t> > filtPartVector( partitionVector );

            std::cout << "Filtering with a radius of "<< filterRadii[i] << "..." <<std::endl;
            treePartition.filterMaxPartitions( filterRadii[i], &filtPartValues, &filtPartVector );

            std::cout << filtPartValues.size() << " Filtered partitions obtained, writing to file..." <<std::endl;
            logFile <<"Filtered partitions:\t"<< filtPartValues.size() <<std::endl;
            std::string outPartFilename( outputFolder + "/filt" + prefix + "parts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outPartFilename += ( "_r" + boost::lexical_cast<std::string>(filterRadii[i]) +  ".txt" );
            treePartition.writePartitionSet(outPartFilename, filtPartValues, filtPartVector);

            std::cout << "Adding filtered partitions to tree and writing..." <<std::endl;

            std::string outTreeFilename( outputFolder + "/" + tree.getName() + "_" + prefix + "parts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outTreeFilename += ( "_r" + boost::lexical_cast<std::string>(filterRadii[i]) +  ".txt" );

            tree.insertPartitions( filtPartVector, filtPartValues );
            tree.writeTree( outTreeFilename, niftiMode );
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
