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

// parallel execution
#include <omp.h>

size_t numComps(0);

bool verbose(false);

// A helper function to simplify the main part.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

typedef enum
{
    PROC_STND,
    PROC_BASE,
    PROC_B2L,
    PROC_CLPS
}
PROCMODE;

int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("partitiontree");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");
        unsigned int threads(0), levelDepth(1), filterRadius(2);

        // program parameters
        std::string treeFilename, outputFolder;

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("verbose,v", "verbose option")
                ("tree-file,t",  boost::program_options::value< std::string >(&treeFilename), "file with the tree file")
                ("output,o",  boost::program_options::value< std::string >(&outputFolder), "output folder where processed tree(s) will be written")
                ("maxgran,m", "obtain only the max granularity partition number")
                ("hoz", "obtain only hoz partitions")

                ("depth,d", boost::program_options::value< unsigned int >(&levelDepth), "matching search depth (default = 1)")
//                ("radius,r", boost::program_options::value< unsigned int >(&filterRadius), "filter kernel radius (default = 2)")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "fnumber of processing threads to run the program in parallel, default: all available")
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
            std::cout << visibleOptions << std::endl;
            exit(0);
        }
        if (variableMap.count("version"))
        {
            std::cout << progName <<", version 1.0"<<std::endl;
            exit(0);
        }
        if (variableMap.count("verbose"))
        {
            std::cout << "verbose output"<<std::endl;
            verbose=true;
        }

        if (variableMap.count("threads")) {
            if (threads==1) {
                std::cout <<"Using a single processor"<< std::endl;
            } else if(threads==0 || threads>=omp_get_num_procs()){
                threads = omp_get_num_procs();
                std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
            } else {
                std::cout <<"Using a maximum of "<< threads <<" processors "<< std::endl;
            }
            omp_set_num_threads( threads );
        } else {
            threads = omp_get_num_procs();
            omp_set_num_threads( threads );
            std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
        }

        if (variableMap.count("tree-file"))
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




        if (variableMap.count("maxgran"))
        {
            WHtree tree(treeFilename);
            std::cout<<tree.getReport( false )<<std::endl;
            WHtreePartition treePartition(&tree);
            std::vector<size_t > maxpart;
            treePartition.level2granularity(&maxpart);
            std::cout<<"maxgranpart size: "<<std::endl<<maxpart.size()<<std::endl;
            return 0;
        }




        if (variableMap.count("output"))
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


        if( levelDepth > 5 )
        {
            std::cout << "Level depth indicated: " << levelDepth << " is too high, setting to a maximum of 5" << std::endl;
            levelDepth = 5;
        }
        std::cout << "Using a search depth of: " << levelDepth << std::endl;

//        if( filterRadius > 10 )
//        {
//            std::cout << "filter radius indicated: " << filterRadius << " is too high, setting to a maximum of 10" << std::endl;
//            filterRadius = 10;
//        }
//        std::cout << "Using a filter radius of: " << filterRadius << std::endl;

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

        WHtree tree(treeFilename);

        logFile << tree.getReport( false ) <<std::endl;
        std::cout<<tree.getReport( false )<<std::endl;

        std::vector< float > partitionValues;
        std::vector< std::vector< size_t> > partitionVector;

        WHtreePartition treePartition(&tree);

        if (variableMap.count("hoz"))
        {

            std::cout <<"getting hoz partitions at all levels..." <<std::endl;
            treePartition.scanHozPartitions( partitionValues, partitionVector );

            std::cout << partitionValues.size() << " Partitions obtained, writing to file..." <<std::endl;
            logFile <<"Initial partitions:\t"<< partitionValues.size() <<std::endl;
            std::string outPartFilename( outputFolder + "/allHozParts.txt" );
            treePartition.writePartitionSet( outPartFilename, partitionValues, partitionVector);


            std::vector< float > filtPartValues( partitionValues );
            std::vector< std::vector< size_t> > filtPartVector( partitionVector );

            std::cout << "Filtering with a radius of 100..." <<std::endl;
            treePartition.filterMaxPartitions( 100, &filtPartValues, &filtPartVector );

            std::cout << filtPartValues.size() << " Filtered partitions obtained, writing to file..." <<std::endl;
            logFile <<"Filtered partitions:\t"<< filtPartValues.size() <<std::endl;
            outPartFilename = ( outputFolder + "/hozFiltParts_r100.txt" );
            treePartition.writePartitionSet(outPartFilename, filtPartValues, filtPartVector);
            return 0;

        }


        std::cout <<"getting optimal partitions at all levels..." <<std::endl;
        treePartition.scanOptimalPartitions( levelDepth, partitionValues, partitionVector );

        std::cout << partitionValues.size() << " Partitions obtained, writing to file..." <<std::endl;
        logFile <<"Initial partitions:\t"<< partitionValues.size() <<std::endl;
        std::string outPartFilename( outputFolder + "/allParts_d" + boost::lexical_cast<std::string>(levelDepth) + ".txt" );
        treePartition.writePartitionSet( outPartFilename, partitionValues, partitionVector);

        {
            // get a tree with the best parttiions for 15 50 100 and 200 clusters
            std::vector< size_t > clustNums;
            clustNums.reserve( 4 );
            clustNums.push_back( 15 );
            clustNums.push_back( 50 );
            clustNums.push_back( 100 );
            clustNums.push_back( 250 );
            size_t cnIndex(0);
            std::vector< float > selPartValues;
            std::vector< std::vector< size_t> > selPartVector;
            selPartValues.reserve( clustNums.size() );
            selPartVector.reserve( clustNums.size() );
            for(size_t i=0; i< partitionVector.size(); ++i)
            {
                if( cnIndex >= clustNums.size())
                {
                    break;
                }

                if ( partitionVector[i].size() < clustNums[cnIndex])
                {
                    continue;
                }
                else
                {
                    selPartValues.push_back(partitionValues[i]);
                    selPartVector.push_back(partitionVector[i]);
                    ++cnIndex;
                }
            }
            tree.insertPartitions( selPartVector, selPartValues );
            std::string outTreeFilename( outputFolder + "/" + tree.getName() + "_parts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outTreeFilename += ( "_num4.txt" );
            tree.writeTree( outTreeFilename );

            selPartVector.erase(selPartVector.end()-1);
            selPartValues.erase(selPartValues.end()-1);
            tree.insertPartitions( selPartVector, selPartValues );
            outTreeFilename = ( outputFolder + "/" + tree.getName() + "_parts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outTreeFilename += ( "_num3.txt" );
            tree.writeTree( outTreeFilename );

        }


        std::vector < unsigned int > filterRadii;
        filterRadii.reserve( 6 );
        //        filterRadii.push_back( 1 );
        //        filterRadii.push_back( 2 );
        //        filterRadii.push_back( 5 );
        //        filterRadii.push_back( 10 );
        //        filterRadii.push_back( 15 );
        //        filterRadii.push_back( 20 );
        filterRadii.push_back( 100 );



        for(size_t i=0; i< filterRadii.size(); ++i)
        {
            std::vector< float > filtPartValues( partitionValues );
            std::vector< std::vector< size_t> > filtPartVector( partitionVector );

            std::cout << "Filtering with a radius of "<< filterRadii[i] << "..." <<std::endl;
            treePartition.filterMaxPartitions( filterRadii[i], &filtPartValues, &filtPartVector );

            std::cout << filtPartValues.size() << " Filtered partitions obtained, writing to file..." <<std::endl;
            logFile <<"Filtered partitions:\t"<< filtPartValues.size() <<std::endl;
            outPartFilename = ( outputFolder + "/filtParts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outPartFilename += ( "_r" + boost::lexical_cast<std::string>(filterRadii[i]) +  ".txt" );
            treePartition.writePartitionSet(outPartFilename, filtPartValues, filtPartVector);

            std::cout << "Adding filtered partitions to tree and writing..." <<std::endl;

            std::string outTreeFilename( outputFolder + "/" + tree.getName() + "_parts_d" + boost::lexical_cast<std::string>(levelDepth) );
            outTreeFilename += ( "_r" + boost::lexical_cast<std::string>(filterRadii[i]) +  ".txt" );

            tree.insertPartitions( filtPartVector, filtPartValues );
            tree.writeTree( outTreeFilename );
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


//    }
//    catch(std::exception& e)
//    {
//        std::cout << e.what() << std::endl;
//        return 1;
//    }
    return 0;
}
