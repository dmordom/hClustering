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


typedef enum
{
    CRSP_DIRECT,
    CRSP_SIMPLE,
    CRSP_MERGE,
    CRSP_RAND
}
CRSP_MODE;

bool verbose(false);

// A helper function to simplify the main part.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

int main( int argc, char *argv[] )
{



//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("matchpartition");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string treeFilename1, treeFilename2, matchTableFilename, outputFolder;
        unsigned int threads(0), levelDepth(1);
        float lambda(0);
        CRSP_MODE compMode(CRSP_DIRECT);
        bool partMatching(false), colorMatching(false), newPartMatching(false);


        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("input,i",  boost::program_options::value< std::string >(&treeFilename1), "file with input partitioned tree")
                ("target,t",  boost::program_options::value< std::string >(&treeFilename2), "file with tree where to find matchign partitions")
                ("matching,m",  boost::program_options::value< std::string >(&matchTableFilename), "file with leaf-matching table")
                ("partm,r",  boost::program_options::value< float >(&lambda), "Partition matching, inster labda coefficient value")
                ("new,n", "New Partition matching")
                //("depth,d", boost::program_options::value< unsigned int >(&levelDepth), "matching search depth (default = 1)")
                ("colorm,c",  "Color matching")
                ("output,o",  boost::program_options::value< std::string >(&outputFolder), "output folder where results will be written")
                ("verbose,v", "verbose option")
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



        if (variableMap.count("input")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename1))) {
                std::cerr << "ERROR: tree1 file \""<<treeFilename1<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tree1 file: "<< treeFilename1 << std::endl;
        } else {
            std::cerr << "ERROR: no tree1 file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("target")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename2))) {
                std::cerr << "ERROR: tree2 file \""<<treeFilename2<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tree2 file: "<< treeFilename2 << std::endl;
        } else {
            std::cerr << "ERROR: no tree1 file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("matching")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(matchTableFilename))) {
                std::cerr << "ERROR: match table file \""<<matchTableFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Match table file: "<< matchTableFilename << std::endl;
        } else {
            std::cerr << "ERROR: no match Table file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("partm"))
        {
            std::cout << "Performing partition and color matching: " << std::endl;
            std::cout <<"Using a lambda factor of "<< lambda << std::endl;
            partMatching = true;
            colorMatching = true;
        }
        else if (variableMap.count("new"))
        {
            std::cout << "Performing new matching (and color matching): " << std::endl;
            newPartMatching = true;
            colorMatching = true;
        }
        else if (variableMap.count("colorm"))
        {
            std::cout << "Performing only color matching: " << std::endl;
            colorMatching = true;
        }
        else
        {
            std::cerr << "ERROR: no matching type selected, select color or partition matching"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

//        if( partMatching || newPartMatching )
//        {
//            if( levelDepth > 5 )
//            {
//                std::cout << "Level depth indicated: " << levelDepth << " is too high, setting to a maximum of 5" << std::endl;
//                levelDepth = 5;
//            }
//            std::cout << "Using a search depth of: " << levelDepth << std::endl;
//        }


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
        logFile <<"Input partitioned tree:\t"<< treeFilename1 <<std::endl;
        logFile <<"Target tree:\t"<< treeFilename2 <<std::endl;
        logFile <<"Matching table:\t"<< matchTableFilename <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
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
            throw std::runtime_error ("ERROR @ compareTrees() datasets have different dimensions");
        }

        if (verbose) {
            std::cout <<"Tree1: "<< tree1.getReport(false) <<std::endl;
            std::cout <<"Tree2: "<< tree2.getReport(false) <<std::endl;
        }


        partitionMatcher matcher(&tree1,&tree2,matchTableFilename, verbose);

        std::cout <<matcher.reportBaseNodes()<<std::endl;
        std::string suffixPart("_pM_l" + string_utils::toString< float >( lambda ) + ".txt");
        std::string suffixNew("_pNew" + string_utils::toString< float >( lambda ) + ".txt");
        std::string suffixColor("_cM.txt");
        bool tree1Changed( false );


        if( partMatching )
        {
            matcher.findMatchingPartitions( lambda );
        }
        else if ( newPartMatching )
        {
            matcher.findMatchingPartitions( -1 );
        }

        if( colorMatching )
        {
            tree1Changed = matcher.matchColors();
        }


        if ( partMatching )
        {
            tree2.writeTree( outputFolder + "/" + tree2.getName() + suffixPart );
        }
        else if ( newPartMatching )
        {
            tree2.writeTree( outputFolder + "/" + tree2.getName() + suffixNew );
        }
        else
        {
            tree2.writeTree( outputFolder + "/" + tree2.getName() + suffixColor );
        }

        if( tree1Changed )
        {
            tree1.writeTree( outputFolder + "/" + tree1.getName() + suffixColor );
        }

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

