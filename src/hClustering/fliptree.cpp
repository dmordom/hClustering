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

        std::string progName("fliptree");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string treeFilename, flippedTreeFolder;
        unsigned int threads(0);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("tree,t",  boost::program_options::value< std::string >(&treeFilename), "tree file")
                ("output,o",  boost::program_options::value< std::string >(&flippedTreeFolder), "output folder where flipped tree will be written")
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

        if (variableMap.count("tree")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename))) {
                std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tree file: "<< treeFilename << std::endl;
        } else {
            std::cerr << "ERROR: no tree file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("output")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(flippedTreeFolder))) {
                std::cerr << "ERROR: output folder \""<<flippedTreeFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "Output folder: "<< flippedTreeFolder << std::endl;
        } else {
            std::cerr << "ERROR: no output folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }

        std::string logFilename(flippedTreeFolder+"/"+progName+"_log.txt");
        std::ofstream logFile(logFilename.c_str());
        if(!logFile) {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }
        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Processors used:\t"<< threads <<std::endl;
        logFile <<"Tree file:\t"<< treeFilename <<std::endl;
        logFile <<"Output folder:\t"<< flippedTreeFolder <<std::endl;
        logFile <<"-------------"<<std::endl;

        /////////////////////////////////////////////////////////////////

        WHtree tree(treeFilename);
        if (verbose) {
            std::cout<<tree.getReport()<<std::endl;
        }
        logFile <<tree.getReport()<<std::endl;

        treeManager manager(&tree, verbose);
        manager.setOutputFolder(flippedTreeFolder);
        manager.flipX();
        manager.writeTree();

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

