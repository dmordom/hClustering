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

        std::string progName("fulltract");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string treeFilename;
        std::string maskFilename;
        std::string singleTractFolder;
        std::string meanTractFolder;
        std::string outputFolder;
        std::vector<size_t> inputNodes, inputLeaves, extraNodes;
        unsigned int threads(0);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("tree-file,t",  boost::program_options::value< std::string >(&treeFilename), "file with the hierarchical tree")
                ("input-nodes,n", boost::program_options::value< std::vector<size_t> >(&inputNodes), "input nodes to compute mean tracts of")
                ("input-leaves,l", boost::program_options::value< std::vector<size_t> >(&inputLeaves), "input leaf to compute full tract of")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "number of processing threads to run the program in parallel, default: all available")
                ("mask-file,m",  boost::program_options::value< std::string >(&maskFilename), "tractogram mask file")
                ("singlet-folder,s", boost::program_options::value< std::string >(&singleTractFolder), "folder with the single-voxel probabilistic tracts")
                ("meant-folder,f", boost::program_options::value< std::string >(&meanTractFolder), "folder with the node mean tracts")
                ("output,o",      boost::program_options::value< std::string >(&outputFolder), "output folder")
                ("verbose,v", "verbose option")
                ;

        // Hidden options, will be allowed both on command line and in config file, but will not be shown to the user.
        boost::program_options::options_description hiddenOptions("Hidden options");
        hiddenOptions.add_options()
                ("rest", boost::program_options::value< std::vector<size_t> >())
                ;
        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions).add(configOptions).add(hiddenOptions);
        boost::program_options::options_description configFileOptions;
        configFileOptions.add(configOptions).add(hiddenOptions);
        boost::program_options::options_description visibleOptions("Allowed options");
        visibleOptions.add(genericOptions).add(configOptions);
        boost::program_options::positional_options_description posOpt; //this arguments do not need to specify the option descriptor when typed in
        posOpt.add("rest", -1);

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
            } else if(threads==0 || threads>=omp_get_max_threads()){
                threads = omp_get_max_threads();
                std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
            } else {
                std::cout <<"Using a maximum of "<< threads <<" processors "<< std::endl;
            }
            omp_set_num_threads( threads );
        } else {
            threads = omp_get_max_threads();
            omp_set_num_threads( threads );
            std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
        }
        if(variableMap.count("input-leaves") && variableMap.count("input-nodes")) {
            std::cerr << "ERROR: only one input option (leaves or nodes) must be stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if(variableMap.count("rest"))
            extraNodes=(variableMap["rest"].as<std::vector<size_t> >());

        if (variableMap.count("input-leaves")) {
            inputLeaves.insert(inputLeaves.begin(),extraNodes.begin(),extraNodes.end());
            std::cout << "Input leaves: "<< inputLeaves << std::endl;
        } else if (variableMap.count("input-nodes")) {
            inputNodes.insert(inputNodes.begin(),extraNodes.begin(),extraNodes.end());
            std::cout << "Input nodes are: "<< inputNodes << std::endl;
        } else {
            std::cerr << "ERROR: no input nodes/leaves stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }



        if (variableMap.count("tree-file")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename))) {
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
        if (variableMap.count("mask-file")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(maskFilename))) {
                std::cerr << "ERROR: mask file \""<<maskFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tractogram mask file: "<< maskFilename << std::endl;
        } else {
            std::cerr << "ERROR: no tract mask file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if (variableMap.count("meant-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(meanTractFolder))) {
                std::cerr << "ERROR: mean tract folder \""<<meanTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "mean tracts folder: "<< meanTractFolder << std::endl;
        }
        if (variableMap.count("singlet-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(singleTractFolder))) {
                std::cerr << "ERROR: single tract folder \""<<singleTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Single tracts folder: "<< singleTractFolder << std::endl;
        }
        if (!variableMap.count("meant-folder") && !variableMap.count("singlet-folder")) {
            std::cerr << "ERROR: no tract folder stated"<<std::endl;
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

        std::string logFilename(outputFolder+"/"+progName+"_log.txt");
        std::ofstream logFile(logFilename.c_str());
        if(!logFile) {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }

        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        if (!inputNodes.empty())
            logFile <<"Input nodes:\t"<< inputNodes <<std::endl;
        if (!inputLeaves.empty())
            logFile <<"Input leaves:\t"<< inputNodes <<std::endl;
        logFile <<"Tree file:\t"<< treeFilename <<std::endl;
        logFile <<"Mask file:\t"<< maskFilename <<std::endl;
        logFile <<"Single tract folder:\t"<< singleTractFolder <<std::endl;
        logFile <<"Mean tract folder:\t"<< meanTractFolder <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Processors used:\t"<< threads <<std::endl;
        logFile <<"-------------"<<std::endl;

        // ==============================================================================

        WHtree tree(treeFilename);
        if (verbose) {
            std::cout<<tree.getReport()<<std::endl;
        }
        logFile <<tree.getReport()<<std::endl;

        treeManager treeMngr(&tree, verbose);
        treeMngr.log(&logFile);

        treeMngr.setMaskFilename(maskFilename);
        if(!meanTractFolder.empty())
            treeMngr.setMeanTractFolder(meanTractFolder);
        if(!singleTractFolder.empty())
            treeMngr.setSingleTractFolder(singleTractFolder);
        treeMngr.setFullTractFolder(outputFolder);

        std::vector<nodeID_t> input;
        for(size_t i(0); i<inputLeaves.size(); ++i)
            input.push_back(std::make_pair(false,inputLeaves[i]));
        for(size_t i(0); i<inputNodes.size(); ++i)
            input.push_back(std::make_pair(true,inputNodes[i]));

        treeMngr.writeFullTract(input);


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



