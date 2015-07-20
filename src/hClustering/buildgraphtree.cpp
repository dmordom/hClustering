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
#include "graphTreeBuilder.h"

size_t numComps(0);

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

        std::string progName("buildgraphtree");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, inputFolder, outputFolder;
        unsigned int selector(0), threads(0);
        tgraph_t graphMethod;

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("roi-file,r",  boost::program_options::value< std::string >(&roiFilename), "file with the seed voxels coordinates")
                ("graph-method,g",  boost::program_options::value< unsigned int >(&selector), "use N graph method (0=single, 1=complete, 2=average, 3=weighted)")
                ("verbose,v", "verbose option")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "fnumber of processing threads to run the program in parallel, default: all available")
                ("input,i",  boost::program_options::value< std::string >(&inputFolder), "input data folder (distance matrix)")
                ("output,o",  boost::program_options::value< std::string >(&outputFolder), "output folder where tree will be written")
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

        if (variableMap.count("roi-file")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(roiFilename))) {
                std::cerr << "ERROR: roi file \""<<roiFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Roi voxels file: "<< roiFilename << std::endl;
        } else {
            std::cerr << "ERROR: no roi file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("input")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(inputFolder))) {
                std::cerr << "ERROR: input folder \""<<inputFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "input folder: "<< inputFolder << std::endl;
        } else {
            std::cerr << "ERROR: no input folder stated"<<std::endl;
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


        if (variableMap.count("graph-method")) {
            if ( (selector!=0)&&(selector!=1)&&(selector!=2)&&(selector!=3) ) {
                std::cerr << "ERROR: invalid graph method"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Graph method. "<< std::flush;
            if (selector==0) {
                graphMethod = TG_SINGLE;
                std::cout<<"Single linkage: D(k,i+j) = min[D(i,k),D(j,k)]"<<std::endl;
            } else if (selector==1) {
                graphMethod = TG_COMPLETE;
                std::cout<<"Complete linkage: D(k,i+j) = MAX[D(i,k),D(j,k)]"<<std::endl;
            } else if (selector==2) {
                graphMethod = TG_AVERAGE;
                std::cout<<"Average linkage: D(k,i+j) = [D(i,k)*Size(i),D(j,k)*size(j)]/[size(i)+size(j)]"<<std::endl;
            } else {
                graphMethod = TG_WEIGHTED;
                std::cout<<"Weighted linkage: D(k,i+j) = [D(i,k)+D(i,k)]/2"<<std::endl;
            }
        } else {
            std::cerr << "ERROR: no graph method stated"<<std::endl;
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
        logFile <<"Roi file:\t"<< roiFilename <<std::endl;
        logFile <<"Input folder:\t"<< inputFolder <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;

        if (graphMethod==TG_SINGLE)
            logFile <<"Method used:\tSingle linkage: D(k,i+j) = min[D(i,k),D(j,k)]"<<std::endl;
        else if (graphMethod==TG_COMPLETE)
            logFile <<"Method used:\tComplete linkage: D(k,i+j) = MAX[D(i,k),D(j,k)]"<<std::endl;
        else if (graphMethod==TG_AVERAGE)
            logFile <<"Method used:\tAverage linkage: D(k,i+j) = [D(i,k)*Size(i),D(j,k)*size(j)]/[size(i)+size(j)]"<<std::endl;
        else if (graphMethod==TG_WEIGHTED)
            logFile <<"Method used:\tWeighted linkage: D(k,i+j) = [D(i,k)+D(i,k)]/2"<<std::endl;
        else {
            std::cerr << "ERROR: unrecognized graph option"<<std::endl;
            exit(-1);
        }

        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Processors used:\t"<< threads <<std::endl;
        logFile <<"-------------"<<std::endl;

        /////////////////////////////////////////////////////////////////

        graphTreeBuilder builder(roiFilename, verbose);
        logFile <<"Roi size:\t"<< builder.roiSize() <<std::endl;

        builder.log(&logFile);
        builder.setInputFolder(inputFolder);
        builder.setOutputFolder(outputFolder);
        builder.buildGraph(graphMethod);

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

