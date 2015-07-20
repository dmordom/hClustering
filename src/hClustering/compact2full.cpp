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

#define NEWBOOST true


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

        std::string progName("compact2full");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string maskFilename;
        std::string outputFilename, inputFilename;
        unsigned int threads(0);
        bool uFloat( false );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("mask-file,m",  boost::program_options::value< std::string >(&maskFilename), "tractogram mask file")
                ("input,i",      boost::program_options::value< std::string >(&inputFilename), "output folder")
                ("output,o",      boost::program_options::value< std::string >(&outputFilename), "output folder")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "number of processing threads to run the program in parallel, default: all available")
                ("verbose,v", "verbose option")
                ("ufloat,f", "use float representation to write tracts")
                ;

        // Hidden options, will be allowed both on command line and in config file, but will not be shown to the user.
        boost::program_options::options_description hiddenOptions("Hidden options");
        hiddenOptions.add_options()
                ;
        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions).add(configOptions).add(hiddenOptions);
        boost::program_options::options_description configFileOptions;
        configFileOptions.add(configOptions).add(hiddenOptions);
        boost::program_options::options_description visibleOptions("Allowed options");
        visibleOptions.add(genericOptions).add(configOptions);
        boost::program_options::positional_options_description posOpt; //this arguments do not need to specify the option descriptor when typed in
        //posOpt.add("tracts", -1);

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

        if (variableMap.count("ufloat")) {
            std::cout << "writing in float"<<std::endl;
            uFloat=true;
        }
        else {
            std::cout << "writing in char"<<std::endl;
            uFloat=false;
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

        if (variableMap.count("input")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(inputFilename))) {
                std::cerr << "ERROR: tract file \""<<inputFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tractogram file: "<< inputFilename << std::endl;
        } else {
            std::cerr << "ERROR: no input tract file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
        }

        if (variableMap.count("output")) {

            std::cout << "Output folder: "<< outputFilename << std::endl;
        } else {
            std::cerr << "ERROR: no output folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }



        // ==============================================================================

        vistaManager vistaMngr( outputFilename );


        vistaMngr.readAsLog();
        vistaMngr.readAsUnThres();
        vistaMngr.loadMask( maskFilename );
        vistaMngr.storeUnzipped();

        if (uFloat)
        {
            vistaMngr.writeInFloat();
        }
        else
        {
            vistaMngr.writeInChar();
        }




            std::string tractFilename( inputFilename );
            std::string tractFilenameZipped( tractFilename );
            std::string extension;

            #if NEWBOOST
                extension = boost::filesystem::path( tractFilename ).extension().string();
            #else
                extension = boost::filesystem::path( tractFilename ).extension();
            #endif

            std::string tractName;

            if ( extension == ".gz" )
            {
                tractFilename.resize(tractFilename.size()-3);

                // Variables for calling system commands
                std::string sysCommand( "gzip -dcf " + tractFilenameZipped + " > " + tractFilename );
                int success(system(sysCommand.c_str()));
                #if NEWBOOST
                    extension = boost::filesystem::path( tractFilename ).extension().string();
                #else
                    extension = boost::filesystem::path( tractFilename ).extension();
                #endif
            }


            if ( extension == ".v" )
            {

            }
            else
            {
                std::cerr << "File \"" << tractFilenameZipped << "\" has no recognized extension (\"." << extension << "\")." << std::endl;
                exit -1;
            }


            compactTract simpleTract;
            vistaMngr.readStringTract(tractFilename,simpleTract );
            vistaMngr.storeFullTract(outputFilename,simpleTract);

            if ( extension == ".gz" )
            {
                // Variables for calling system commands
                std::string sysCommand( "rm -f " + tractFilename );
                int success(system(sysCommand.c_str()));
            }

        // ==============================================================================

        // save and print total time
        time_t programEndTime(time(NULL));
        int totalTime( difftime(programEndTime,programStartTime) );
        std::cout << "Program Finished, total time: " << totalTime/3600 <<"h "<<  (totalTime%3600)/60 <<"' "<< ((totalTime%3600)%60) <<"\""<< std::endl;


//    }
//    catch(std::exception& e)
//    {
//        std::cout << e.what() << std::endl;
//        return 1;
//    }
    return 0;
}



