// std librabry
#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <stdexcept>
#include <fstream>
#include <algorithm>

// boost library
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

// parallel execution
#include <omp.h>

// classes
#include "WFileParser.h"
#include "WHcoord.h"
#include "vistaManager.h"
#include "compactTract.h"

bool verbose(false);




int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("randtracts");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, outputFolder;
        unsigned int threads(0);
        size_t dimension(0);
        bool method2(false);
        size_t randSeed(0);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("verbose,v", "verbose")
                ("method-2,2", "use alternative method")
                ("roi-file,r",  boost::program_options::value< std::string >(&roiFilename), "file with the seed voxels coordinates")
                ("dimension,d",  boost::program_options::value< size_t >(&dimension), "dimension of the random tractograms")
                ("seed,s",  boost::program_options::value< size_t >(&randSeed), "random number generator seed")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "fnumber of processing threads to run the program in parallel, default: all available")
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

        if (variableMap.count("method-2")) {
            std::cout << "Using alternative method"<<std::endl;
            method2=true;
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

        if (variableMap.count("dimension")) {
            std::cout << "Tractograms dimension: "<< dimension << std::endl;
        } else {
            std::cerr << "ERROR: tractogram dimension stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }

        std::cout << "Seed: "<< randSeed << std::endl;


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
        logFile <<"Roi file:\t"<< roiFilename <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Tractogram dimension:\t"<< dimension << std::endl;
        logFile <<"Seed:\t"<< randSeed << std::endl;
        if(method2)
            logFile <<"Using method 2"<<std::endl;
        else
            logFile <<"Using method 1"<<std::endl;

        logFile <<"-------------"<<std::endl;


        //====================== Read Roi ============================

        std::vector<WHcoord> roivect;

        {
            WFileParser parser( roiFilename );
            if( !parser.readFile() ) {
                std::cerr<<"ERROR Parser error when reading roi"<<std::endl;
                return false;
            }
            std::vector<std::string>lines = parser.getRawLines();
            if( lines.size() == 0 ) {
                std::cerr<<"ERROR roi file is empty"<<std::endl;
                return false;
            }

            {
                std::vector< std::vector< std::string> >coordStrings = parser.getLinesForTagSeparated( "roi" );
                if( coordStrings.size() == 0 ) {
                    std::cerr<<"ERROR @ treeBuilder::readRoi(): no roi coordinates in roi file (lacking #roi tag?)"<<std::endl;
                    return false;
                }
                roivect.reserve(coordStrings.size());
                for( size_t i = 0; i < coordStrings.size(); ++i ) {
                    WHcoord tempCoord(boost::lexical_cast<coord_t>(coordStrings[i][0]) ,
                                         boost::lexical_cast<coord_t>(coordStrings[i][1]) ,
                                         boost::lexical_cast<coord_t>(coordStrings[i][2]) );
                    roivect.push_back(tempCoord);
                }
            }
            std::sort(roivect.begin(),roivect.end());

            if(verbose)
                std::cout<<"Roi loaded, "<< roivect.size() <<" seed voxels"<<std::endl;
        }


        //==================== Generate Tracts ======================

        time_t loopStart(time(NULL)), lastTime(time(NULL));
        volatile size_t threadCount(0);
        vistaManager vMngr(outputFolder);
        vMngr.writeInChar();

        // This is the underlying integer random number generator
        boost::mt19937 igen;
        igen.seed(randSeed);

        // normally distributed random number generator mean 0 variance 1
        boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  randnPrng(igen, boost::normal_distribution<>(0,1));

        // uniformly distributed random number generator between 0 and 1
        boost::variate_generator<boost::mt19937, boost::uniform_01<> >  u01Prng(igen, boost::uniform_01<>());


        if (verbose)
            std::cout <<"Generating tractograms uniformly distributed over the hypersphere surface, Tract size: "<<dimension<<"..."<< std::endl;

        for (int i=0; i<roivect.size(); ++i) {

            std::vector<float> randTract;
            randTract.reserve(dimension);

            if (!method2) {

                float radius(0);
                // first create a vector of normally distributed values
                for (int j=0; j<dimension; ++j ){
                    float value(fabs(randnPrng()));
                    randTract.push_back((value));
                    radius+=value*value;
                }
                //calculate the lenght of the vector
                radius=sqrt(radius);
                // divide the vector by its radius and make absolute value to normalize it to the first quadrant of the unit hypersphere
                for (int j=0; j<randTract.size(); ++j ){
                    randTract[j]=(randTract[j]/radius);
                }
                // now the vector is part of a uniformly distributed set of vector over the positive quadrant surface of the unit hypersphere

            } else {

                bool notFound(true);
                while (notFound) {

                    randTract.clear();
                    float norm(0);

                    //first create a vector uniformly distributed over the hyperspace
                    for (int j=0; j<dimension; ++j ){
                        float value(u01Prng());
                        randTract.push_back((value));
                        norm+=value*value;
                    }
                    norm=sqrt(norm);
                    // if the vector lieas outside the unit hyperphere, discard it, if its inside, project it to the sphere
                    if (norm<=1) {
                        for (int j=0; j<dimension; ++j ){
                            randTract[j]=(randTract[j]/norm);
                        }
                        notFound=false;
                    }
                }
            }

            compactTract tract(randTract);
            vMngr.writeLeafTract(roivect[i],tract);

            if(verbose) {
                time_t currentTime(time(NULL));
                if(currentTime-lastTime>1) {
                    lastTime=currentTime;
                    float progress ((i*100.)/(roivect.size()));
                    std::string coutString(boost::lexical_cast<std::string>((int)progress)+" % of tractograms written ("+boost::lexical_cast<std::string>(i)+"). Expected remaining time: ");
                    if (progress>0) {
                        int expectedRemain(difftime(time(NULL),loopStart)*((100.-progress)/progress));
                        coutString+=boost::lexical_cast<std::string>(expectedRemain/3600) +"h "+ boost::lexical_cast<std::string>((expectedRemain%3600)/60) +"' "+ boost::lexical_cast<std::string>((expectedRemain%3600)%60) +"\"  ";
                    }
                    std::cout<<"\r"<<coutString<<std::flush;
                }
            } // end verbose

        } // end big loop

        while(threadCount!=0)  // wait until all threads have finished
            boost::this_thread::sleep(boost::posix_time::microseconds(100));

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

