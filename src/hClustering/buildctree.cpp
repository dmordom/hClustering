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
#include "cnbTreeBuilder.h"
#include "randCnbTreeBuilder.h"

size_t numComps(0);


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

        std::string progName("buildtree");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, inputFolder, outputFolder, meanTractFolder;
        float memory(2), maxNbDist( 0.1 );
        unsigned int nbLevel(0), threads(0);
        bool keepDiscarded(false), randomtracts(false), verbose(false);
        TC_GROWTYPE growType(TC_GROWOFF);
        size_t baseSize( 0 );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("roi-file,r",  boost::program_options::value< std::string >(&roiFilename), "file with the seed voxels coordinates")
                ("cnbhood,c",  boost::program_options::value< unsigned int >(&nbLevel), "use centroid method with C neighborhood level")
                ("keep-discarded,k", "keep discarded voxels in the discarded section of the tree")
                ("verbose,v", "verbose option")
                ("maxnbdist,d",  boost::program_options::value< float >(&maxNbDist), "maximum distance a voxel can have to not be discarded")
                ("basesize,s",  boost::program_options::value< size_t >(&baseSize), "grow homogeneous base nodes of size S")
                ("basenum,n",  boost::program_options::value< size_t >(&baseSize), "grow N homogeneous base nodes")
                ("rand", "small randmom tracts")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "fnumber of processing threads to run the program in parallel, default: all available")
                ("input,i",  boost::program_options::value< std::string >(&inputFolder), "input data folder (single tractograms)")
                ("meantract-folder,t",  boost::program_options::value< std::string >(&meanTractFolder), "mean tract folder location for centroid method")
                ("output,o",  boost::program_options::value< std::string >(&outputFolder), "output folder where tree will be written")
                ("cache-memory,m",  boost::program_options::value< float >(&memory), "maximum of memory (in GBytes) to use for tractogram cache memory")
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
        if (variableMap.count("keep-discarded")) {
            std::cout << "keep discarded voxels"<<std::endl;
            keepDiscarded=true;
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

        if (variableMap.count("cnbhood")) {
            if ( (nbLevel!=6)&&(nbLevel!=18)&&(nbLevel!=26)&&(nbLevel!=32)&&(nbLevel!=92)&&(nbLevel!=124) ) {
                std::cerr << "ERROR: invalid nbhood level, only (6,18,26,32,92,124) are accepted"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "Centroid method. "<< nbLevel <<" neighborhood"<< std::endl;
        } else {
            std::cerr << "ERROR: no output neighborhood level"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if ( ( variableMap.count("basesize") && variableMap.count("basenum") ) )
        {
            std::cerr << "ERROR: options --basesize --basenum are mutually exclusive"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("rand")) {
            std::cout << "random tracts"<<std::endl;
            randomtracts=true;
        }
        if (variableMap.count("basesize")) {
            std::cout << "Growing homogeneous base nodes of size: "<< baseSize <<std::endl;
            growType = TC_GROWSIZE;
        }
        if (variableMap.count("basenum")) {
            std::cout << "Growing "<< baseSize <<" homogeneous base nodes" <<std::endl;
            growType = TC_GROWNUM;
        }


        if (variableMap.count("meantract-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(meanTractFolder))) {
                std::cerr << "ERROR: meantract folder \""<<meanTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "meanTractFolder folder: "<< meanTractFolder << std::endl;
        } else if(!randomtracts) {
            std::cerr << "ERROR: no meanTractFolder folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (memory<0.1 || memory>50) {
            std::cerr << "ERROR: cache size must be a positive float between 0.1 and 50"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        std::cout <<"Tractogram cache memory: "<< memory <<" GBytes"<< std::endl;



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
        logFile <<"Method used:\tCentroid "<<nbLevel<<" nbhood" <<std::endl;
        logFile <<"Mean tract folder:\t"<< meanTractFolder <<std::endl;
        logFile <<"Cache size:\t"<< memory <<" GB"<<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Processors used:\t"<< threads <<std::endl;
        logFile <<"-------------"<<std::endl;

        /////////////////////////////////////////////////////////////////


        switch(growType)
        {
        case TC_GROWOFF:
            logFile << "Region growing: None" <<std::endl;
            break;
        case TC_GROWSIZE:
            logFile << "Region growing: Size: "<< baseSize <<std::endl;
            break;
        case TC_GROWNUM:
            logFile << "Region growing: Number: "<< baseSize <<std::endl;
            break;
        }

        if( randomtracts )
        {
            randCnbTreeBuilder builderRand(roiFilename, verbose);
            logFile <<"Roi size:\t"<< builderRand.roiSize() <<std::endl;

            builderRand.log(&logFile);
            builderRand.setInputFolder(inputFolder);
            builderRand.setOutputFolder(outputFolder);
            builderRand.buildCentroiRand(nbLevel,memory,keepDiscarded,growType,baseSize);
        }
        else
        {

            cnbTreeBuilder builder(roiFilename, verbose);
            logFile <<"Roi size:\t"<< builder.roiSize() <<std::endl;

            builder.log(&logFile);
            builder.setInputFolder(inputFolder);
            builder.setOutputFolder(outputFolder);
            builder.setMaxNbDist( maxNbDist );
            builder.buildCentroid(nbLevel,memory,meanTractFolder,keepDiscarded,growType,baseSize);
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

