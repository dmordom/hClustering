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
#include "../common/treeManager.h"


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

        std::string progName("fliptracts");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, inputTractFolder, flippedTractFolder, maskFilename, treeFilename;
        unsigned int threads(0);
        std::vector<size_t> inNodes;
        bool roiMode(false), treeBaseMode(false), overWrite( false );
        size_t loopLength( 0 );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("verbose,v", "verbose option")
                ("overwrite,o", "tracts will be overwritten even if they already exist")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "fnumber of processing threads to run the program in parallel, default: all available")
                ("roi,r",  boost::program_options::value< std::string >(&roiFilename), "roi file")
                ("bases,b",  boost::program_options::value< std::string >(&treeFilename), "tree file to recover base nodes")
                ("input-folder,i",  boost::program_options::value< std::string >(&inputTractFolder), "folder with the input tractograms")
                ("flipped-folder,f",  boost::program_options::value< std::string >(&flippedTractFolder), "folder with the flipped tractograms will be written")
                ("mask,m",  boost::program_options::value< std::string >(&maskFilename), "tract mask file")
                ;

        // Hidden options, will be allowed both on command line and in config file, but will not be shown to the user.
        boost::program_options::options_description hiddenOptions("Hidden options");
        hiddenOptions.add_options()
                ( "mean-tracts", boost::program_options::value< std::vector<size_t> >())
                ;


        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions).add(configOptions).add(hiddenOptions);
        boost::program_options::options_description configFileOptions;
        configFileOptions.add(configOptions).add(hiddenOptions);
        boost::program_options::options_description visibleOptions("Allowed options");
        visibleOptions.add(genericOptions).add(configOptions);
        boost::program_options::positional_options_description posOpt; //this arguments do not need to specify the option descriptor when typed in
        posOpt.add("mean-tracts", -1);



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
        if (variableMap.count("overwrite")) {
            std::cout << "overwrite enabled"<<std::endl;
            overWrite=true;
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

        if (variableMap.count("roi")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(roiFilename))) {
                std::cerr << "ERROR: roi file \""<<roiFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Roi file: "<< roiFilename << std::endl;
            roiMode = true;
            if (variableMap.count("mean-tracts"))
            {
                std::cout << "WARNING: single tract roi mode, mean tract ids inserted will be ignored"<< inputTractFolder << std::endl;
            }
        }
        else  if (variableMap.count("bases")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename))) {
                std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tree file: "<< treeFilename << std::endl;
            treeBaseMode = true;
            if (variableMap.count("mean-tracts"))
            {
                std::cout << "WARNING: single tract roi mode, mean tract ids inserted will be ignored"<< inputTractFolder << std::endl;
            }
        }
        else if (variableMap.count("mean-tracts"))
        {
            inNodes= (variableMap["mean-tracts"].as<std::vector<size_t> >());
            std::cout << "flipping " << inNodes.size() << " node mean tracts" << std::endl;
            loopLength = inNodes.size();
        }
        else
        {
            std::cerr << "ERROR: no roi file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("mask")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(maskFilename))) {
                std::cerr << "ERROR: mask file \""<<maskFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Mask file: "<< maskFilename << std::endl;
        } else {
            std::cerr << "ERROR: no mask file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("input-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(inputTractFolder))) {
                std::cerr << "ERROR: folder for input tractograms \""<<inputTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "input tractogram folder: "<< inputTractFolder << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tract folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }

        if (variableMap.count("flipped-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(flippedTractFolder))) {
                std::cerr << "ERROR: folder for single tractograms \""<<flippedTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "flipped tractogram folder: "<< flippedTractFolder << std::endl;
        } else {
            std::cerr << "ERROR: no flipped tract folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }


        std::string logFilename(flippedTractFolder+"/"+progName+"_log.txt");
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
        logFile <<"Mask file:\t"<< maskFilename <<std::endl;
        logFile <<"Single tracts folder:\t"<< inputTractFolder <<std::endl;
        logFile <<"Flipped tracts folder:\t"<< flippedTractFolder <<std::endl;
        logFile <<"-------------"<<std::endl;

        /////////////////////////////////////////////////////////////////

        std::vector<WHcoord> roi;
        std::vector<size_t> bNodes;
        WHcoord datasetSize;

        if (roiMode)
        {
            // ============ Read roi =============

            HC_GRID datasetGrid;

            WFileParser parser( roiFilename );
            if( !parser.readFile() ) {
                std::cerr<<"ERROR @ treeBuilder::readRoi(): Parser error"<<std::endl;
                return false;
            }
            std::vector<std::string>lines = parser.getRawLines();
            if( lines.size() == 0 ) {
                std::cerr<<"ERROR @ treeBuilder::readRoi(): File is empty"<<std::endl;
                return false;
            }

            {
                std::vector< std::vector< std::string> >datasetStrings = parser.getLinesForTagSeparated( "imagesize" );
                if( datasetStrings.size() == 0 ) {
                    std::cerr<<"ERROR @ treeBuilder::readRoi(): Dataset size was not found in tree file"<<std::endl;
                    return false;
                }
                if( datasetStrings.size() > 1 ) {
                    std::cerr<<"ERROR @ treeBuilder::readRoi(): Dataset attribute had multiple lines"<<std::endl;
                    return false;
                }
                WHcoord dataSize(boost::lexical_cast<coord_t>(datasetStrings[0][0]) ,
                                       boost::lexical_cast<coord_t>(datasetStrings[0][1]) ,
                                       boost::lexical_cast<coord_t>(datasetStrings[0][2]) );
                std::string gridString(datasetStrings[0][3]);
                if (gridString==getGridString(HC_VISTA)) {
                    datasetGrid=HC_VISTA;
                } else if (gridString==getGridString(HC_NIFTI)) {
                    std::cerr<<"ERROR @ treeBuilder::readRoi(): "<< gridString <<" format not supported, only "<< getGridString(HC_VISTA) <<" format supported"<<std::endl;
                    return false;
                } else {
                    std::cerr<<"ERROR @ treeBuilder::readRoi(): Dataset grid type string \""<< gridString <<"\" could not be identified"<<std::endl;
                    return false;
                }
                datasetSize = dataSize;
            }
            {
                std::vector< std::vector< std::string> >coordStrings = parser.getLinesForTagSeparated( "roi" );
                if( coordStrings.size() == 0 ) {
                    std::cerr<<"ERROR @ treeBuilder::readRoi(): no roi coordinates in roi file (lacking #roi tag?)"<<std::endl;
                    return false;
                }
                roi.reserve(coordStrings.size());
                for( size_t i = 0; i < coordStrings.size(); ++i ) {
                    WHcoord tempCoord(boost::lexical_cast<coord_t>(coordStrings[i][0]) ,
                                         boost::lexical_cast<coord_t>(coordStrings[i][1]) ,
                                         boost::lexical_cast<coord_t>(coordStrings[i][2]) );
                    roi.push_back(tempCoord);
                }
            }
            std::sort(roi.begin(),roi.end());
            if(verbose) {
                std::cout<<"Roi loaded, "<< roi.size() <<" seed voxels"<<std::endl;
                std::cout<<"Dataset size is: "+datasetSize.getNameString()+" in "+getGridString(datasetGrid)+" format";
            }

            if (datasetGrid != HC_VISTA) {
                    std::cerr<<"Roi coordinates are not in vista format. Only vista is supported for this operation"<<std::endl;
                    return 1;
            }

            loopLength=roi.size();

        }
        else if( treeBaseMode )
        {
            WHtree thisTree( treeFilename );
            bNodes = thisTree.getRootBaseNodes();
            loopLength=bNodes.size();

        }

        // ======= flip and write tracts ==========

        if (verbose)
            std::cout<<std::endl<<"Saving flipped single tracts to folder: "<< flippedTractFolder <<"..." << std::endl;



        time_t lastTime(time(NULL)), startTime(time(NULL));
        size_t tractProg(0), realDone(0);


        #pragma omp parallel for schedule(guided)
        for (size_t i=0; i<loopLength; ++i) {


            vistaManager tractReader(inputTractFolder);
            tractReader.readAsLog();
            tractReader.readAsUnThres();
            tractReader.loadMask(maskFilename);

            vistaManager tractWriter(flippedTractFolder);
            if (roiMode)
            {
                tractWriter.writeInChar();
                tractWriter.storeUnzipped();

            }
            else
            {
                tractWriter.writeInFloat();
                tractWriter.storeZipped();

            }

            std::string tractFilename;
            compactTract tract;


            #pragma omp atomic
                ++tractProg;

            if (roiMode)
            {
                WHcoord flippedCoord(roi[i]);
                flippedCoord.m_x= datasetSize.m_x-1-(flippedCoord.m_x);
                tractWriter.getTractFilename(flippedCoord,tractFilename);



                if( !overWrite )
                {
                    // if tract already exists dont overwrite1
                    if(boost::filesystem::is_regular_file(boost::filesystem::path(tractFilename)))
                        continue;
                }

                tractReader.readLeafTract(roi[i],tract);
                tractReader.flipXtract(tract);
                tractWriter.writeLeafTract(flippedCoord,tract);
            }
            else if(treeBaseMode)
            {
                tractWriter.getTractFilename(bNodes[i],tractFilename);

                if( !overWrite )
                {
                    // if tract already exists dont overwrite1
                    if(boost::filesystem::is_regular_file(boost::filesystem::path(tractFilename)))
                        continue;
                }

                tractReader.readNodeTract(bNodes[i],tract);
                tractReader.flipXtract(tract);
                tractWriter.writeNodeTract(bNodes[i],tract);
            }
            else
            {
                tractWriter.getTractFilename(inNodes[i],tractFilename);

                if( !overWrite )
                {
                    // if tract already exists dont overwrite1
                    if(boost::filesystem::is_regular_file(boost::filesystem::path(tractFilename)))
                        continue;
                }

                tractReader.readNodeTract(inNodes[i],tract);
                tractReader.flipXtract(tract);
                tractWriter.writeNodeTract(inNodes[i],tract);
            }


            #pragma omp atomic
                ++realDone;

            #pragma omp single nowait // only one thread executes output
            if (verbose) {
                time_t currentTime(time(NULL));
                if(currentTime-lastTime>1) {
                    // output progress every 2 seconds
                    lastTime=currentTime;
                    size_t currentCount(tractProg);
                    float progress = (currentCount)*100./(loopLength);
                    int expected_remain(difftime(currentTime,startTime)*((100.-progress)/progress));
                    std::cout<<"\r"<<(int)progress<<" % Completed ("<< currentCount <<" flipped tracts)" <<". Expected remaining time: ";
                    std::cout<<expected_remain/3600 <<"h "<<  (expected_remain%3600)/60 <<"' "<< ((expected_remain%3600)%60) <<"\"  "<< std::flush;
                }
            }

        } // end for

        if (verbose) {
            std::cout<<"\r100 % Completed ("<< tractProg <<" flipped tracts)" << std::endl;
            std::cout<<"( Only "<< realDone <<" computed now)" << std::endl;

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

