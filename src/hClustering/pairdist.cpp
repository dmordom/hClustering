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
#include "WHcoord.h"
#include "WHnode.h"
#include "WHtree.h"
#include "compactTract.h"
#include "compactTractChar.h"
#include "vistaManager.h"
#include "distBlock.h"


const float thresValue(0.4);


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
        //boost::filesystem::path workingDir( boost::filesystem::current_path());

        // ========== PROGRAM PARAMETERS ==========

        std::string progName("pairdist");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string treeFilename;
        std::string singleTractFolder, meanTractFolder;
        std::string distMatrixFolder;
        std::pair<nodeID_t,nodeID_t> inNodes;
        std::pair<WHcoord,WHcoord> inCoords;
        bool byTract(false), byMatrix(false), byTree(false), byCoords(false), byNodes(false), verbose(false);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("input-nodes,n", boost::program_options::value<bool> (), "input nodeID pair (each consisting of a boolean leaf/node flag value followed by an id number) to compute distance of, if used, input-coords option will be ignored")
                ("input-coords,c", boost::program_options::value<coord_t> (), "input leaf coordinate pair to compute distance of")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("tree-file,t",  boost::program_options::value< std::string >(&treeFilename), "file with the hierarchical tree")
                ("dist-folder,d",  boost::program_options::value< std::string >(&distMatrixFolder), "folder with the distance matrix files")
                ("singlet-folder,s", boost::program_options::value< std::string >(&singleTractFolder), "folder with the single-voxel probabilistic tracts")
                ("meant-folder,f", boost::program_options::value< std::string >(&meanTractFolder), "folder with the node mean tracts")
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


        ////////////////////////////
# if 0
        if (variableMap.count("input-coords")) {
            coord_t mode(variableMap["input-coords"].as<coord_t>());
            WHcoord coor(10,10,10);
            WHcoord sizeG(160,200,160);
            std::cout<<mode<<" nbhood of coordinate "<<coor<<std::endl;
            std::vector<WHcoord> nbs(coor.getPhysNbs(sizeG,mode));
            std::cout<<"Total nbs: "<<nbs.size()<<std::endl;
            for (size_t i=0; i<nbs.size();++i) {
                std::cout<<nbs[i]<<std::endl;
            }
            return 0;
        }

#endif
        ////////////////////////////////////





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
        if (variableMap.count("tree-file")) {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(treeFilename))) {
                std::cerr << "ERROR: tree file \""<<treeFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Input tree file: "<< treeFilename << std::endl;
            byTree=true;
        } else if (variableMap.count("input-nodes")) {
            std::cout << "WARNING: input-nodes option will be ignored as no tree file was stated" << std::endl;
        }
        if (variableMap.count("dist-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(distMatrixFolder))) {
                std::cerr << "ERROR: distance matrix folder \""<<distMatrixFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Distance matrix folder: "<< distMatrixFolder << std::endl;
            byMatrix=true;
        }
        if (variableMap.count("singlet-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(singleTractFolder))) {
                std::cerr << "ERROR: single tract folder \""<<singleTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Single tracts folder: "<< singleTractFolder << std::endl;
            byTract=true;
        }
        if (variableMap.count("meant-folder")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(meanTractFolder))) {
                std::cerr << "ERROR: mean tract folder \""<<meanTractFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Mean tracts folder: "<< meanTractFolder << std::endl;
            byTract=true;
        }
        if(!byTree && !byMatrix && !byTract) {
            std::cerr << "ERROR: at least one source must be entered (tract-folder, mdist-folder, or tree-file) "<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("input-nodes") && byTree) {

            if (!variableMap.count("rest")) {
                std::cerr << "ERROR: 2 full node IDs must be entered"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::vector<size_t> inNodesVector(variableMap["rest"].as<std::vector<size_t> >());

            if(inNodesVector.size()!=3) {
                std::cerr << "ERROR: 2 full node IDs must be entered"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            } else {

                inNodes=std::make_pair(std::make_pair(variableMap["input-nodes"].as<bool>(),inNodesVector[0]), std::make_pair(inNodesVector[1],inNodesVector[2]));
                std::cout << "Input nodes are: "<< inNodes.first.first <<"-"<< inNodes.first.second <<" and "<< inNodes.second.first <<"-"<< inNodes.second.second << std::endl;
                byNodes=true;

                if  (inNodes.first.first||inNodes.second.first) {
                    if ((byTract)&&(!variableMap.count("meant-folder")) ) {
                        std::cerr << "ERROR: When using node ids, mean tract folder must be stated"<<std::endl;
                        std::cerr << visibleOptions << std::endl;
                        exit(-1);
                    }
                    if (byMatrix) {
                        std::cout << "WARNING: dist-folder option will be ignored" << std::endl;
                        byMatrix=false;
                    }
                }
                if ( (!inNodes.first.first||!inNodes.second.first)&&(byTract)&&(!variableMap.count("singlet-folder")) ) {
                    std::cerr << "ERROR: When using leaf ids, single tract folder must be stated"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }

                if (byTree && variableMap.count("input-coords"))
                    std::cout << "WARNING: input-coords option will be ignored" << std::endl;
            }

        } else if (variableMap.count("input-coords")){

            if (!variableMap.count("rest")) {
                std::cerr << "ERROR: 2 full coordinate sets must be entered"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::vector<coord_t> inCoordVector;
            inCoordVector.push_back(variableMap["input-coords"].as<coord_t>());
            {
                std::vector<size_t> tempVector(variableMap["rest"].as<std::vector<size_t> >());
                inCoordVector.reserve(tempVector.size());
                for(std::vector<size_t>::iterator iter(tempVector.begin()); iter!=tempVector.end(); ++iter) {
                    inCoordVector.push_back(*iter);
                }
            }


            if(inCoordVector.size()!=6) {
                std::cerr << "ERROR: a total of 6 coordinates must be entered (3 per leaf)"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            } else {
                WHcoord coord1(inCoordVector[0],inCoordVector[1],inCoordVector[2]);
                WHcoord coord2(inCoordVector[3],inCoordVector[4],inCoordVector[5]);
                inCoords=std::make_pair(coord1,coord2);
                std::cout << "Input coordinates are: "<< coord1 <<" and "<< coord2 << std::endl;
                byCoords=true;

            }
        } else {
            std::cerr << "ERROR: no input nodes stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        // ========== OBTAIN DISTANCES ==========



        if(byTree) {
            WHtree tree(treeFilename);
            std::cout<<tree.getReport()<<std::endl<<std::endl;
            dist_t copheneticDist(0);
            if(byNodes) {
                if (inNodes.first.first) {
                    std::cout<<"Node "<< inNodes.first.first <<"-"<< inNodes.first.second <<std::endl;
                } else {
                    inCoords.first=tree.getCoordinate4leaf(inNodes.first.second);
                    std::cout<<"Leaf "<< inNodes.first.first <<"-"<< inNodes.first.second <<" : "<<inCoords.first<<std::endl;
                }
                if (inNodes.second.first) {
                    std::cout<<"Node "<< inNodes.second.first <<"-"<< inNodes.second.second <<std::endl;
                } else {
                    inCoords.second=tree.getCoordinate4leaf(inNodes.second.second);
                    std::cout<<"Leaf "<< inNodes.second.first <<"-"<< inNodes.second.second <<" : "<<inCoords.second<<std::endl;
                }
                copheneticDist =  tree.getDistance(inNodes.first,inNodes.second);
            } else {
                copheneticDist = tree.getDistance(inCoords.first, inCoords.second);
            }
            std::cout<<"Cophenetic distance:\t"<<boost::lexical_cast<std::string>(copheneticDist)<<std::endl;
        }

        if(byTract) {
            compactTract tract1, tract2;
            compactTractChar charTract1, charTract2;
            if(byNodes) {
                if (inNodes.first.first) {
                    vistaManager vMngr(meanTractFolder);
                    vMngr.readAsLog();
                    vMngr.readAsUnThres();
                    vMngr.readNodeTract(inNodes.first.second,tract1);                    
                } else {
                    vistaManager vMngr(singleTractFolder);
                    vMngr.readAsLog();
                    vMngr.readAsUnThres();
                    vMngr.readLeafTract(inCoords.first,tract1);
                    vMngr.readLeafTract(inCoords.first,charTract1);
                }
                if (inNodes.second.first) {
                    vistaManager vMngr(meanTractFolder);
                    vMngr.readAsLog();
                    vMngr.readAsUnThres();
                    vMngr.readNodeTract(inNodes.second.second,tract2);
                } else {
                    vistaManager vMngr(singleTractFolder);
                    vMngr.readAsLog();
                    vMngr.readAsUnThres();
                    vMngr.readLeafTract(inCoords.second,tract2);
                    vMngr.readLeafTract(inCoords.second,charTract2);
                }
            } else {
                vistaManager vMngr(singleTractFolder);
                vMngr.readAsLog();
                vMngr.readAsUnThres();
                vMngr.readLeafTract(inCoords.first,tract1);
                vMngr.readLeafTract(inCoords.first,charTract1);
                vMngr.readLeafTract(inCoords.second,tract2);
                vMngr.readLeafTract(inCoords.second,charTract2);
            }
            dist_t directDist(0);
            tract1.threshold(thresValue);
            tract2.threshold(thresValue);
            charTract1.threshold(thresValue);
            charTract2.threshold(thresValue);
            tract1.getNorm();
            tract2.getNorm();
            charTract1.getNorm();
            charTract2.getNorm();

            directDist = tract1.tractDistance(tract2);
            std::cout<<"Direct distance:\t"<<boost::lexical_cast<std::string>(directDist)<<std::endl;

            if( charTract1.size() > 0 )
            {
                dist_t charDist(0);
                charDist = charTract1.tractDistance(tract2);
                std::cout<<"Direct (char-float) distance:\t"<<boost::lexical_cast<std::string>(charDist)<<std::endl;
            }
            if( charTract2.size() > 0 )
            {
                dist_t charDist(0);
                charDist = tract1.tractDistance(charTract2);
                std::cout<<"Direct (float-char) distance:\t"<<boost::lexical_cast<std::string>(charDist)<<std::endl;
            }

            if( charTract1.size() > 0 && charTract2.size() > 0 )
            {
                dist_t charDist(0);
                charDist = charTract1.tractDistance(charTract2);
                std::cout<<"Direct (char-char) distance:\t"<<boost::lexical_cast<std::string>(charDist)<<std::endl;
            }
        }

        if(byMatrix) {
            dist_t matrixDist(0);
            distBlock distanceBlock(distMatrixFolder);
            distanceBlock.loadBlock(inCoords.first,inCoords.second);
            matrixDist = distanceBlock.getDistance(inCoords.first,inCoords.second);
            std::cout<<"Matrix distance:\t"<<boost::lexical_cast<std::string>(matrixDist)<<std::endl;
        }
        std::cout<<std::endl;


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

