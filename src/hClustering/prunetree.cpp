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
#include "WHtreeProcesser.h"

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

        std::string progName("prunetree");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string treeFilename, outputFolder, treeName, basesFilename;
        PROCMODE procMode(PROC_STND);
        float sizePruneRatio(0), distPruneLevel(0), flatGap(0);
        size_t randPruneNumber(0), safeSize(0), joinSize(0);
        unsigned int coarseRatio(0);
        bool keepBases( false);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("verbose,v", "verbose option")
                ("tree-file,t",  boost::program_options::value< std::string >(&treeFilename), "file with the tree file")
                ("tree-name,n",  boost::program_options::value< std::string >(&treeName), "name of the tree (for output naming)")
                ("output,o",  boost::program_options::value< std::string >(&outputFolder), "output folder where processed tree(s) will be written")
                ("standard,s", "standard tree processing, inserted tree must be non-monotonic raw tree")
                ("collapse,c", boost::program_options::value< float >(&flatGap), "perform linear collapse, enter gap value")
                ("keepbases,k", "keep base nodes")
                ("base,b", boost::program_options::value< std::string >(&basesFilename),
                 "homogeneous base tree processing, file with list of bases must follow, inserted tree must be non-monotonic raw tree")
                ("bases2leaves,l", "convert bases to leaves and eliminate coordinate information, inserted tree must be a base tree")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
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

        if (variableMap.count("tree-name"))
        {
            std::cout << "Tree name: "<< treeName << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no tree name stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
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

        if ( ( variableMap.count("standard") + variableMap.count("base") +
               variableMap.count("bases2leaves") + variableMap.count("collapse") ) > 1 )
        {
            std::cerr << "ERROR: multiple processing options chosen"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("standard"))
        {
            std::cout << "Standard tree processing:"<<std::endl;
            procMode=PROC_STND;
        }

        if (variableMap.count("collapse"))
        {
            std::cout << "Linear collapse of nodes:"<<std::endl;
            procMode=PROC_CLPS;
        }


        if (variableMap.count("base"))
        {
            std::cout << "Base tree processing, base file: "<< basesFilename <<std::endl;
            procMode=PROC_BASE;
        }

        if (variableMap.count("bases2leaves"))
        {
            std::cout << "Bases to leaves tree processing: " <<std::endl;
            procMode=PROC_B2L;
        }


        /////////////////////////////////////////////////////////////////


        std::string outFilename;

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


        switch( procMode )
        {

        case PROC_STND:
        {

            WHtreeProcesser treePrcssr(&tree);
            WHtree treeUp(tree);
            WHtreeProcesser treePrcssrUp(&treeUp);
            WHtree treeDown(tree);
            WHtreeProcesser treePrcssrDown(&treeDown);


            std::cout<< "forcing monotonicity up." << std::endl;
            logFile<< "forcing monotonicity up." << std::endl;
            treePrcssrUp.forceMonotonicityUp();
            logFile << treeUp.getReport( false ) <<std::endl;
            std::cout<<treeUp.getReport( false )<<std::endl;
            std::cout<< "Debinarizing." << std::endl;
            treePrcssrUp.debinarize( false );
            logFile << treeUp.getReport( false ) <<std::endl;
            std::cout << treeUp.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + "_monoUp.txt";
            treeUp.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;

            std::cout<< "forcing monotonicity down." << std::endl;
            logFile<< "forcing monotonicity down." << std::endl;
            treePrcssrDown.forceMonotonicityDown();
            logFile << treeDown.getReport( false ) <<std::endl;
            std::cout<<treeDown.getReport( false )<<std::endl;
            std::cout<< "Debinarizing." << std::endl;
            treePrcssrDown.debinarize( false );
            logFile << treeDown.getReport( false ) <<std::endl;
            std::cout << treeDown.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + "_monoDown.txt";
            treeDown.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;

            std::cout<< "forcing monotonicity weighted." << std::endl;
            logFile<< "forcing monotonicity weighted." << std::endl;
            treePrcssr.forceMonotonicity();
            logFile << tree.getReport( false ) <<std::endl;
            std::cout<<tree.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + "_bin.txt";
            tree.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;
            std::cout<< "Debinarizing." << std::endl;
            treePrcssr.debinarize( false );
            logFile << tree.getReport( false ) <<std::endl;
            std::cout << tree.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + ".txt";
            tree.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;
        }
            break;

        case PROC_CLPS:
        {
            if (variableMap.count("keepbases"))
            {
                std::cout << "keeping base nodes"<<std::endl;
                keepBases=true;
            }
            WHtreeProcesser treePrcssr(&tree);
            std::cout<< "Linear node collapse, gap: "<< flatGap << std::endl;
            logFile<< "Linear node collapse, gap: "<< flatGap << std::endl;
            treePrcssr.collapseTreeLinear(0.05 , keepBases);
            logFile << tree.getReport( false ) <<std::endl;
            std::cout << tree.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + ".txt";
            tree.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;

        }
            break;
        case PROC_BASE:
        {

            WFileParser parser( basesFilename );
            std::vector< size_t > baseVector;
            std::vector< size_t > prunedVector;

            if( !parser.readFile() )
            {
                std::cerr << "ERROR: Parser error when reading bases" << std::endl;
                exit(-1);
            }
            std::vector<std::string> lines = parser.getRawLines();
            if( lines.size() == 0 )
            {
                std::cerr << "ERROR: bases file is empty" << std::endl;
                exit(-1);
            }
            {
                std::vector< std::vector< std::string> >baseStrings = parser.getLinesForTagSeparated( "bases" );
                baseVector.reserve( baseStrings.size() );
                for( size_t i = 0; i < baseStrings.size(); ++i )
                {
                    if( baseStrings[i].size() != 1 )
                    {
                        std::cerr << "ERROR: multiple base IDs in the same line, check format" << std::endl;
                        exit(-1);
                    }
                    baseVector.push_back( boost::lexical_cast<size_t>( baseStrings[i][0] ) );
                }

                std::vector< std::vector< std::string> >prunedStrings = parser.getLinesForTagSeparated( "pruned" );
                prunedVector.reserve( prunedStrings.size() );
                for( size_t i = 0; i < prunedStrings.size(); ++i )
                {
                    if( prunedStrings[i].size() != 1 )
                    {
                        std::cerr << "ERROR: multiple pruned IDs in the same line, check format" << std::endl;
                        exit(-1);
                    }
                    prunedVector.push_back( boost::lexical_cast<size_t>( prunedStrings[i][0] ) );
                }
            }

            WHtreeProcesser treePrcssr( &tree );
            treePrcssr.flagLeaves( prunedVector );
            WHtree treeUp(tree);
            WHtreeProcesser treePrcssrUp(&treeUp);
            WHtree treeDown(tree);
            WHtreeProcesser treePrcssrDown(&treeDown);

            std::cout<< "forcing monotonicity up." << std::endl;
            logFile<< "forcing monotonicity up." << std::endl;
            treePrcssrUp.forceMonotonicityUp();
            logFile << treeUp.getReport( false ) <<std::endl;
            std::cout<<treeUp.getReport( false )<<std::endl;
            std::cout<< "Flattening base nodes and pruning out unconnected voxels." << std::endl;
            treePrcssrUp.flattenSelection( baseVector, false );
            logFile << treeUp.getReport( false ) <<std::endl;
            std::cout << treeUp.getReport( false )<<std::endl;
            std::cout<< "Debinarizing." << std::endl;
            treePrcssrUp.debinarize( true );
            logFile << treeUp.getReport( false ) <<std::endl;
            std::cout << treeUp.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + "_monoUp.txt";
            treeUp.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;

            std::cout << std::endl << "forcing monotonicity down." << std::endl;
            logFile<< "forcing monotonicity down." << std::endl;
            treePrcssrDown.forceMonotonicityDown();
            logFile << treeDown.getReport( false ) <<std::endl;
            std::cout<<treeDown.getReport( false )<<std::endl;
            std::cout<< "Flattening base nodes and pruning out unconnected voxels." << std::endl;
            treePrcssrDown.flattenSelection( baseVector, false );
            logFile << treeDown.getReport( false ) <<std::endl;
            std::cout << treeDown.getReport( false )<<std::endl;
            std::cout<< "Debinarizing." << std::endl;
            treePrcssrDown.debinarize( true );
            logFile << treeDown.getReport( false ) <<std::endl;
            std::cout << treeDown.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + "_monoDown.txt";
            treeDown.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;

            std::cout << std::endl << "forcing monotonicity weighted." << std::endl;
            logFile<< "forcing monotonicity weighted." << std::endl;
            treePrcssr.forceMonotonicity();
            logFile << tree.getReport( false ) <<std::endl;
            std::cout<<tree.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + "_bin.txt";
            tree.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;
            std::cout<< "Flattening base nodes and pruning out unconnected voxels." << std::endl;
            treePrcssr.flattenSelection( baseVector, false );
            logFile << tree.getReport( false ) <<std::endl;
            std::cout << tree.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + "_bases";
            tree.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;
            std::cout<< "Debinarizing." << std::endl;
            treePrcssr.debinarize( true );
            logFile << tree.getReport( false ) <<std::endl;
            std::cout << tree.getReport( false )<<std::endl;
            outFilename = outputFolder + "/" + treeName + ".txt";
            tree.writeTree(outFilename);
            std::cout << "written to: "<< outFilename <<std::endl;
            logFile << "written to: "<< outFilename <<std::endl;


            if( tree.testRootBaseNodes() )
            {
                std::vector< size_t> newBaseVector = tree.getRootBaseNodes();
                std::sort( newBaseVector.begin(), newBaseVector.end() );
                outFilename = outputFolder + "/baselist.txt";
                std::ofstream newBaseFile( outFilename.c_str() );
                if( !newBaseFile )
                {
                    std::cerr << "ERROR: unable to open out file: \"" << newBaseFile << "\"" << std::endl;
                    exit( -1 );
                }
                newBaseFile << "#bases" << std::endl;
                for( std::vector<size_t>::const_iterator nodeIter( newBaseVector.begin() ) ; nodeIter != newBaseVector.end() ; ++nodeIter )
                {
                    newBaseFile << *nodeIter << std::endl;
                }
                newBaseFile << "#endbases" << std::endl << std::endl;

                std::cout << "Final base list written in: "<< outFilename << std::endl;
                logFile << "Final base list written in: "<< outFilename << std::endl;
            }
            else
            {
                std::cout << "Final tree is not a pure basenode tree" << std::endl;
                logFile << "Final tree is not a pure basenode tree" << std::endl;
            }


        }
            break;

        case PROC_B2L:
        {
            WHtreeProcesser treePrcssr(&tree);
            std::cout << "Collapsing base nodes into leaves..." << std::endl;
            logFile << "Collapsing base nodes into leaves" << std::endl;
            treePrcssr.baseNodes2Leaves();
            std::cout<<"Done. "<<tree.getReport()<<std::endl;
            logFile <<tree.getReport() <<std::endl;
            std::cout << tree.getReport() <<std::endl;
            outFilename = outputFolder + "/" + treeName + "_baset.txt";
            tree.writeTreeSimple(outFilename);
        }
            break;

        default:
            break;
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


#if 0

// Declare a group of options that will be allowed only on command line
boost::program_options::options_description genericOptions("Generic options");
genericOptions.add_options()
        ("version,V", "print version string")
        ("help,h", "produce help message")
        ("monotonic,m", "force tree monotonicity")
        ("tree-file,t",  boost::program_options::value< std::string >(&treeFilename), "file with the tree file")
        ("prune-size,s", boost::program_options::value< size_t >(&joinSize), "prune branches smaller than safe-size that join a node equal or bigger than N")
        ("prune-ratio,r", boost::program_options::value< float >(&sizePruneRatio), "prune branches smaller than safe-size whose size ratio to its parent node is smaller than N")
        ("prune-level,l", boost::program_options::value< float >(&distPruneLevel), "prune branches smaller than safe-size whose parent its at a distance level grater than N")
        ("prune-random,r", boost::program_options::value< size_t >(&randPruneNumber), "prune N randomly selected leaves")
        ("flatten,f", boost::program_options::value< float >(&flatGap), "flatten all those nodes whose distance to their respective parents/children are smaller than N")
        ("coarse,c", boost::program_options::value< unsigned int >(&coarseRatio), "decimate tree results to fit an N-times coarser anatomical grid")
        ("smoothen,h", "performs tree smoothing, requires prune-size, prune-level and safe-size parameters, other options will be ignored")
        ("base-tree,b", "converts base nodes into leaves and writes the simplififed tree")
        ("verbose,v", "verbose option")
        ;


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
    std::cout << "Roi voxels file: "<< treeFilename << std::endl;
} else {
    std::cerr << "ERROR: no tree file stated"<<std::endl;
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

/////////////////////////////////////////////////////////////////


std::string outFile;

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


WHtreeProcesser treePrcssr(&tree);

if (variableMap.count("monotonic") )
{
    WHtree treeUp(tree);
    WHtreeProcesser treePrcssrUp(&treeUp);
    std::cout<< "forcing monotonicity up." << std::endl;
    logFile<< "forcing monotonicity up." << std::endl;
    treePrcssrUp.forceMonotonicityUp();
    logFile << treeUp.getReport( false ) <<std::endl;
    std::cout<<treeUp.getReport( false )<<std::endl;
    outFile = outputFolder + "/" + treeUp.getName() + "_monoUp.txt";
    treeUp.writeTree(outFile);
    std::cout << "written to: "<< outFile <<std::endl;
    logFile << "written to: "<< outFile <<std::endl;

    WHtree treeDown(tree);
    WHtreeProcesser treePrcssrDown(&treeDown);
    std::cout<< "forcing monotonicity down." << std::endl;
    logFile<< "forcing monotonicity down." << std::endl;
    treePrcssrDown.forceMonotonicityDown();
    logFile << treeDown.getReport( false ) <<std::endl;
    std::cout<<treeDown.getReport( false )<<std::endl;
    outFile = outputFolder + "/" + treeDown.getName() + "_monoDown.txt";
    treeDown.writeTree(outFile);
    std::cout << "written to: "<< outFile <<std::endl;
    logFile << "written to: "<< outFile <<std::endl;

    std::cout<< "forcing monotonicity weighted." << std::endl;
    logFile<< "forcing monotonicity weighted." << std::endl;
    treePrcssr.forceMonotonicity();
    logFile << tree.getReport( false ) <<std::endl;
    std::cout<<tree.getReport( false )<<std::endl;
    outFile = outputFolder + "/" + tree.getName() + "_monoW.txt";
    tree.writeTree(outFile);
    std::cout << "written to: "<< outFile <<std::endl;
    logFile << "written to: "<< outFile <<std::endl;
}


if (variableMap.count("smoothen") )
{
    if (variableMap.count("safe-size") && variableMap.count("prune-level") && variableMap.count("prune-size")) {
        std::cout << "Smoothing... safe-size: "<<safeSize<<". smooth-size: "<<joinSize<<". smooth-gap: "<<distPruneLevel<<std::endl;
        logFile <<"Smoothing => safe-size: "<<safeSize<<". smooth-size: "<<joinSize<<". smooth-gap: "<<distPruneLevel<<std::endl;

        if (variableMap.count("prune-ratio") || variableMap.count("prune-random") || variableMap.count("flatten")) {
            std::cerr << "WARNING: prune-ratio, prune-random and flatten options will be ignored"<<std::endl;
        }

        treePrcssr.smoothTree(safeSize,joinSize,distPruneLevel);

        std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
        logFile <<tree.getReport(false) <<std::endl;
        outFile = outputFolder + "/" + tree.getName() + "_smooth.txt";
        tree.writeTree(outFile);
        outFile = outputFolder + "/" + tree.getName() + "_smooth_debug.txt";
        tree.writeTreeDebug(outFile);


    }
    else
    {
        std::cerr << "ERROR: smoothing requires the additional parameters: safe-size, prune-level and prune-size"<<std::endl;
        exit(-1);
    }
}
else
{


    if (variableMap.count("prune-ratio") || variableMap.count("prune-level") || variableMap.count("prune-size"))
    {
        if (variableMap.count("prune-random"))
        {
            std::cerr << "ERROR: pruned random may not be used simultaneously with other pruning methods"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        std::string safeString;
        if(variableMap.count("safe-size"))
        {
            safeString = ("smaller than "+boost::lexical_cast<std::string>(safeSize));
            logFile <<"Pruning-safe size:\t"<< safeSize <<std::endl;

            if (variableMap.count("prune-size"))
            {
                std::cout << "Pruning by minimum size if a cluster smaller than "<< safeSize <<" joins a node bigger than more than "<< joinSize <<". it will be pruned"<< std::endl;
                logFile <<"Min-size pruning:\t"<< safeSize << ":"<< joinSize <<std::endl;

                treePrcssr.pruneTree(joinSize,safeSize,HTPR_JOINSIZE);

                std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
                logFile <<tree.getReport(false) <<std::endl;
                outFile = outputFolder + "/" + tree.getName() + "_prunedS.txt";
                tree.writeTree(outFile);
                outFile = outputFolder + "/" + tree.getName() + "_prunedS_debug.txt";
                tree.writeTreeDebug(outFile);

            }
        }


        if (variableMap.count("prune-ratio"))
        {
            std::cout << "Pruning by size ratio: if a cluster "<< safeString <<" joins a node more than "<< sizePruneRatio <<" times its size, it will be pruned"<< std::endl;
            logFile <<"Size-pruning ratio:\t"<< sizePruneRatio <<std::endl;

            treePrcssr.pruneTree(sizePruneRatio,safeSize,HTPR_SIZERATIO);

            std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
            logFile <<tree.getReport(false) <<std::endl;
            outFile = outputFolder + "/" + tree.getName() + "_prunedR.txt";
            tree.writeTree(outFile);
            outFile = outputFolder + "/" + tree.getName() + "_prunedR_debug.txt";
            tree.writeTreeDebug(outFile);

        }
        else if (variableMap.count("prune-level"))
        {
            if (safeSize==0)
            {
                std::cerr << "ERROR: when pruning by level a positive safe size must be stated"<<std::endl;
                exit(-1);
            }
            std::cout << "Pruning by distance level: if a cluster "<< safeString <<" joins a node at a level higher than "<< distPruneLevel <<" times its size, it will be pruned"<< std::endl;
            logFile <<"Level-pruning distance:\t"<< distPruneLevel <<std::endl;

            treePrcssr.pruneTree(distPruneLevel,safeSize,HTPR_JOINLEVEL);

            std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
            logFile <<tree.getReport(false) <<std::endl;
            outFile = outputFolder + "/" + tree.getName() + "_prunedL.txt";
            tree.writeTree(outFile);
            outFile = outputFolder + "/" + tree.getName() + "_prunedL_debug.txt";
            tree.writeTreeDebug(outFile);
        }

    }

    if (variableMap.count("prune-random"))
    {
        std::cout << "Randomly pruning "<< randPruneNumber <<" leaves from the tree"<< std::endl;
        logFile <<"Randomly pruned leaves:\t"<< randPruneNumber <<std::endl;


        std::cout << tree.getReport(false) <<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        treePrcssr.pruneRandom(randPruneNumber);

        std::cout<<"Done. "<<tree.getReport(false)<<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_prunedRand.txt";
        tree.writeTree(outFile);

        std::cout << "written to: "<< outFile <<std::endl;
        logFile << "written to: "<< outFile <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_prunedRand_debug.txt";
        tree.writeTreeDebug(outFile);




    }

    if (variableMap.count("flatten"))
    {
        std::cout << "Flattening branches with intervals smaller than "<< flatGap << std::endl;
        logFile <<"Flattening gap:\t"<< flatGap <<std::endl;

        std::cout << tree.getReport(false) <<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        treePrcssr.collapseTree(flatGap,1);

        std::cout << tree.getReport(false) <<std::endl;
        logFile <<tree.getReport(false) <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_flat.txt";
        tree.writeTree(outFile);

        std::cout << "written to: "<< outFile <<std::endl;
        logFile << "written to: "<< outFile <<std::endl;

        outFile = outputFolder + "/" + tree.getName() + "_flat_debug.txt";
        tree.writeTreeDebug(outFile);


    }

}

if (variableMap.count("coarse"))
{
    std::cout << "Decimating tree to a grid "<< coarseRatio <<" times coarser" << std::endl;
    logFile <<"Coarsing factor:\t"<< coarseRatio <<std::endl;

    treePrcssr.coarseTree(coarseRatio);

    std::cout<<"Done. "<<tree.getReport()<<std::endl;
    logFile <<tree.getReport() <<std::endl;
    outFile = outputFolder + "/" + tree.getName() + "_coarse.txt";
    tree.writeTree(outFile);
    outFile = outputFolder + "/" + tree.getName() + "_coarse_debug.txt";
    tree.writeTreeDebug(outFile);

}

if (variableMap.count("base-tree"))
{
    logFile <<tree.getReport() <<std::endl;
    std::cout << tree.getReport() <<std::endl;

    std::cout << "Collapsing base nodes into leaves..." << std::endl;
    logFile << "Collapsing base nodes into leaves" << std::endl;
    treePrcssr.baseNodes2Leaves();

    std::cout<<"Done. "<<tree.getReport()<<std::endl;

    logFile <<tree.getReport() <<std::endl;
    std::cout << tree.getReport() <<std::endl;

    outFile = outputFolder + "/" + tree.getName() + "_baset.txt";
    tree.writeTreeSimple(outFile);
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



#endif

