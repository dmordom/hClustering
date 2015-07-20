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
#include "treeComparer.h"

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

        std::string progName("comparetrees");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string treeFilename1, treeFilename2, meanTractFolder1, meanTractFolder2, singleTractFolder1, singleTractFolder2, outputFolder, matrixFilename;
        unsigned int threads(0);
        bool meanTractsFromFile( false ), coordsFromFile( false ), nocomp( false ), onlyCpcc( false );
        CRSP_MODE compMode(CRSP_DIRECT);
        bool matchNoise( false );
//        float noiseAlpha( 0 );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version,V", "print version string")
                ("help,h", "produce help message")
                ("cl", "direct leaf-wise correspondence")
                ("cs", boost::program_options::value< std::string >(&matrixFilename), "simple base node correspondence, file with the similarity matrix/file where to write it")
                ("cm", "merge-up base node correspondence")
                ("cr", "random base node correspondence")
                ("t1",  boost::program_options::value< std::string >(&treeFilename1), "file with first tree")
                ("t2",  boost::program_options::value< std::string >(&treeFilename2), "file with second tree")
                ("m1",  boost::program_options::value< std::string >(&meanTractFolder1), "folder with the mean tracts for the first tree")
                ("m2",  boost::program_options::value< std::string >(&meanTractFolder2), "folder with the mean tracts for the second tree")
                ("output,o",  boost::program_options::value< std::string >(&outputFolder), "output folder where results will be written")
                ("nocomp", "only obtain tree correspondence, not the ")
                ("cpcc", "only obtain the cpcc, not the triples")
                ("verbose,v", "verbose option")
//                ("noise,n", boost::program_options::value< float >(&noiseAlpha), "matching noise correction. insert alpha value (0-1)")
                ("noise,n", "loop alpha values of matching noise correction")

                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("s1",  boost::program_options::value< std::string >(&singleTractFolder1), "folder with the single voxel tractograms for first tree")
                ("s2",  boost::program_options::value< std::string >(&singleTractFolder2), "folder with the single voxel tractograms for second tree")
                ("threads,p",  boost::program_options::value< unsigned int >(&threads), "fnumber of processing threads to run the program in parallel, default: all available")
                ("inter-subject,i", "comparison is between different brains (mean tracts and coordinates taken from file)")
                ("tract-files,f", "mean tracts will be taken from file")
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




        if (variableMap.count("t1")) {
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

        if (variableMap.count("t2")) {
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


        if (variableMap.count("m1")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(meanTractFolder1))) {
                std::cerr << "ERROR: Mean tract folder for tree1 \""<<meanTractFolder1<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "Mean tract folder for tree1: "<< meanTractFolder1 << std::endl;
        } else {
            std::cerr << "ERROR: no mean tract folder for tree1 stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }

        if (variableMap.count("m2")) {
            if(!boost::filesystem::is_directory(boost::filesystem::path(meanTractFolder2))) {
                std::cerr << "ERROR: Mean tract folder for tree2 \""<<meanTractFolder2<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            std::cout << "Mean tract folder for tree2: "<< meanTractFolder2 << std::endl;
        } else {
            std::cerr << "ERROR: no mean tract folder for tree2 stated"<<std::endl;
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

        if (variableMap.count("inter-subject"))
        {
            std::cout << "Inter subject study, tracts and coords are taken from file"<<std::endl;
            meanTractsFromFile = true;
            coordsFromFile = true;
        }
        else if (variableMap.count("tract-files"))
        {
            std::cout << "Tracts and coords are taken from file"<<std::endl;
            meanTractsFromFile = true;
        }
        else
        {
            if (variableMap.count("s1")) {
                if(!boost::filesystem::is_directory(boost::filesystem::path(singleTractFolder1))) {
                    std::cerr << "ERROR: folder for single tractograms for tree 1 \""<<singleTractFolder1<<"\" is not a directory"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }
                std::cout << "Single tractogram folder for tree 1: "<< singleTractFolder1 << std::endl;
            } else {
                std::cerr << "ERROR: no single tract folder stated for tree 1"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if (variableMap.count("s2")) {
                if(!boost::filesystem::is_directory(boost::filesystem::path(singleTractFolder2))) {
                    std::cerr << "ERROR: folder for single tractograms for tree 2 \""<<singleTractFolder2<<"\" is not a directory"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }
                std::cout << "Single tractogram folder for tree 2: "<< singleTractFolder2 << std::endl;
            } else {
                std::cerr << "ERROR: no single tract folder stated for tree 2"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
        }

        if (variableMap.count("nocomp"))
        {
            std::cout << "Obtaining only tree correspondence"<<std::endl;
            nocomp = true;
        }
        else if (variableMap.count("cpcc"))
        {
            std::cout << "Obtaining only cpcc and not triples"<<std::endl;
            onlyCpcc = true;
        }


        if (variableMap.count("noise"))
        {
            std::cout << "Using matching noise correction. ";
            matchNoise = true;
//            if( noiseAlpha < 0)
//            {
//                std::cout << "alpha value insterted is invalid (negative), not using matching noise correction";
//                noiseAlpha = 0;
//            }
//            else if (noiseAlpha > 1)
//            {
//                std::cout << "alpha value insterted is too high ( >1 ), using an alpha of 1";
//                noiseAlpha = 1;
//            }
//            else
//            {
//                std::cout << "matching noise alpha value: "<<noiseAlpha;
//            }
//            std::cout<<std::endl;
        }
        else
        {
            std::cout << "NOT Using matching noise correction. ";
            matchNoise = false;
        }


        size_t countComp(0);
        if (variableMap.count("cl")) {
            std::cout << "Direct leaf-wise correspondence"<<std::endl;
            compMode=CRSP_DIRECT;
            ++countComp;
        }
        if (variableMap.count("cs")) {
            std::cout << "Simple baseNode-wise correspondence"<<std::endl;
            compMode=CRSP_SIMPLE;
            ++countComp;
        }
        if (variableMap.count("cm")) {
            std::cout << "Merge-Up baseNode-wise correspondence"<<std::endl;
            compMode=CRSP_MERGE;
            ++countComp;
        }
        if (variableMap.count("cr")) {
            std::cout << "Random baseNode-wise correspondence"<<std::endl;
            compMode=CRSP_RAND;
            ++countComp;
        }

        if (countComp==0) {
            std::cerr << "ERROR: no comparison mode stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if (countComp!=1) {
            std::cerr << "ERROR: More than one comparison mode stated, only one mode allowed"<<std::endl;
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
        logFile <<"Tree1 file:\t"<< treeFilename1 <<std::endl;
        logFile <<"Tree2 file:\t"<< treeFilename2 <<std::endl;
        logFile <<"Mean tracts folder for tree1:\t"<< meanTractFolder1 <<std::endl;
        logFile <<"Mean tracts folder for tree2:\t"<< meanTractFolder2 <<std::endl;
        logFile <<"Single tracts folder for tree 1:\t"<< singleTractFolder1 <<std::endl;
        logFile <<"Single tracts folder for tree 2:\t"<< singleTractFolder2 <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Tracts from file:\t"<< meanTractsFromFile <<std::endl;
        logFile <<"Coordinates from file:\t"<< coordsFromFile <<std::endl;
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
            throw std::runtime_error ("ERROR @ compareTrees() datasets to compare have different dimensions");
        }

        if (verbose) {
            std::cout <<"Tree1: "<< tree1.getReport(false) <<std::endl;
            std::cout <<"Tree2: "<< tree2.getReport(false) <<std::endl;
        }


        treeComparer comparer(&tree1,&tree2, verbose);
        comparer.setCoordsFromFile( coordsFromFile );
        comparer.setMeanTractsFromFile( meanTractsFromFile );

#if 0

        std::cout <<comparer.reportBaseNodes()<<std::endl;

        size_t numred(100);
        std::cout <<"Reducing to "<< numred <<" base nodes..."<<std::endl;

        comparer.reduceBaseNodes(numred);
        if (verbose) {
            std::cout <<"Tree1: "<< tree1.getReport(false) <<std::endl;
            std::cout <<"Tree2: "<< tree2.getReport(false) <<std::endl;
        }
        std::cout <<comparer.reportBaseNodes()<<std::endl;


#endif


        if (compMode!=CRSP_DIRECT) {
            std::cout <<comparer.reportBaseNodes()<<std::endl;
            if (!comparer.areRealBaseNodes()) {
                throw std::runtime_error ("ERROR @ compareTrees(): base nodes are of mixed type (contain both leaves and nodes)");
            }
        }

        comparer.setSingleTractFolder1(singleTractFolder1);
        comparer.setMeanTractFolder1(meanTractFolder1);
        comparer.setSingleTractFolder2(singleTractFolder2);
        comparer.setMeanTractFolder2(meanTractFolder2);
        comparer.log(&logFile);
        comparer.setMaxPhysDist(20);


        if (compMode == CRSP_RAND) {
            logFile <<"Correspondance mode:\t Random"<<std::endl;

            const size_t randRepetitions(100);

            std::cout<<std::endl<<"100 rep loop: "<<std::endl;
            std::string outCpctFilename(outputFolder+"/randCpct.txt");
            std::string outStriplesFilename(outputFolder+"/randStriples.txt");
            std::ofstream outCpctFile(outCpctFilename.c_str());
            if(!outCpctFile) {
                std::cerr << "ERROR: unable to open out file: \""<<outCpctFilename<<"\""<<std::endl;
                exit(-1);
            }
            std::ofstream outStriplestFile(outStriplesFilename.c_str());
            if(!outStriplestFile) {
                std::cerr << "ERROR: unable to open out file: \""<<outStriplesFilename<<"\""<<std::endl;
                exit(-1);
            }

            outCpctFile <<  "weighted-tCPCC    simple-sCPCC   %usedPairs   effectiveGranularity" << std::endl;

            for(int i=0; i<randRepetitions; ++i) {
                comparer.randomCorrespondence();
                std::pair< std::pair< float, float >, std::pair< float, float > >  tCPCC(comparer.doTcpcc());
                outCpctFile <<  tCPCC.first.first << " " <<  tCPCC.first.second << " " << std::flush;
                outCpctFile <<  tCPCC.second.first*100 << " " <<  tCPCC.second.second << std::endl;


                if( !onlyCpcc )
                {
                    std::pair< float, float >sTriplets(comparer.simpleTriplets(3));
                    outStriplestFile << sTriplets.first << " " << sTriplets.second << std::endl;
                }
            }

        } else {

            size_t sampleFreq(1);

            bool redoCoords( true );
            std::string outputFileName(outputFolder+"/compValues.txt");
            std::ofstream outFile(outputFileName.c_str());
            if(!outFile) {
                std::cerr << "ERROR: unable to open out file: \""<<outputFileName<<"\""<<std::endl;
                exit(-1);
            }

            std::vector<float> correspRates;

            switch (compMode) {
            case CRSP_DIRECT:
                logFile <<"Correspondance mode:\t Direct"<<std::endl;
                if (verbose) {
                    std::cout<<"Equalizing leaves.."<<std::endl;
                }
                if(!comparer.leafCorrespondence()) {
                    std::cout<<"Tree coordinates match, no changes made."<<std::endl;
                }

                sampleFreq = 20;
                break;

            case CRSP_SIMPLE:
                logFile <<"Correspondance mode:\t Simple"<<std::endl;
                if( boost::filesystem::is_regular_file( boost::filesystem::path( matrixFilename ) ) )
                {
                    if (verbose)
                    {
                        std::cout <<"Getting zipped similarity matrix from file: "<< matrixFilename << " ..." << std::flush;
                    }
                    comparer.readBaseDistMatrix( matrixFilename );
                    if (verbose)
                    {
                        std::cout <<"Done"<< std::endl;
                    }
                    logFile <<"Similarity matrix read from file: " << matrixFilename << std::endl;

                }
                else if ( boost::filesystem::is_regular_file( boost::filesystem::path( matrixFilename + ".gz" ) ) )
                {
                    if (verbose)
                    {
                        std::cout <<"Getting zipped similarity matrix from file: "<<  matrixFilename + ".gz" << " ..." << std::flush;
                    }
                    comparer.readBaseDistMatrix( matrixFilename + ".gz" );
                    if (verbose)
                    {
                        std::cout <<"Done"<< std::endl;
                    }
                    logFile <<"Similarity matrix read from file: " << matrixFilename << ".gz" << std::endl;
                }
                else
                {
                    if (verbose)
                    {
                        std::cout <<"Similarity matrix file does not exist, it will be computed and saved in file"<< std::endl;
                    }
                    comparer.getBaseDistMatrix();
                    comparer.writeBaseDistMatrix( matrixFilename );
                    if (verbose)
                    {
                        std::cout << "Similarity matrix saved in file: " << matrixFilename << std::endl;
                    }
                    logFile <<"Similarity matrix saved in file: " << matrixFilename << std::endl;

                }
                if (nocomp)
                {
                    redoCoords = false;
                }


                comparer.simpleCorrespondence( 0.9, redoCoords );
                comparer.writeFullCorrespondence( outputFolder + "/fullCorresp.txt" );
                comparer.writeCorrespondence( outputFolder + "/correspTable.txt" );

                correspRates = ( comparer.rateCorrespondence() );
                if( correspRates.size() != 6 )
                {
                    std::cerr << "ERROR: correspondance ratings vector has wrong size" << std::endl;
                }
                else
                {

                    outFile << "Size-Match Correlation:\t " << correspRates[0] << std::endl;
                    logFile << "Size-Match Correlation:\t " << correspRates[0] << std::endl;

                    outFile << "Mean Match Distance:\t " << correspRates[1] << std::endl;
                    logFile << "Mean Match Distance:\t " << correspRates[1] << std::endl;

                    outFile << "Size-Weighted Match Distance:\t " << correspRates[2] << std::endl;
                    logFile << "Size-Weighted Match Distance:\t " << correspRates[2] << std::endl;

                    outFile << "% of matches:\t " << 100.0 * correspRates[3] << std::endl;
                    logFile << "% of matches:\t " << 100.0 * correspRates[3] << std::endl;

                    outFile << "Mean Euclidean Distance:\t " << correspRates[4] << std::endl;
                    logFile << "Mean Euclidean Distance:\t " << correspRates[4] << std::endl;

                    outFile << "Size-Weighted Euclidean Distance:\t " << correspRates[5] << std::endl;
                    logFile << "Size-Weighted Euclidean Distance:\t " << correspRates[5] << std::endl;
                }

                break;

            case CRSP_MERGE:
                logFile <<"Correspondance mode:\t MergeUp"<<std::endl;
                comparer.mergedUpCorrespondence();
                break;

            default:
                std::cerr << "ERROR: unrecognized correspondance mode"<<std::endl;
                exit(-1);
                break;
            } // end switch



            if( matchNoise )
            {

                std::string outpuAlphatFileName(outputFolder+"/alphaTest.txt");
                std::ofstream outAlpha(outpuAlphatFileName.c_str());
                if(!outAlpha) {
                    std::cerr << "ERROR: unable to open out file: \""<<outpuAlphatFileName<<"\""<<std::endl;
                    exit(-1);
                }

                outAlpha << "Alpha gran1 gran2 meanGran tCPCC %pairs effGran: "<<std::endl;
                if(verbose)
                {
                    std::cout<< "Starting alpha loop "<< std::endl<< std::endl;
                }


                for( float noiseAlpha = 0; noiseAlpha <= 1.01; noiseAlpha+= 0.05 )
                {

                    outAlpha << noiseAlpha << " " << std::flush;
                    WHtree treeloop1(tree1);
                    WHtree treeloop2(tree2);
                    treeComparer comparerLoop(&treeloop1, &treeloop2, comparer);



                    if( matchNoise )
                    {

                        outFile << "matching noise alpha value: "<<noiseAlpha<<std::endl;
                        logFile << "matching noise alpha value: "<<noiseAlpha<<std::endl;


                        std::pair< float, float > matchGran( comparerLoop.applyNoiseBaseline( noiseAlpha ) );

                        outFile << "Size of maxgran part 1: " << matchGran.first << std::endl;
                        outFile << "Size of maxgran part 2: " << matchGran.second << std::endl;
                        outFile << "Average maxgran size: "<< ( ( matchGran.first + matchGran.second )  / 2.0 ) << std::endl;

                        logFile << "Size of maxgran part 1: " << matchGran.first << std::endl;
                        logFile << "Size of maxgran part 2: " << matchGran.second << std::endl;
                        logFile << "Average maxgran size: "<< ( ( matchGran.first + matchGran.second )  / 2.0 ) << std::endl;

                        outAlpha << matchGran.first << " " << matchGran.second << " " << ( matchGran.first + matchGran.second )  / 2.0  << " " << std::flush;

                    }

                    if( !nocomp )
                    {
                        std::pair< std::pair< float, float >, std::pair< float, float > > tCPCC(comparerLoop.doTcpcc());

                        outFile << "weighted tCPCT:\t " << tCPCC.first.first << std::endl;
                        logFile << "weighted tCPCC:\t " << tCPCC.first.first << std::endl;
                        outFile << "simple CPCC:\t " << tCPCC.first.second << std::endl;
                        logFile << "simple CPCC:\t " << tCPCC.first.second << std::endl;
                        outFile << "% used pairs:\t " << 100*tCPCC.second.first << std::endl;
                        logFile << "% used pairs:\t " << 100*tCPCC.second.first << std::endl;
                        outFile << "effect. Granularity:\t " << tCPCC.second.second << std::endl;
                        logFile << "effect. Granularity:\t " << tCPCC.second.second << std::endl;

                        outAlpha << tCPCC.first.first << " " << 100*tCPCC.second.first << " " << tCPCC.second.second  << " " << std::flush;


                        if( !onlyCpcc )
                        {
                            std::pair< float, float > sTriplets(comparerLoop.simpleTriplets(sampleFreq));

                            outFile << "Simple Triplets Unweighted:\t " << sTriplets.first << std::endl;
                            logFile << "Simple Triplets Unweighted:\t " << sTriplets.first << std::endl;

                            outFile << "Simple Triplets Size-Weighted:\t " << sTriplets.second << std::endl;
                            logFile << "Simple Triplets Size-Weighted:\t " << sTriplets.second << std::endl;
                        }
                    }


                    logFile << "Final Tree 1: " << treeloop1.getReport( false ) << std::endl;
                    logFile << "Final Tree 2: " << treeloop2.getReport( false ) << std::endl;
                    logFile << comparerLoop.reportBaseNodes() << std::endl;



                    if (verbose) {
                        std::cout << "Final Tree 1: " << treeloop1.getReport( false ) << std::endl;
                        std::cout << "Final Tree 2: " << treeloop2.getReport( false ) << std::endl;
                        std::cout << comparerLoop.reportBaseNodes() << std::endl;
                    }

                    if(verbose)
                    {
                        std::cout<< std::endl;
                    }

                    outAlpha << std::endl;

                } // end alpha for
            } // end matchnoise if
            else
            {
                if( !nocomp )
                {
                    std::pair< std::pair< float, float >, std::pair< float, float > > tCPCC(comparer.doTcpcc());

                    outFile << "weighted tCPCT:\t " << tCPCC.first.first << std::endl;
                    logFile << "weighted tCPCC:\t " << tCPCC.first.first << std::endl;
                    outFile << "simple CPCC:\t " << tCPCC.first.second << std::endl;
                    logFile << "simple CPCC:\t " << tCPCC.first.second << std::endl;
                    outFile << "% used pairs:\t " << 100*tCPCC.second.first << std::endl;
                    logFile << "% used pairs:\t " << 100*tCPCC.second.first << std::endl;
                    outFile << "effect. Granularity:\t " << tCPCC.second.second << std::endl;
                    logFile << "effect. Granularity:\t " << tCPCC.second.second << std::endl;


                    if( !onlyCpcc )
                    {
                        std::pair< float, float > sTriplets(comparer.simpleTriplets(sampleFreq));

                        outFile << "Simple Triplets Unweighted:\t " << sTriplets.first << std::endl;
                        logFile << "Simple Triplets Unweighted:\t " << sTriplets.first << std::endl;
                        outFile << "Simple Triplets Size-Weighted:\t " << sTriplets.second << std::endl;
                        logFile << "Simple Triplets Size-Weighted:\t " << sTriplets.second << std::endl;
                    }
                }// end matchnoise else


                logFile << "Final Tree 1: " << tree1.getReport( false ) << std::endl;
                logFile << "Final Tree 2: " << tree2.getReport( false ) << std::endl;
                logFile << comparer.reportBaseNodes() << std::endl;



                if (verbose) {
                    std::cout << "Final Tree 1: " << tree1.getReport( false ) << std::endl;
                    std::cout << "Final Tree 2: " << tree2.getReport( false ) << std::endl;
                    std::cout << comparer.reportBaseNodes() << std::endl;
                }

                if(verbose)
                {
                    std::cout<< std::endl;
                }
            }

            tree1.writeTree( outputFolder + "/treeCompared1.txt" );
            tree2.writeTree( outputFolder + "/treeCompared2.txt" );

        } // end no-rand if

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

