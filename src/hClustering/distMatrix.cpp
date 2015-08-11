// std librabry
#include <vector>
#include <set>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <ctime>
#include <cstring>
#include <cstdlib>


#include "distMatComputer.h"

#include "fileManagerFactory.h"

float threshold (0.4);

void PrintUsage( const char* name ) {
    std::cerr << std::endl;
    std::cerr << "Usage: " << name << " [options]" << std::endl;
    std::cerr << "[-path] : working path. If specified all other file/directory names will be relative to that path"<< std::endl;
    std::cerr << " -roi   : roi file with seed coordinates" << std::endl;
    std::cerr << " -tracd : tractogram folder" << std::endl;
    std::cerr << " -out   : output folder" << std::endl;
    std::cerr << " -mem   : maximum memory available (in Mb, or Gb if <=10)" << std::endl;
    std::cerr << " -rand : use random distance information" << std::endl;
    std::cerr << " -nothres : dont threshold tracts" << std::endl;
    std::cerr << "[-bsize]: block size (in thousands, default: 5 -> 5000x5000" << std::endl;
    std::cerr << "[-n]   : number of seeds on which to compute each correlation run, if absent maximum for selected memory will be used" << std::endl;
    std::cerr << "[-nth] : number of maximum threads on wich to run the program. Default(0): all available." << std::endl;
    std::cerr << "[-v]   : be verbose. Default: off" << std::endl;
    std::cerr << "[-vv]  : be very verbose. Default: off" << std::endl;


    std::cerr << std::endl;


    exit( -1 );
}

void missingVal (const char* name, const char* opt) {
    std::cerr <<"Error: Value for option \""<< opt <<"\" is missing!"<< std::endl;
    PrintUsage( name );
}

int main( int argc, char *argv[] )
{

    time_t program_start(time(NULL));

    // store command line parameters
    std::string file_path, roiFilename, inputFolder, outputFolder;

    size_t n_samples(0);
    float mem(0);
    size_t threads(0), bsize(5);
    bool rand_mode(false);
    bool verbose(false), veryvb(false);

    int cc = 1;
    while ( cc < argc )
    {
        if ( 0 == strcmp( "-path", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-path");
            file_path = argv[cc++];
        } else if ( 0 == strcmp( "-roi", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-roi");
            roiFilename = argv[cc++];
        } else if ( 0 == strcmp( "-tracd", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-tracd");
            inputFolder = argv[cc++];
        } else if ( 0 == strcmp( "-out", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-out");
            outputFolder = argv[cc++];
        } else if ( 0 == strcmp( "-mem", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-mem");
            if ( (EOF == sscanf( argv[cc], "%f", &mem)) || (0 == sscanf( argv[cc], "%f", &mem)) ) {
                std::cerr << "Error: Value of parameter " << argv[cc] <<" from option -mem cannot be interpreted!" << std::endl;
                PrintUsage( argv [0] );
            }
            ++cc;
        } else if ( 0 == strcmp( "-n", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-n");
            if ( (EOF == sscanf( argv[cc], "%d", &n_samples)) || (0 == sscanf( argv[cc], "%d", &n_samples)) ) {
                std::cerr << "Error: Value of parameter " << argv[cc] <<" from option -clevel cannot be interpreted!" << std::endl;
                PrintUsage( argv [0] );
            }
            ++cc;
        } else if ( 0 == strcmp( "-bsize", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-bsize");
            if ( (EOF == sscanf( argv[cc], "%d", &bsize )) || (0 == sscanf( argv[cc], "%d", &bsize )) ) {
                std::cerr << "Error: Value of parameter " << argv[cc] << " cannot be interpreted!" << std::endl;
                PrintUsage( argv [0] );
            }
            ++cc;
        } else if ( 0 == strcmp( "-nth", argv[cc] ) ) {
            if ( ++cc >= argc )
                missingVal(argv[0],"-nth");
            if ( (EOF == sscanf( argv[cc], "%d", &threads )) || (0 == sscanf( argv[cc], "%d", &threads )) ) {
                std::cerr << "Error: Value of parameter " << argv[cc] << " cannot be interpreted!" << std::endl;
                PrintUsage( argv [0] );
            }
            ++cc;
        } else if ( 0 == strcmp( "-rand", argv[cc] ) ) {
            ++cc;
            rand_mode = true;
        } else if ( 0 == strcmp( "-v", argv[cc] ) ) {
            ++cc;
            verbose = true;
        } else if ( 0 == strcmp( "-vv", argv[cc] ) ) {
            ++cc;
            veryvb = true;
        } else if ( 0 == strcmp( "-nothres", argv[cc] ) ) {
            ++cc;
            threshold = 0;
        } else {
            std::cerr << "Error: Unknown parameter " << argv[cc] << std::endl;
            PrintUsage( argv[0] );
        }
    }

    bool ok( true );
    if ( roiFilename.empty() ) {
        std::cerr << "Missing parameter -btd" << std::endl;
        ok = false;
    }
    if ( inputFolder.empty() ) {
        std::cerr << "Mtreeissing parameter -tracd" << std::endl;
        ok = false;
    }
    if (outputFolder.empty() ) {
        std::cerr << "Missing parameter -out" << std::endl;
        ok = false;
    }
    if ( (mem <= 0)  ) {
        std::cerr << "Error using option -mem: value must be positive" << std::endl;
        ok = false;
    }
    if ( (bsize < 0)  ) {
        std::cerr << "Error using option -bsize: value must be positive" << std::endl;
        ok = false;
    }
    if ( (n_samples < 0)  ) {
        std::cerr << "Error using option -n: value must be positive" << std::endl;
        ok = false;
    }
    if ( !ok )
        PrintUsage( argv[0] );

    if (!file_path.empty()) {
        roiFilename = file_path + "/" + roiFilename;
        inputFolder = file_path + "/" + inputFolder;
    }
    if (veryvb == true)
        verbose = true;





    // ========== Get Binary Tree Data ==========

    std::cout <<"Threshold: "<< threshold << std::endl;



    if (rand_mode) {
        std::cout << "Random option selected. A random distance matrix with the same dimensions of the seed mask introduced will be created" << std::endl;
        srand (1);

    }

    // ========== set number of threads ==========

    if (threads==1) {
        std::cout <<"Using a single processor"<< std::endl;
    } else if (threads && threads < omp_get_max_threads()) {
        std::cout <<"Using a maximum of "<< threads <<" processors "<< std::endl;
    } else {
        threads = omp_get_max_threads();
        std::cout <<"Using all available processors ("<< threads <<")." << std::endl;
    }

    // ========== Decide size of sub-block size ==========



    // if memory value introduced  is 10 or less its interpreted as GBytes, otherwise, its interpreted as MBytes
    std::cout <<"Maximum memory to be used: "<<std::flush;
    std::cout<< mem <<" GBytes "<< std::endl;

    fileManagerFactory fileMF;
    fileMF.setNifti();

    distMatComputer distMatComp( roiFilename, 0, verbose, false );
    distMatComp.setInputFolder( inputFolder );
    distMatComp.setOutputFolder( outputFolder );
    distMatComp.setBlockSize( mem, bsize );

/*
    // ========== write down lookup list (for other programs to locate easily the distances between two seeds in the filess to be written) ==========


    //initialize containers
    std::map<coordinate,size_t> roimap; //map to contain the seed voxel coordinates for fast access by coordinates
    std::vector<coordinate> roivect;    //vector to contain the seed voxel coordinates for fast access by identifier
    std::vector<std::pair<size_t,size_t> >roi_block_index; // vector to store to which block index belongs each seed


    // insert block index in lookup list
    for (int row=0 ; row < num_b ; ++row) {

        // get subsets of seeds
        size_t first_row_seed(row*bsize), postlast_row_seed((row+1)*bsize);
        if (postlast_row_seed>roivect.size())
            postlast_row_seed = roivect.size();

        std::vector<coordinate> seed_set_row(&roivect[first_row_seed], &roivect[postlast_row_seed]);
        for (std::vector<coordinate>::iterator iter(seed_set_row.begin()); iter!=seed_set_row.end();++iter)
            roi_block_index.push_back(std::make_pair(row,iter-seed_set_row.begin()));
    }

    write_output(out_dir, roivect, roi_block_index);

    // ========== Calculate distances and write down blocks ==========

    do_dist_blocks(tract_dir, out_dir, roivect, bsize, tract_block_size, tract_length, threads, rand_mode, verbose, veryvb);
*/

    time_t program_end(time(NULL));
    int total_time( difftime(program_end,program_start) );
    std::cout << "Program Finished, total time: " << total_time/3600 <<"h "<<  (total_time%3600)/60 <<"' "<< ((total_time%3600)%60) <<"\""<< std::endl;


    return 0;
}


