//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
//
// For more reference on the underlying algorithm and research they have been used for refer to:
// - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014).
//   A hierarchical method for whole-brain connectivity-based parcellation.
//   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528
// - Moreno-Dominguez, D. (2014).
//   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography.
//   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig.
//   ISBN 978-3-941504-45-5
//
// hClustering is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// http://creativecommons.org/licenses/by-nc/3.0
//
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------
//
//  randtracts
//
//  Generates a set of tractograms matching in number and name to those of an existing roi file but that will generate a randomly uniform dissimilarity matrix when compouting their distance.
//   This is intended to be able to establish a baseline for tree quality values, by testing the quality measures on trees built over tractograms with no similarity structure (random).
//   For this purpose the tractograms must be vectors uniformly distributed on the positive quadrant surface of the unit hypersphere.
//   The main method achieves this by first generating normally distributed vectors in the hyperspace, then normalizing these vectors to the unit hyperphere.
//    However, for this approximation to be accurate dimension must be relatively low, therefore a dimension=10 is recommended (default).
//   The alternative method generates uniformly distributed vectors in the hyperspace, then filters out elements outside the unit hypersphere, and normalizes the remaining elements to the surface.
//    This method has higher accuracy than the main method at higher dimensions, but computing time also increases exponentially, as most elements must be filetered out.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -r --roi:        Roi file with leaf coordinates/trackIDs of tractograms to generate.
//
//   -O --outputf:    Output folder where results be written.
//
//  [-d --dim]:       Desired dimension of the random tractograms. Default: 10.
//
//  [-s --seeds]:     Random number generator seed. Change in order to obtain a different set of results. Same seed will always be reproducible. Default: 0.
//
//  [--alt]:          Use alternative method: filtered uniform hyperspace vector. Better approximation at higher dimension values but much more time-consuming.
//
//  [-v --verbose]:   verbose output (recommended).
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files.
//
//  [-z --zip]:       zip output files.
//
//  [-F --ufloat]:    use float32 representation to write output tracts (default is uint8).
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  * Usage example:
//
//   randtracts -r roi.txt -O results/ -d 10 -s 0 -v
//
//
//  * Outputs (in output folder defined at option -O):
//
//   - 'probtract_X.cmpct' (default)- (where X is a tract ID) artificial compact tractograms that would yield a uniformly random distance matrix.
//.  - 'connect_X_Y_Z.v' (--vista option)- (where XYZ are tract seed voxel coordinates) artificial compact tractograms in vista format that would yield a uniformly random distance matrix.
//   - 'randtracts_log.txt' - A text log file containing the parameter details and in-run and completion information of the program.
//
//---------------------------------------------------------------------------

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
#include "roiLoader.h"
#include "WHcoord.h"
#include "fileManagerFactory.h"
#include "compactTract.h"





int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("randtracts");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, outputFolder;
        unsigned int threads( 0 );
        size_t dimension( 10 ), randSeed( 0 );
        bool alternative( false );
        bool verbose( false ), niftiMode( true ), doZip( false ), useFloat( false );


        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "roi,r",  boost::program_options::value< std::string >( &roiFilename ), "File with the seed voxel coordinates" )
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "output folder" )
                ( "dim,d",  boost::program_options::value< size_t >(&dimension)->implicit_value(10), "[opt] Desired dimension of the random tractograms. Default: 10")
                ( "seed,s",  boost::program_options::value< size_t >(&randSeed)->implicit_value(0), "[opt] Random number generator seed. Default: 0")
                ( "alt", "[opt] Use alternative method: filtered uniform hyperspace vector. Better approximation at higher dimension values but much more time-consuming.")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "verbose,v", "[opt] verbose output." )
                ( "vista", "[opt] use vista file format (default is nifti)." )
                ( "zip,z", "[opt] zip output files.")
                ( "ufloat,F", "[opt] use float32 representation to write tracts (default is uint8)")
                ( "pthreads,p",  boost::program_options::value< unsigned int >(&threads), "[opt] number of processing cores to run the program in. Default: all available." )
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



        if ( variableMap.count( "help" ) )
        {
            std::cout << "---------------------------------------------------------------------------" << std::endl;
            std::cout << std::endl;
            std::cout << " Project: hClustering" << std::endl;
            std::cout << std::endl;
            std::cout << " Whole-Brain Connectivity-Based Hierarchical Parcellation Project" << std::endl;
            std::cout << " David Moreno-Dominguez" << std::endl;
            std::cout << " d.mor.dom@gmail.com" << std::endl;
            std::cout << " moreno@cbs.mpg.de" << std::endl;
            std::cout << " www.cbs.mpg.de/~moreno" << std::endl;
            std::cout << std::endl;
            std::cout << " For more reference on the underlying algorithm and research they have been used for refer to:" << std::endl;
            std::cout << " - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014)." << std::endl;
            std::cout << "   A hierarchical method for whole-brain connectivity-based parcellation." << std::endl;
            std::cout << "   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528" << std::endl;
            std::cout << " - Moreno-Dominguez, D. (2014)." << std::endl;
            std::cout << "   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography." << std::endl;
            std::cout << "   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig." << std::endl;
            std::cout << "   ISBN 978-3-941504-45-5" << std::endl;
            std::cout << std::endl;
            std::cout << " hClustering is free software: you can redistribute it and/or modify" << std::endl;
            std::cout << " it under the terms of the GNU Lesser General Public License as published by" << std::endl;
            std::cout << " the Free Software Foundation, either version 3 of the License, or" << std::endl;
            std::cout << " (at your option) any later version." << std::endl;
            std::cout << " http://creativecommons.org/licenses/by-nc/3.0" << std::endl;
            std::cout << std::endl;
            std::cout << " hClustering is distributed in the hope that it will be useful," << std::endl;
            std::cout << " but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
            std::cout << " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
            std::cout << " GNU Lesser General Public License for more details." << std::endl;
            std::cout << std::endl;
            std::cout << "---------------------------------------------------------------------------" << std::endl << std::endl;
            std::cout << "randtracts" << std::endl << std::endl;
            std::cout << "Generates a set of tractograms matching in number and name to those of an existing roi file but that will generate a randomly uniform dissimilarity matrix when compouting their distance." << std::endl;
            std::cout << " This is intended to be able to establish a baseline for tree quality values, by testing the quality measures on trees built over tractograms with no similarity structure (random)." << std::endl;
            std::cout << " For this purpose the tractograms must be vectors uniformly distributed on the positive quadrant surface of the unit hypersphere." << std::endl;
            std::cout << " The main method achieves this by first generating normally distributed vectors in the hyperspace, then normalizing these vectors to the unit hyperphere." << std::endl;
            std::cout << "  However, for this approximation to be accurate dimension must be relatively low, therefore a dimension=10 is recommended (default)." << std::endl;
            std::cout << " The alternative method generates uniformly distributed vectors in the hyperspace, then filters out elements outside the unit hypersphere, and normalizes the remaining elements to the surface." << std::endl;
            std::cout << "  This method has higher accuracy than the main method at higher dimensions, but computing time also increases exponentially, as most elements must be filetered out." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -r --roi:        Roi file with leaf coordinates/trackIDs of tractograms to generate." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where results be written." << std::endl << std::endl;
            std::cout << "[-d --dim]:       Desired dimension of the random tractograms. Default: 10." << std::endl << std::endl;
            std::cout << "[-s --seeds]:     Random number generator seed. Change in order to obtain a different set of results. Same seed will always be reproducible. Default: 0." << std::endl << std::endl;
            std::cout << "[--alt]:          Use alternative method: filtered uniform hyperspace vector. Better approximation at higher dimension values but much more time-consuming." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   Verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files." << std::endl << std::endl;
            std::cout << "[-z --zip]:       Zip output files." << std::endl << std::endl;
            std::cout << "[-F --ufloat]:    Use float32 representation to write output tracts (default is uint8)." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " randtracts -r roi.txt -O results/ -d 10 -s 0 -v" << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " - 'probtract_X.cmpct' (default)- (where X is a tract ID) artificial compact tractograms that would yield a uniformly random distance matrix." << std::endl;
            std::cout << " - 'connect_X_Y_Z.v' (--vista option)- (where XYZ are tract seed voxel coordinates) artificial compact tractograms in vista format that would yield a uniformly random distance matrix." << std::endl;
            std::cout << " - 'randtracts_log.txt' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
            std::cout << std::endl;
            exit(0);
        }


        if ( variableMap.count( "version" ) )
        {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }
        if (variableMap.count( "verbose" ) )
        {
            std::cout << "verbose output"<<std::endl;
            verbose=true;
        }

        if ( variableMap.count( "pthreads" ) )
        {
            if ( threads == 1 )
            {
                std::cout << "Using a single processor" << std::endl;
            }
            else if( threads == 0 || threads >= omp_get_num_procs() )
            {
                threads = omp_get_num_procs();
                std::cout << "Using all available processors ( " << threads << " )." << std::endl;
            }
            else
            {
                std::cout << "Using a maximum of " << threads << " processors " << std::endl;
            }
            omp_set_num_threads( threads );
        }
        else
        {
            threads = omp_get_num_procs();
            omp_set_num_threads( threads );
            std::cout << "Using all available processors ( " << threads << " )." << std::endl;
        }

        if ( variableMap.count( "vista" ) )
        {
            if( verbose )
            {
                std::cout << "Using vista format" << std::endl;
            }
            fileManagerFactory fmf;
            fmf.setVista();
            niftiMode = false;
        }
        else
        {
            if( verbose )
            {
                std::cout << "Using nifti format" << std::endl;
            }
            fileManagerFactory fmf;
            fmf.setNifti();
            niftiMode = true;
        }


        if ( variableMap.count( "roi" ) )
        {
            if( !boost::filesystem::is_regular_file(boost::filesystem::path( roiFilename ) ) )
            {
                std::cerr << "ERROR: roi file \"" << roiFilename << "\" is not a regular file" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            std::cout << "Roi voxels file: "<< roiFilename << std::endl;
        }
        else
        {
            std::cerr << "ERROR: no roi file stated" <<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if ( variableMap.count( "outputf" ) )
        {
            if( !boost::filesystem::is_directory( boost::filesystem::path( outputFolder ) ) )
            {
                std::cerr << "ERROR: output folder \"" <<outputFolder<< "\" is not a directory" << std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            else if( verbose )
            {
                std::cout << "Output folder: " << outputFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no output folder stated" << std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }



        if( verbose )
        {
            std::cout << "Tractograms dimension: "<< dimension << std::endl;
            std::cout << "Seed: "<< randSeed << std::endl;
        }


        if ( variableMap.count( "alt" ) )
        {
            if( verbose )
            {
                std::cout << "Using alternative method: filtered uniform hyperspace vector." << std::endl;
            }
            alternative=true;
        }
        else
        {
            if( verbose )
            {
                std::cout << "Using main method: normalized uniform hypersphere vector." << std::endl;
            }
            alternative=false;
        }


        if (variableMap.count("ufloat"))
        {
            if( verbose )
            {
                std::cout << "writing in float"<<std::endl;
            }
            useFloat=true;
        }
        else
        {
            if( verbose )
            {
                std::cout << "writing in char"<<std::endl;
            }
            useFloat=false;
        }

        if (variableMap.count("zip"))
        {
            if( verbose )
            {
                std::cout << "zipping output files"<<std::endl;
            }
            doZip=true;
        }
        else
        {
            doZip=false;
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
        logFile <<"Roi file:\t"<< roiFilename <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Tractogram dimension:\t"<< dimension << std::endl;
        logFile <<"Seed:\t"<< randSeed << std::endl;
        if(alternative)
        {
            logFile <<"Using alternative method: filtered uniform hyperspace vector"<<std::endl;
        }
        else
        {
            logFile <<"Using main method: normalized uniform hypersphere vector"<<std::endl;
        }
        if( niftiMode )
        {
            logFile << "Using nifti file format" << std::endl;
        }
        else
        {
            logFile << "Using vista file format" << std::endl;
        }
        if( useFloat )
        {
            logFile << "Writing in float32" << std::endl;
        }
        else
        {
            logFile << "Writing in uint8" << std::endl;
        }
        if( doZip )
        {
            logFile << "Zipping output tracts" << std::endl;
        }
        logFile <<"-------------"<<std::endl;


        //====================== Read Roi ============================

        std::vector<WHcoord> roivect;
        std::vector<size_t> trackids;
        HC_GRID datasetGrid;
        WHcoord datasetSize;
        size_t numStreamlines( 0 );

        RoiLoader roiLoader( niftiMode );
        roiLoader.readRoi( roiFilename, &datasetGrid, &datasetSize, &numStreamlines, &roivect, &trackids );

        if(verbose)
        {
            std::cout << "Roi loaded, " << roivect.size() << " seed voxels" << std::endl;
        }


        //==================== Generate Tracts ======================

        time_t loopStart(time(NULL)), lastTime(time(NULL));
        volatile size_t threadCount(0);
        fileManagerFactory fileMF( outputFolder );
        fileManager& fileMngr( fileMF.getFM() );
        if( useFloat )
        {
            fileMngr.writeInFloat();
        }
        else
        {
            fileMngr.writeInChar();
        }

        if( doZip )
        {
            fileMngr.storeZipped();
        }
        else
        {
            fileMngr.storeUnzipped();
        }

        // This is the underlying integer random number generator
        boost::mt19937 igen;
        igen.seed(randSeed);

        // normally distributed random number generator mean 0 variance 1
        boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  randnPrng(igen, boost::normal_distribution<>(0,1));

        // uniformly distributed random number generator between 0 and 1
        boost::variate_generator<boost::mt19937, boost::uniform_01<> >  u01Prng(igen, boost::uniform_01<>());


        if (verbose)
        {
            std::cout <<"Generating tractograms uniformly distributed over the hypersphere surface, Tract size: "<<dimension<<"..."<< std::endl;
        }

        for (size_t i=0; i<roivect.size(); ++i)
        {

            std::vector<float> randTract;
            randTract.reserve(dimension);

            if (!alternative)
            {

                float radius(0);
                // first create a vector of normally distributed values
                for (size_t j=0; j<dimension; ++j )
                {
                    float value(fabs(randnPrng()));
                    randTract.push_back((value));
                    radius+=value*value;
                }
                //calculate the lenght of the vector
                radius=sqrt(radius);
                // divide the vector by its radius and make absolute value to normalize it to the first quadrant of the unit hypersphere
                for (int j=0; j<randTract.size(); ++j )
                {
                    randTract[j]=(randTract[j]/radius);
                }
                // now the vector is part of a uniformly distributed set of vector over the positive quadrant surface of the unit hypersphere

            }
            else
            {

                bool notFound(true);
                while (notFound)
                {

                    randTract.clear();
                    float norm(0);

                    //first create a vector uniformly distributed over the hyperspace
                    for (size_t j=0; j<dimension; ++j )
                    {
                        float value(u01Prng());
                        randTract.push_back((value));
                        norm+=value*value;
                    }
                    norm=sqrt(norm);
                    // if the vector lieas outside the unit hyperphere, discard it, if its inside, project it to the sphere
                    if (norm<=1)
                    {
                        for (size_t j=0; j<dimension; ++j )
                        {
                            randTract[j]=(randTract[j]/norm);
                        }
                        notFound=false;
                    }
                }
            }

            compactTract tract(randTract);
            fileMngr.writeLeafTract(i,trackids,roivect,tract);

            if(verbose)
            {
                time_t currentTime(time(NULL));
                if(currentTime-lastTime>1)
                {
                    lastTime=currentTime;
                    float progress ((i*100.)/(roivect.size()));
                    std::string coutString(boost::lexical_cast<std::string>((int)progress)+" % of tractograms written ("+boost::lexical_cast<std::string>(i)+"). Expected remaining time: ");
                    if (progress>0)
                    {
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

