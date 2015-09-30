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
//  surfprojection
//
//  Performs an interpolated projection from the roi seed voxels to the vertices of a (freesurfer) surface.
//   Options include nearest neighbor, averaging, and gaussian interpolation.
//
//  * Arguments:
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//  [-k --kradius]:  Kernel radius (in voxel dimension units). Use 0 for nearest neighbor interpolation (default) and > 0 for average interpolation.
//
//  [-g --gauss]:    [use only with -k and radius > 0] Use gaussian smoothing instead of average. indicate full-width half-maximum (in voxel dimension units).
//
//   -r --roifile:    File with the seed voxels coordinates.
//
//   -s --surffile:   File with the surface vertex coordinates.
//
//   -I --inputf:     Input tractogram folder (leaf tractograms).
//
//   -O --outputf:    Output folder where resulting tracts will be written.
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
//   surfprojection -k 6 -g 3 -r roi.txt -s surf.txt -I leaftracts/ -O output/ -v
//
//  * Outputs (in output folder defined at option -O):
//
//   - 'compact_X.cmpct(.v)' (where X is the corresponding surface vertex ID): A compact tractogram with the mean tractogram projected to vertex X.
//   - 'surfprojection_log.txt' - A text log file containing the parameter details and in-run and completion information of the program.
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

// parallel execution
#include <omp.h>

// boost library
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// classes
#include "surfProjecter.h"




int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        boost::filesystem::path workingDir( boost::filesystem::current_path());


        // ========== PROGRAM PARAMETERS ==========

        std::string progName("surfprojection");
        std::string configFilename("../../config/"+progName+".cfg");

        // program parameters
        std::string roiFilename, surfFilename, inputFolder, outputFolder;
        float kernelRadius(0),fwhm(0);
        unsigned int threads(0);
        bool gauss(false), verbose(false), useFloat( false ), doZip( false ), niftiMode( true );

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "kradius,k",  boost::program_options::value< float >(&kernelRadius), "[opt] Kernel radius (in voxel dimension units). Use 0 for nearest neighbor interpolation (default) and > 0 for average interpolation.")
                ( "gauss,g", boost::program_options::value< float >(&fwhm),  "[opt | use only with -k and radius > 0] Use gaussian smoothing instead of average. indicate fwhm (in voxel dimension units).")
                ( "roifile,r",  boost::program_options::value< std::string >(&roiFilename), "File with the seed voxels coordinates")
                ( "surffile,s",  boost::program_options::value< std::string >(&surfFilename), "File with the surface vertex coordinates")
                ( "inputf,I",  boost::program_options::value< std::string >(&inputFolder), "Input tractogram folder (leaf tractograms)")
                ( "outputf,O",  boost::program_options::value< std::string >(&outputFolder), "Output folder where resulting tracts will be written")
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




        if (variableMap.count("help"))
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
            std::cout << "surfprojection" << std::endl << std::endl;
            std::cout << "Performs an interpolated projection from the roi seed voxels to the vertices of a (freesurfer) surface." << std::endl;
            std::cout << " Options include nearest neighbor, averaging, and gaussian interpolation." << std::endl << std::endl;
            std::cout << "* Arguments:" << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << "[-k --kradius]:  Kernel radius (in voxel dimension units). Use 0 for nearest neighbor interpolation (default) and > 0 for average interpolation." << std::endl << std::endl;
            std::cout << "[-g --gauss]:    [use only with -k and radius > 0] Use gaussian smoothing instead of average. indicate full-width half-maximum (in voxel dimension units)." << std::endl << std::endl;
            std::cout << " -r --roifile:    File with the seed voxels coordinates." << std::endl << std::endl;
            std::cout << " -s --surffile:   File with the surface vertex coordinates." << std::endl << std::endl;
            std::cout << " -I --inputf:     Input tractogram folder (leaf tractograms)." << std::endl << std::endl;
            std::cout << " -O --outputf:    Output folder where resulting tracts will be written." << std::endl << std::endl;
            std::cout << "[-v --verbose]:   Verbose output (recommended)." << std::endl << std::endl;
            std::cout << "[--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files." << std::endl << std::endl;
            std::cout << "[-z --zip]:       Zip output files." << std::endl << std::endl;
            std::cout << "[-F --ufloat]:    Use float32 representation to write output tracts (default is uint8)." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Usage example:" << std::endl << std::endl;
            std::cout << " surfprojection -k 6 -g 3 -r roi.txt -s surf.txt -I leaftracts/ -O output/ -v" << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "* Outputs (in output folder defined at option -O):" << std::endl << std::endl;
            std::cout << " - 'compact_X.cmpct(.v)' (where X is the corresponding surface vertex ID): A compact tractogram with the mean tractogram projected to vertex X." << std::endl;
            std::cout << " - 'surfprojection_log.txt' - A text log file containing the parameter details and in-run and completion information of the program." << std::endl;
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

        if (variableMap.count("roifile"))
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(roiFilename)))
            {
                std::cerr << "ERROR: roi file \""<<roiFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Roi voxels file: "<< roiFilename << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no roi file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count("surffile"))
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(surfFilename)))
            {
                std::cerr << "ERROR: surf file \""<<surfFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "surface vertex coords file: "<< surfFilename << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no surf file stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("kradius"))
        {
            if(kernelRadius < 0)
            {
                std::cerr << "ERROR: kernel size must be positive value"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            else if (kernelRadius == 0)
            {
                if( verbose )
                {
                    std::cout << "Nearest neighbor interpolation" << std::endl;
                }
            }
            else
            {
                if( verbose )
                {
                    std::cout << "Kernel Size: "<< kernelRadius << std::endl;
                }
                if (variableMap.count("gauss"))
                {
                    if( verbose )
                    {
                        std::cout << "Gaussian smoothing kernel" << std::endl;
                    }
                    gauss = true;
                }
                else
                {
                    if( verbose )
                    {
                        std::cout << "Average smoothing kernel" << std::endl;
                    }
                    gauss = false;
                }
            }
        }

        if ( variableMap.count( "gauss" ) && kernelRadius == 0 )
        {
            std::cerr << "WARNING: Nearest neighbor interpolation chosen, -g option will be ignored"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if ( variableMap.count( "inputf" ) )
        {
            if(!boost::filesystem::is_directory(boost::filesystem::path(inputFolder))) {
                std::cerr << "ERROR: input tractogram folder \""<<inputFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);

            }
            if( verbose )
            {
                std::cout << "Input tractogram folder: "<< inputFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no input tract folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

        }

        if ( variableMap.count( "outputf" ) )
        {
            if( !boost::filesystem::is_directory(boost::filesystem::path( outputFolder ) ) )
            {
                std::cerr << "ERROR: output folder \""<<outputFolder<<"\" is not a directory"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            if( verbose )
            {
                std::cout << "Output folder: "<< outputFolder << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: no output folder stated"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);

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
        if(!logFile)
        {
            std::cerr << "ERROR: unable to open log file: \""<<logFilename<<"\""<<std::endl;
            exit(-1);
        }


        logFile <<"Start Time:\t"<< ctime(&programStartTime) <<std::endl;
        logFile <<"Working directory:\t"<< workingDir.string() <<std::endl;
        logFile <<"Roi file:\t"<< roiFilename <<std::endl;
        logFile <<"kernel size:\t"<< kernelRadius <<std::endl;
        if (gauss)
        {
            logFile <<"Gaussian kernel" <<std::endl;
        }
        else
        {
            logFile <<"Square kernel (mean)" <<std::endl;
        }
        logFile <<"Surf file:\t"<< surfFilename <<std::endl;
        logFile <<"Input folder:\t"<< inputFolder <<std::endl;
        logFile <<"Output folder:\t"<< outputFolder <<std::endl;
        logFile <<"Verbose:\t"<< verbose <<std::endl;
        logFile <<"Processors used:\t"<< threads <<std::endl;
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

        /////////////////////////////////////////////////////////////////

        surfProjecter sdmatch( roiFilename, surfFilename, verbose );

        if( kernelRadius > 0 )
        {
            if ( gauss )
            {
                sdmatch.kernelGauss( kernelRadius, fwhm );
            }
            else
            {
                sdmatch.kernelMean( kernelRadius );
            }
            sdmatch.matchCoordsKernel();
        }
        else
        {
            sdmatch.matchCoordsNearestNb();
        }

        sdmatch.writeMeanTracts( inputFolder, outputFolder, useFloat, doZip );


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

