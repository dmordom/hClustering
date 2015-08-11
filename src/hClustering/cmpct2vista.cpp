//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
// This file is also part of OpenWalnut ( http://www.openwalnut.org ).
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
//---------------------------------------------------------------------------
//
//  cmpct2vista
//
//  Transform a 1D .cmpct compact tract vector into a 1D .v vista image vector.
//
//  WARNING: this program will port the vector as-is, and is meant simply to be able to easily visualize the contents and facilitate conversion to other formats.
//            a .cmpct vector compacted from nifti coordinates and then tranformed to .v with this program, will not prodcue a correct image if then blown to full 3D with compact2full.
//            to transform a .cmpct tract to a fully corresponding .v tract, firstly blow to 3D .nii image, then convert to vista with vnifti2image and then compact with full2compact.
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -i --input:      [mutually exclusive with -f] Input .cmpct tractogram to be converted into vista 1D vector, multiple inputs allowed separated by spaces.
//
//   -f --filenames:  [mutually exclusive with -i] Text file with a list of multiple input filenames.
//
//  [-o --output]:    Output file or folder to write vista 1D tracts.
//
//  [-z --zip]:       zip output files.
//
//  [-F --ufloat]:    use float32 representation to write output tracts (default is uint8).
//
//  [-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors.
//
//
//  example:
//
//  cmpct2vista -i tract.cmpct -o tract.v
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

// boost library
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// classes
#include "fileManagerFactory.h"
#include "fileManager.h"

#define SUFFIX ""

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

        std::string progName("cmpct2vista");
        std::string configFilename("/home/raid2/moreno/Code/hClustering/config/"+progName+".cfg");

        // program parameters
        std::string outputFilename, inputListFilename;
        std::vector< std::string > inputFilenameVector;
        unsigned int threads(0);
        bool useFloat( false ), outIsFolder( false ), doZip(true);

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ("version", "print version string")
                ("help,h", "produce help message")
                ("input,i",      boost::program_options::value< std::vector<std::string > >(&inputFilenameVector)->multitoken(), "[xor with -f] input file(s)")
                ("filenames,f",      boost::program_options::value< std::string >(&inputListFilename), "[xor with -i] text file with a list of input filenames")
                ("output,o",      boost::program_options::value< std::string >(&outputFilename), "output filename or output directory")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ("zip,z", "[opt] zip output files")
                ("ufloat,F", "[opt] use float32 representation to write tracts (default is uint8)")
                ("pthreads,p",  boost::program_options::value< unsigned int >(&threads), "[opt] number of processing threads to run the program in parallel, default: all available")
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
            std::cout << "cmpct2vista" << std::endl << std::endl;
            std::cout << "Transform a 1D .cmpct compact tract vector into a 1D .v vista image vector." << std::endl << std::endl;
            std::cout << "WARNING: this program will port the vector as-is, and is meant simply to be able to easily visualize the contents and facilitate conversion to other formats." << std::endl;
            std::cout << "          a .cmpct vector compacted from nifti coordinates and then tranformed to .v withth this program, will not prodcue a correct image if then blown to full 3D with compact2full." << std::endl;
            std::cout << "          to transform a .cmpct tract to a fully corresponding .v tract, firstly blow to 3D .nii image, then convert to vista with vnifti2image and then compact with full2compact." << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       Produce extended program help message." << std::endl << std::endl;
            std::cout << " -i --input:      [mutually exclusive with -f]  Input .cmpct tractogram to be converted into vista 1D vector, multiple inputs allowed separated by spaces." << std::endl << std::endl;
            std::cout << " -f --filenames:  [mutually exclusive with -i] Text file with a list of multiple input filenames." << std::endl << std::endl;
            std::cout << "[-o --output]:    Output file or folder to write vista 1D tracts." << std::endl << std::endl;
            std::cout << "[-z --zip]:       zip output files." << std::endl << std::endl;
            std::cout << "[-F --ufloat]:    use float32 representation to write output tracts (default is uint8)." << std::endl << std::endl;
            std::cout << "[-p --pthreads]:  Number of processing threads to run the program in parallel. Default: use all available processors." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "example:" << std::endl << std::endl;
            std::cout << "cmpct2vista -i tract.cmpct -o tract.v" << std::endl << std::endl;
            exit(0);
        }
        if (variableMap.count("version")) {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }


        if (variableMap.count("ufloat"))
        {
            std::cout << "writing in float"<<std::endl;
            useFloat=true;
        }
        else
        {
            std::cout << "writing in char"<<std::endl;
            useFloat=false;
        }

        if (variableMap.count("zip"))
        {
            std::cout << "zipping output files"<<std::endl;
            doZip=true;
        }
        else
        {
            doZip=false;
        }


        if (variableMap.count("pthreads")) {
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

        if ( !variableMap.count("input") && !variableMap.count("filenames") )
        {
            std::cerr << "ERROR: no input tract file or filenames stated, please use either option -i or -f"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if ( variableMap.count("input") && variableMap.count("filenames") )
        {
            std::cerr << "ERROR: please use only either option -i or -f"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }

        if (variableMap.count("input"))
        {
            std::cout << "Tractogram files: "<<  std::flush;

            for(size_t i=0; i<inputFilenameVector.size(); ++i )
            {
                std::string inputFilename(inputFilenameVector[i]);
                if(!boost::filesystem::is_regular_file(boost::filesystem::path(inputFilename)))
                {
                    std::cerr << std::endl << "ERROR: tract file \""<<inputFilename<<"\" is not a regular file"<<std::endl;
                    std::cerr << visibleOptions << std::endl;
                    exit(-1);
                }
                std::cout << inputFilename << "_" << std::flush;
            }
            std::cout<<std::endl;
        }

        if (variableMap.count("filenames"))
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path(inputListFilename)))
            {
                std::cerr << std::endl << "ERROR: tract filenames file \""<<inputListFilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tractogram filenames file: "<< inputListFilename << std::endl;
        }


        if (variableMap.count("output"))
        {
            if(boost::filesystem::is_directory(boost::filesystem::path(outputFilename)))
            {
                std::cout << "output folder: "<< outputFilename << std::endl;
                outIsFolder = true;
            }
            else
            {
                std::cout << "output file: "<< outputFilename << std::endl;
                outIsFolder = false;
            }
        }
        else
        {
            std::cerr << "ERROR: missing output file/folder"<<std::endl;
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        // ==============================================================================


        std::vector<std::string> finalInputs;
        if( !inputFilenameVector.empty() )
        {
            finalInputs = inputFilenameVector;
        }
        else if ( !inputListFilename.empty() )
        {
            std::ifstream ifs( inputListFilename.c_str(), std::ifstream::in );
            std::string line;

            if( !ifs.is_open() )
            {
                std::cerr << "ERROR: could not open filenames list file"<<std::endl;
                ifs.close();
                exit(-1);
            }

            while( !ifs.eof() )
            {
                getline( ifs, line );
                finalInputs.push_back( line );
            }

            ifs.close();
            std::cout << finalInputs.size() << " input filenames read from file" <<std::endl;
        }
        else
        {
            std::cerr << "ERROR: no input files??"<<std::endl;
            exit(-1);
        }


        if( finalInputs.size() > 1 && !outIsFolder )
        {
            std::cerr << "ERROR: multiple input files but output is a filename, for multiple inputs please indicate an output directory"<<std::endl;
            exit(-1);
        }

        //////////////////////////////

        // file io classes
        fileManagerFactory ioFMFactory;
        ioFMFactory.isNifti();
        fileManager& inputFM(ioFMFactory.getFM());
        ioFMFactory.isVista();
        fileManager& outputFM(ioFMFactory.getFM());


        ValueType repValueType;
        inputFM.readAsUnThres();
        inputFM.readAsLog();

        if( useFloat )
        {
            outputFM.writeInFloat();
            repValueType = VTFloat32;
        }
        else
        {
            outputFM.writeInChar();
            repValueType = VTUINT8;
        }
        if( doZip )
        {
            outputFM.storeZipped();
        }
        else
        {
            outputFM.storeUnzipped();
        }

        #pragma omp parallel for schedule( static )
        for( size_t index = 0; index < finalInputs.size(); ++index )
        {
            std::string thisInput(finalInputs[index]);

            std::string extension( boost::filesystem::path(thisInput).extension().string() );
            std::string stem( boost::filesystem::path(thisInput).stem().string() );
            if ( extension == ".gz" )
            {
                extension = boost::filesystem::path(stem).extension().string();
                stem = boost::filesystem::path(stem).stem().string();
            }


            if( extension != COMPACT_EXT )
            {
                std::cerr << "ERROR: input files dont have .cmpct extension (.cmpct.gz IS also allowed) "<<std::endl;
                exit(-1);
            }



            compactTract compactVect;
            inputFM.readTract( thisInput,  &compactVect);

            if(!outputFilename.empty())
            {
                std::string outPath;
                if(outIsFolder)
                {
                    outPath = outputFilename + "/" + stem + SUFFIX + VISTA_EXT;
                }
                else
                {
                    outPath = outputFilename;
                }
                std::cout<< "writing file: "<< outPath <<std::endl;
                outputFM.writeTract( outPath, compactVect );
            }
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



