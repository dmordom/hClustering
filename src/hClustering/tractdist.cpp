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
//  tractdist
//
//  Compute the distance (dissimilarity) value between two tractograms read from file.
//
//   --version:       Program version.
//
//   -h --help:       Produce extended program help message.
//
//   -a --tracta:     Filename of first tractogram.
//
//   -b --tractb:     Filename of second tractogram.
//
//  [-t --threshold]: Threshold to apply directly to the tractogram values before computing the dissimilarity (in order to avoid tractography noise affect the result).
//                     Unlike in other hClusterign commands, this threshold value is an absolute value to apply to the tractogram data as is, not a relative threshold.
//                     Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired.
//
//  [--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files.
//
//
//  example:
//
//  tractdist -a tractA.nii -b tractB.nii -t 0.2
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
#include "compactTract.h"
#include "compactTractChar.h"
#include "fileManagerFactory.h"
#include "WStringUtils.h"


int main( int argc, char *argv[] )
{
//    try {

        time_t programStartTime(time(NULL));
        //boost::filesystem::path workingDir( boost::filesystem::current_path());

        // ========== PROGRAM PARAMETERS ==========

        std::string progName("tractdist");
        std::string configFilename("../../config/"+progName+".cfg");

        std::string tractAfilename, tractBfilename;
        float thresValue(0);
        // program parameters

        // Declare a group of options that will be allowed only on command line
        boost::program_options::options_description genericOptions("Generic options");
        genericOptions.add_options()
                ( "version", "Program version" )
                ( "help,h", "Produce extended program help message" )
                ( "tracta,a", boost::program_options::value<std::string> (&tractAfilename), "filename of first tract")
                ( "tractb,b", boost::program_options::value<std::string> (&tractBfilename), "filename of second tract")
                ( "threshold,t", boost::program_options::value<float> (&thresValue)->implicit_value(0), "[opt] threshold to apply before dissimilarity computation. Default 0 (no threshold)")
                ;

        // Declare a group of options that will be allowed both on command line and in config file
        boost::program_options::options_description configOptions("Configuration");
        configOptions.add_options()
                ( "vista", "[opt] use vista file format (default is nifti)." )
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
        //posOpt.add("rest", -1);


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
            std::cout << "tractdist" << std::endl << std::endl;
            std::cout << "Compute the distance (dissimilarity) value between two tractograms read from file." << std::endl << std::endl;
            std::cout << " --version:       Program version." << std::endl << std::endl;
            std::cout << " -h --help:       produce extended program help message." << std::endl << std::endl;
            std::cout << " -a --tracta:     Filename of first tractogram." << std::endl << std::endl;
            std::cout << " -b --tractb:     Filename of second tractogram." << std::endl << std::endl;
            std::cout << "[-t --threshold]: Threshold to apply directly to the tractogram values before computing the dissimilarity (in order to avoid tractography noise affect the result)." << std::endl;
            std::cout << "                   Unlike in other hClusterign commands, this threshold value is an absolute value to apply to the tractogram data as is, not a relative threshold." << std::endl;
            std::cout << "                   Valid values: [0,1) Use a value of 0 (default) if no thresholding is desired." << std::endl << std::endl;
            std::cout << "[--vista]:        Read/write vista (.v) files [default is nifti (.nii) and compact (.cmpct) files]." << std::endl << std::endl;
            std::cout << std::endl;
            std::cout << "example:" << std::endl << std::endl;
            std::cout << "tractdist -a tractA.nii -b tractB.nii -t 0.2" << std::endl << std::endl;
            exit(0);
        }

        if (variableMap.count( "version" ) )
        {
            std::cout << progName <<", version 2.0"<<std::endl;
            exit(0);
        }

        if (variableMap.count( "tracta" ) )
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path( tractAfilename ) ) )
            {
                std::cerr << "ERROR: tract a file \""<<tractAfilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "tract a  file: "<< tractAfilename << std::endl;
        }
        else
        {
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }
        if ( variableMap.count( "tractb" ) )
        {
            if(!boost::filesystem::is_regular_file(boost::filesystem::path( tractBfilename ) ) )
            {
                std::cerr << "ERROR: tract b file \""<<tractBfilename<<"\" is not a regular file"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "tract a  file: "<< tractBfilename << std::endl;
        }
        else
        {
            std::cerr << visibleOptions << std::endl;
            exit(-1);
        }


        if (variableMap.count( "threshold" ) )
        {
            if(thresValue < 0 || thresValue >= 1)
            {
                std::cerr << "ERROR: tthreshold must be [0,1)"<<std::endl;
                std::cerr << visibleOptions << std::endl;
                exit(-1);
            }
            std::cout << "Tractogram threshold: "<< thresValue << std::endl;
        }
        else
        {
            std::cout << "No tractogram threshold will be applied" << std::endl;
        }



        if ( variableMap.count( "vista" ) )
        {
            std::cout << "Using vista format" << std::endl;
            fileManagerFactory fmf;
            fmf.setVista();
        }
        else
        {
            std::cout << "Using nifti format" << std::endl;
            fileManagerFactory fmf;
            fmf.setNifti();
        }

        // ========== OBTAIN DISTANCES ==========

        fileManagerFactory fileMF;
        fileManager& fileMgr( fileMF.getFM() );
        fileMgr.readAsLog();
        fileMgr.readAsUnThres();
        compactTract tractA, tractB;
        compactTractChar charTractA, charTractB;
        fileMgr.readTract(tractAfilename, &tractA);
        fileMgr.readTract(tractAfilename, &charTractA);
        fileMgr.readTract(tractBfilename, &tractB);
        fileMgr.readTract(tractBfilename, &charTractB);



        double directDist(0);
        tractA.threshold(thresValue);
        tractB.threshold(thresValue);
        charTractA.threshold(thresValue);
        charTractB.threshold(thresValue);
        tractA.computeNorm();
        tractB.computeNorm();
        charTractA.computeNorm();
        charTractB.computeNorm();

        directDist = tractA.tractDistance(tractB);
        std::cout<<"Direct distance:\t"<< string_utils::toString(directDist) <<std::endl;

        if( charTractA.size() > 0 && charTractB.size() > 0 )
        {
            double charDist(0);
            charDist = charTractA.tractDistance(charTractB);
            std::cout<<"Direct (char-char) distance:\t"<<string_utils::toString(charDist)<<std::endl;
        }


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

