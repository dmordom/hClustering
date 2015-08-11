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

#ifndef FILEMANAGERFACTORY_H
#define FILEMANAGERFACTORY_H

// std library
#include <string>

// hClustering
#include "fileManager.h"
#include "niftiManager.h"
#include "vistaManager.h"

/**
 * This class handles the creation of an appropiate file manager working either nifti or vista format files depending on the value of a static member
 */
class fileManagerFactory
{
public:
    /**
     * Constructor
     * \param ioFolderInit input-output folder for the file manager
     */
    fileManagerFactory( std::string ioFolderInit=""): m_niftiMngr(ioFolderInit), m_vistaMngr(ioFolderInit) {}

    //! Destructor
    ~fileManagerFactory() {}

    // === PUBLIC MEMBER FUNCTIONS ===

    /**
      * queries whether the nifti mode is currently set
      * \return if true current file mode is nifti
      */
    bool isNifti() { return m_isNifti; }

    /**
      * queries whether the vista mode is currently set
      * \return if true current file mode is vista
      */
    bool isVista() { return !m_isNifti; }

    /**
      * prepares the class to return a reference to a nifti file manager
      */
    void setNifti() { m_isNifti=true; }

    /**
      * prepares the class to return a reference to a vista file manager
      */
    void setVista() { m_isNifti=false; }

    /**
      * returns a refererence to the appropiate file manager
      * \return reference to a fileManger object
      */
    fileManager& getFM();

private:
    // === PRIVATE DATA MEMBERS ===

    static bool m_isNifti; //!< flag to store whether we are in vista or nifti mode
    niftiManager m_niftiMngr; //!< file manager to work with nifti files
    vistaManager m_vistaMngr;  //!< file manager to work with vista files
};

#endif // FILEMANAGERFACTORY_H
