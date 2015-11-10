# hClustering
Hierarchical clustering algorithms for whole-brain connectivity similarity analysis
 
 David Moreno-Dominguez. d.mor.dom@gmail.com, moreno@cbs.mpg.de, www.cbs.mpg.de/alumni/moreno-11266

The hClustering repository groups a set of algorithms created to generate fully hierarchical characterization of human whole-brain high-res dMRI-based anatomical connectivity data (although with minor source code modifications it could be used with functional or any other kind of data that can be represented as a “fingerprint” vector).

These algorithms are delivered in the form of easy-to-use Linux command-line tools with comprehensive command help information. The source code is fully documented using Doxygen, and several additional manuals are included such as a quick-start guide and descriptions of the custom file-formats used.

The original code was developed in the framework of the Whole-Brain Hierarchical clustering project carried out in the Max Planck Institute for Human Cognitive and Brain Sciences in Leipzig, as part of the PhD Thesis:  “Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography” by David Moreno-Dominguez. For this release, the code has been documented, upgraded and streamlined to be more easy to use and accept both the neuroimaging de-facto standard Nifti (.nii) and Lipsia neuroimaging processing package Vista (.v) formats.

* Doxygen-generated documentation can be found in the "doc" directory (open html/index.html file for browsing).
* User guides can be found in the "manuals" folder (in .pdf format).

For further information on the underlying algorithms and research done with this code refer to:

 - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014).
   A hierarchical method for whole-brain connectivity-based parcellation.
   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528

 - Moreno-Dominguez, D. (2014).
   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography.
   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig.
   ISBN 978-3-941504-45-5.



###	Prior requirements

Install if necessary the following library packages:

- Boost [libboost-all]:	Tested with version 1.48, higher versions should also work.
- via [libviaio]:		     Provided as part of the repository. For vista files handling.
- nifti-1 [libniftiio]: For nifti files handling.
- lznz [libznz]: 		     Required by nifti library.

Clone and compile the hClustering  code from GitHub: 
https://github.com/dmordom/hClustering.git.

The minimum data required to use the hClustering commands is listed below:
- 	A 3D white matter mask (the one used as tracking target space).
- 	A set of probabilistic tractograms in 1D compact form (as many elements as voxels in the white matter mask). vdconnect tractography directly outputs tractograms in compact form (in vista format). If available tracts are 3D nifti images they can be transformed to compact form with the full2compact command from this repository.
-	An ASCII roi file with the seed coordinates, image size and streamline number information (see file format guide for details).
These (and more) data can be automatically generated by the dMRI preprocessing pipeline available at https://github.com/dmordom/dmri_prepro_nipype.git. For more information refer to the pipeline repository readme and the quick-start guide.

Visualization and interactive exploration of the hierarchical tree files generated with this repository is possible through the "Hierarchical Clustering" module developed and integrated in the OpenWalnut framework. This program can be downloaded from www.openwalnut.com. For more information refer to the manuals and the OpenWalnut website.

###____________________________________________

 Whole-Brain Connectivity-Based Hierarchical Parcellation Project
 
 David Moreno-Dominguez
 d.mor.dom@gmail.com
 moreno@cbs.mpg.de
 www.cbs.mpg.de/alumni/moreno-11266

 hClustering is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 http://creativecommons.org/licenses/by-nc/3.0

 hClustering is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.


