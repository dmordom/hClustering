//---------------------------------------------------------------------------
//
// Project: hCLustering
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
// hCLustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------

#ifndef TREECOMPARER_H
#define TREECOMPARER_H

// parallel execution
#include <omp.h>

// std library
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <sstream>

// boost library
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

// hClustering
#include "compactTract.h"
#include "WHcoord.h"
#include "WHtree.h"
#include "listedCache.hpp"
#include "fileManagerFactory.h"
#include "treeManager.h"
#include "WHtreeProcesser.h"


#define TREE1 true
#define TREE2 false

/**
 * this class implements the main functionalities required for performing meta-leaf matching across trees and calculating tree comparison algorithms across matched-trees
 */
class treeComparer
{
public:
    /**
     * Constructor
     * \param tree1 a pointer to the first tree being matched/compared
     * \param tree2 a pointer to the second tree being matched/compared
     * \param verbose the verbose output flag
     */
    treeComparer( WHtree* const tree1, WHtree* const tree2, bool verbose );
    /**
     * Constructor
     * \param tree1 a pointer to the first tree being matched/compared
     * \param tree2 a pointer to the second tree being matched/compared
     * \param comparer the treeComparer object with the datamembers to be copied
     */
    treeComparer( WHtree* const tree1, WHtree* const tree2, treeComparer &comparer );

    //! Destructor
    ~treeComparer() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * Sets the folder path to find leaf tracts corresponding to tree 1
     * \param singleTractFolder1 the path to the folder with tree 1 leaf tracts
     */
    inline void setSingleTractFolder1( std::string singleTractFolder1 ) { m_singleTractFolder1 = singleTractFolder1; }

    /**
     * Sets the folder path to find leaf tracts corresponding to tree 2
     * \param singleTractFolder2 the path to the folder with tree 2 leaf tracts
     */
    inline void setSingleTractFolder2( std::string singleTractFolder2 ) { m_singleTractFolder2 = singleTractFolder2; }

    /**
     * Sets the folder path to find node tracts corresponding to tree 1
     * \param meanTractFolder1 the path to the folder with tree 1 node tracts
     */
    inline void setMeanTractFolder1( std::string meanTractFolder1 ) { m_meanTractFolder1 = meanTractFolder1; }

    /**
     * Sets the folder path to find node tracts corresponding to tree 2
     * \param meanTractFolder2 the path to the folder with tree 2 node tracts
     */
    inline void setMeanTractFolder2( std::string meanTractFolder2 ) { m_meanTractFolder2 = meanTractFolder2; }

    /**
     * Sets the maximum euclidean distance allowed between centres of matched meta-leaves (across trees). If a negative value is entered it will be set to 0
     * \param maxPhysDist the max euclidean distance allowed between matched meta-leaf centres
     */
    inline void setMaxPhysDist( size_t maxPhysDist ) { m_maxPhysDist = ( maxPhysDist > 0 ) ? maxPhysDist : 0; }

    /**
     * sets output file stream for the program log file
     * \param logfile a pointer to the output log file stream
     */
    inline void log( std::ofstream* const logfile ) { m_logfile = logfile; }

    /**
     * when set the meta-leaves contained seed voxel coordinates will be read from cluster mask files instead of computed
     * \param coordsFromFile the boolean value to set the m_coordsFromFile flag to
     */
    inline void setCoordsFromFile( bool coordsFromFile ) { m_coordsFromFile = coordsFromFile; }

    /**
     * when set the meta-leaves mean tracts will be read from files instead of computed
     * \param meanTractsFromFile the boolean value to set the m_meanTractsFromFile flag to
     */
    inline void setMeanTractsFromFile( bool meanTractsFromFile ) { m_meanTractsFromFile = meanTractsFromFile; }

    /**
     * returns the flag indicating whether both loaded trees have valid base-nodes (meta-leaves) with only seed leaf children
     * \return true if both trees have valid meta-leaves, false otherwise
     */
    inline bool areRealBaseNodes() const { return m_realBaseNodes; }


    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * sets the threshold (relative to the number of streamlines that were generated per seed voxel on the leaf tracts) that will be applied to
     * meta-leaf tracts before computing similarity for leaf matching (in order to avoid tracking noise/artifacts to affect similarities).
     * normalized thershold value will be automatically calculated.
     * \param relThres the relative threshold value to use
     */
    void setRelativeThreshold( float relThres );

    /**
     * Computes the dissimilarity matrix from the meta-leaves mean-tracts of one tree to those of the other tree
     * mean tracts should have been previosuly computed and properly transformed to a common space
     */
    void getBaseDistMatrix();

    /**
     * Saves a precomputed cross-tree meta-leaves mean-tracts dissimilarity matrix to file
     * \param matrixFilename the filepath to save the matrix data to
     */
    void readBaseDistMatrix( std::string matrixFilename );

    /**
     * Reads a previously saved cross-tree meta-leaves mean-tracts dissimilarity matrix from file
     * \param matrixFilename the filepath to read the matrix data from
     */
    void writeBaseDistMatrix( std::string matrixFilename );

    /**
     * Computes the tree cpcc comparison value between two matched trees whose meta-leaf dissimilarity matrix has been previosuly computed.
     */
    std::pair< std::pair< float, float >, std::pair< float, float > > doTcpcc() const;

    /**
     * Computes the simple triplets comparison value between two matched trees whose meta-leaf dissimilarity matrix has been previosuly computed.
     * \param sampleFreq a sample frequency in case not all the points wish to be analyzed (the process  can be quite time consuming).
     *        default value (1) uses all possible points and provides the real exact value. Higher frequencies provide faster computation at the price of approximated value.
     */
    std::pair< float, float > simpleTriplets( size_t sampleFreq = 1 ) const;

    /**
     * Pefrorm leaf-wise correspondence between trees and get a matching table as output.
     * Use only when trees have been built on the same set of seed voxel tractograms
     */
    bool leafCorrespondence();

    /**
     * Pefrorm base-nodes greedy correspondence between trees and get a matching table as output.
     * \param dissimThreshold A dissimilarity limit for the matching, pairs that match with a dissimilairty higher than the threshold will be regarded as having no valid match
     * \param redoCoords if set base-node clusters mean-coordinate information will be recomputed from cluster mask files
     */
    void greedyCorrespondence( float dissimThreshold, bool redoCoords = true );

    /**
     * Pefrorm a random base-node matching between trees (euclidean cluster distance restrictions will still be applied).
     * used in order to obtain a baseline for the comparison algorithm values to asses robustness and meaningfulness
     */
    void randomCorrespondence();

    /**
     * Write correspondence table (matchign of base-nodes across trees) to file
     * Base nodes will be identified with relative values, that is, base node position in the base node vector
     * \param outFilename the filepath to the write the table to
     */
    void writeCorrespondence( std::string outFilename );

    /**
     * Write correspondence table (matchign of base-nodes across trees) to file
     * Base nodes will be identified with absolute values values, that is, node IDs within the tree structure
     * \param outFilename the filepath to the write the table to
     */
    void writeFullCorrespondence( std::string outFilename );

    /**
     * Performs an analyiss of the base node matching and gives several informational outputs to evaluate the quality of the matching
     * Among them specially the matchign quality defined as the size-normalized tractogram similarity between matched nodes averaged across all basenodes
     * \return a vector with the matching evaluation outputs in this order: mean matchign dissimilarity; weighted matchign dissimilarity (matchign quality);
     *         fraction of total leaves matched; mean euclidean distance between matched centers; weighted euclidean dist between matched centers.
     */
    std::vector<float> rateCorrespondence();

    /**
     * Calculates the noise level of the trees defined as the matching distance of each base node to the corresponding one of the matched tree
     * this "noise level" is then used to eliminate the hierarchical structure below that level (all granularities higher than that) under the assumption that
     * there is insufficient matching quality to use that information for tree comparison
     * \param noiseAlpha noise alpha coefficient. Controls the weight of the matching quality on the noise level.
     *        When 0 no noise restriction will be applied, when 1 the npoise level will be the averaged matchting quality level of the contained base nodes.
     * \return a pair structutre with the number of clusters remaining at the highest granularity level for each tree after application of the noise simulation
     */
    std::pair< float, float > applyNoiseBaseline( const float noiseAlpha = 0 );

    /**
     * Analyzes the number and size of the base nodes of each tree and returns a report message with the results
     * \return the report message
     */
    std::string reportBaseNodes() const;

    /**
     * Loads the base node ID information of both trees into the corresponding data members and calculates the mean coordinate of each base-node clusterif so indicated
     * \param doGetCoords tif true mean coordinates for each base node cluster will be computed from cluster mask averaging
     */
    bool fetchBaseNodes(bool doGetCoords = true);


private:
    // === PRIVATE DATA MEMBERS ===

    WHtree &m_tree1;                    //!< A reference to the first tree to be matched and compared
    WHtree &m_tree2;                    //!< A reference to the second tree to be matched and compared

    std::string m_singleTractFolder1;   //!< the folder path to find leaf tractograms from tree 1 seed voxels
    std::string m_singleTractFolder2;   //!< the folder path to find leaf tractograms from tree 2 seed voxels
    std::string m_meanTractFolder1;     //!< the folder path to find node tractograms from tree 1 warped to common space
    std::string m_meanTractFolder2;     //!< the folder path to find node tractograms from tree 2 warped to common space
    float m_maxPhysDist;                //!< the maximum euclidean distance between matched cluster centres to be allowed as a match
    float m_tractThreshold1;            //!< the threshold to apply to normalized tracts of tree 1 before computing tract dissimilarity (to avoid noise artifacts)
    float m_tractThreshold2;            //!< the threshold to apply to normalized tracts of tree 2 before computing tract dissimilarity (to avoid noise artifacts)

    std::ofstream *m_logfile;           //!< A pointer to the output log file stream

    std::vector< size_t > m_baseNodes1;         //!< the base node list from tree 1 after tree matching
    std::vector< size_t > m_originalBaseNodes1; //!< the original base node list from tree 1 before tree matching
    std::vector< WHcoord > m_baseCoords1;       //!< the mean coordinate of each matched base node cluster of tree 1
    std::vector< float > m_noiseLevels1;        //!< the matching noise level of each matched base node of tree 1

    std::vector< size_t > m_baseNodes2;         //!< the base node list from tree 2 after tree matching
    std::vector< size_t > m_originalBaseNodes2; //!< the original base node list from tree 2 before tree matching
    std::vector< WHcoord > m_baseCoords2;       //!< the mean coordinate of each matched base node cluster of tree 2
    std::vector< float > m_noiseLevels2;        //!< the matching noise level of each matched base node of tree 2

    std::pair< size_t, size_t > m_initialSizes; //!< original number of leaves for each tree before matching (for comparison computation purposes, unmatched leaves are eliminated after matching)
    std::vector< std::vector< dist_t > > m_baseDistMatrix;  //!< the distance matrix between base node tracts of tree 1 and tree 2
    bool m_realBaseNodes;                       //!< if true both trees have valid base nodes
    bool m_coordsFromFile;                      //! if true cluster coordinate masks will be read from file (instead of computed on the fly)
    bool m_meanTractsFromFile;                  //! if true mean node tractograms will be read from file (instead of computed on the fly)
    bool m_verbose;                             //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.

    std::vector< size_t > m_fullCorrespondence; //!< correspondence table between matched base nodes of each tree using absolute node IDs (within each trees tructure)
    std::vector< size_t > m_newCorrespondence;  //!< updated relative correspondence table between matched base nodes across trees after elimination of unmatched nodes (IDs are relative to base node positiioon in base node vector)
    std::vector< size_t > m_newCorrespReverse;  //! the correspondence table described below but in the reverse direction (from tree 2 base nodes to tree 1 base nodes)
    std::vector< std::pair< float, float > > m_correspDistances; //!< correspondence distances (tractogram distance, euclidean distance of cluster centres)


    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * Load or compute the seed voxels contained in each base node cluster of each tree and calculate the average coordinate from all of them, saving the results in the appropiate data members
     */
    void getBaseCoords();

    /**
     * Apply the noise correction (structure simplification below noise level) for the desired tree
     * \param treeCode a boolean value identifying whether to apply the correction to tree 1 or tree 2
     * \param noiseAlpha noise alpha coefficient. Controls the weight of the matching quality on the noise level.
     *        When 0 no noise restriction will be applied, when 1 the npoise level will be the averaged matchting quality level of the contained base nodes.
     * \return the highest granularity remainign after noise correction (maximum number of clusters in the structure)
     */
    size_t noiseBaseline( const bool treeCode, const float noiseAlpha );

    /**
     * Returns the relative noise ID of a base node (position in the base node vector)
     * \param absoluteID the absolute ID of the base node in the tree structure
     * \param baseNodes a reference to the vector containing the base node IDs
     * \return the highest granularity remainign after noise correction (maximum number of clusters in the structure)
     */
    size_t findRelativeBasenodeID( size_t absoluteID, const std::vector< size_t > &baseNodes ) const;

    /**
     * Finds the minimum euclidean distance between a coordinate and a set of coordinates
     * (in practice between the mean coordinate of a cluster and the seed voxel coordinates of a second cluster)
     * \param thisCoord the coordinate to compute euclidean distance from
     * \param nodeCoords the set of seed voxels from which to find the closest one to the coordinate and calculate euclidean distance from
     * \return the esulting euclidean distance value
     */
    float minPhysDist( const WHcoord thisCoord, std::vector< WHcoord > nodeCoords );
};

#endif  // TREECOMPARER_H
