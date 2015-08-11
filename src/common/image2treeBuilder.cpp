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



// std library
#include <vector>
#include <list>
#include <string>
#include <map>
#include <utility>
#include <algorithm>

#include <boost/filesystem.hpp>

#include "WStringUtils.h"
#include "image2treeBuilder.h"


#define DEBUG false

image2treeBuilder::image2treeBuilder( std::string imageFilename, std::string treeFilename, bool verbose, std::string baseFilename ):
    m_partTreeReady(false), m_logfile( 0 ), m_tree(treeFilename), m_verbose( verbose )
{
    fileManagerFactory fMFtestFormat;
    m_niftiMode = fMFtestFormat.isNifti();

    m_treeLoaded = (m_tree.isLoaded());
    m_tree.m_treeName = boost::filesystem::path( imageFilename ).stem().string();

    m_roi = m_tree.getRoi();
    m_datasetSize = m_tree.getDataSize();
    m_datasetGrid = m_tree.getDataGrid();
    m_tree.m_discarded.clear();
    m_imageLoaded = loadImage( imageFilename );
    m_basesLoaded =  loadBases( baseFilename );
}

bool image2treeBuilder::loadBases( std::string baseFilename )
{
    m_baseVector.clear();
    if( baseFilename.empty() )
    {
        if( !m_tree.testRootBaseNodes() )
        {
            std::cerr << "ERROR @ image2treeBuilder::loadBases: no bases file stated and tree has no valid meta-leaves (leaves-only base-nodes)" << std::endl;
            return false;
        }
        m_baseVector = m_tree.getRootBaseNodes();
    }
    else
    {
        WFileParser parser( baseFilename );

        if( !parser.readFile() )
        {
            std::cerr << "ERROR @ image2treeBuilder::loadBases: Parser error when reading bases file" << std::endl;
            return false;
        }
        std::vector<std::string> lines = parser.getRawLines();
        if( lines.size() == 0 )
        {
            std::cerr << "ERROR @ image2treeBuilder::loadBases: bases file is empty" << std::endl;
            return false;
        }
        {
            std::vector< std::vector< std::string> >baseStrings = parser.getLinesForTagSeparated( "bases" );
            if( baseStrings.empty() )
            {
                std::cerr << "ERROR @ image2treeBuilder::loadBases: no entires found in the #bases tag check bases-file contents and/or format" << std::endl;
                return false;
            }
            m_baseVector.reserve( baseStrings.size() );
            for( size_t i = 0; i < baseStrings.size(); ++i )
            {
                if( baseStrings[i].size() != 1 )
                {
                    std::cerr << "ERROR @ image2treeBuilder::loadBases: multiple base IDs in the same line, check  bases-file format" << std::endl;
                    return false;
                }
                m_baseVector.push_back( string_utils::fromString<size_t>( baseStrings[i][0] ) );
            }
        }
    }


    if( m_verbose )
    {
        std::cout << "Bases loaded, " << m_baseVector.size() << " base nodes" << std::endl;
    }
    return true;
} // end treeBuilder::loadBases() -------------------------------------------------------------------------------------


bool image2treeBuilder::loadImage( std::string imageFilename )
{
    fileManagerFactory fileMF;
    fileManager& fManager(fileMF.getFM());

    std::vector<std::vector<std::vector<float> > > partImage;
    ValueType imageValueType = fManager.readImage( imageFilename, &partImage );

    if (imageValueType==VTError)
    {
        std::cerr<< "ERROR @ image2treeBuilder::loadImage(): there was an error calling readIamge()" <<std::endl;
        return false;
    }
    if (imageValueType!=VTUINT8)
    {
        std::cerr<< "ERROR @ image2treeBuilder::loadImage(): image must be of type UINT8"<<std::endl;
        return false;
    }

    const size_t dimz   = partImage[0][0].size();
    const size_t dimy    = partImage[0].size();
    const size_t dimx = partImage.size();

    m_partImage.clear();

    {
        std::vector<size_t> zvect(dimz,0);
        std::vector<std::vector<size_t> >yzmatrix(dimy,zvect);
        m_partImage.resize(dimx,yzmatrix);
    }
    for (int i=0 ; i<dimx ; ++i)
    {
        for (int j=0 ; j<dimy ; ++j)
        {
            for (int k=0 ; k<dimz ; ++k)
            {
                m_partImage[i][j][k] = partImage[i][j][k];
            }
        }
    }


    if( m_verbose )
        std::cout << "Partition Image loaded, Dimensions" << m_partImage.size() << "x" << m_partImage.front().size() << "x" << m_partImage.front().front().size() << std::endl;
    return true;
} // end treeBuilder::loadImage() -------------------------------------------------------------------------------------

void image2treeBuilder::importImagePart()
{

    //// TEST IF READY ////

    if( !inReady() )
    {
        throw std::runtime_error("ERROR: not all files were properly loaded" );
    }

    //// TEST IF IMAGE DIMENSIONS AND TREE MATCH ////

    if( m_partImage.size() != m_datasetSize.m_x || m_partImage.front().size() != m_datasetSize.m_y || m_partImage.front().front().size() != m_datasetSize.m_z )
    {
        std::cerr << "ERROR: Image and tree dimensions dont match. ";
        std::cerr << "Tree: " << m_datasetSize.m_x << "x" << m_datasetSize.m_y << "x" << m_datasetSize.m_z << ". ";
        std::cerr << "Image: " << m_partImage.size() << "x" << m_partImage.front().size() << "x" << m_partImage.front().front().size() << std::endl;
        throw std::runtime_error( "ERROR: Image and tree dimensions dont match" );
    }

    //// TEST IF NON-0 VOXELS AND ROI SIZE MATCH ////

    size_t voxelCount(0);
    for(size_t i(0); i< m_partImage.size(); ++i)
    {
        for(size_t j(0); j< m_partImage[i].size(); ++j)
        {
            for(size_t k(0); k< m_partImage[i][j].size(); ++k)
            {
                if( m_partImage[i][j][k] != 0 )
                {
                    ++voxelCount;
                }
            }

        }

    }
    if( voxelCount != m_roi.size() )
    {
        std::cerr << "ERROR: Image and tree roi sizes dont match. Tree: " << m_roi.size() << ". Image: " << voxelCount << std::endl;
        throw std::runtime_error( "ERROR: Image and tree roi sizes dont match" );
    }


    if (!m_baseVector.empty())
    {
        //// PRUNE ALL NON BASE NODES AND RESET BASES TO 0.1 ////

        for(size_t i(0); i < m_tree.getNumNodes(); ++i)
        {
            m_tree.fetchNode( i )->setFlag( true );
        }

        for(size_t i(0); i < m_baseVector.size(); ++i)
        {
            WHnode*  thisNode(m_tree.fetchNode( m_baseVector[i] ) );
            thisNode->setFlag( false );
            thisNode->setDistLevel(0.1);
        }

        m_tree.fetchRoot()->setFlag( false );

        m_tree.cleanup();


        //// ASSIGN EACH BASE TO THAT OF LABEL MOST CONTAINED ////

        std::vector<size_t> bases( m_tree.getRootBaseNodes() );
        std::vector<size_t> bestLabels;
        bestLabels.reserve( bases.size() );

        for(size_t i(0); i < bases.size(); ++i)
        {
            // obtain voxels for each base node
            std::vector<WHcoord> baseVoxels( m_tree.getCoordinates4node( bases[i] ) );
            std::vector<size_t> labels;
            std::vector<size_t> labelCount;
            labels.reserve( baseVoxels.size() );
            labelCount.reserve( baseVoxels.size() );

            // obtain labels for each voxel
            for(size_t j(0); j < baseVoxels.size(); ++j)
            {
                size_t thisLabel( m_partImage[baseVoxels[j].m_x][baseVoxels[j].m_y][baseVoxels[j].m_z] );
                if( thisLabel == 0 )
                {
                    std::cerr << "ERROR: at i="<< i << " j="<< j <<"  Label is 0 for voxel " << baseVoxels[j].m_x << "," << baseVoxels[j].m_y << "," << baseVoxels[j].m_z << std::endl;
                    throw std::runtime_error( "ERROR: Label is 0 for this voxel" );
                }

                std::vector<size_t>::iterator findLabelIter( std::find( labels.begin(), labels.end(), thisLabel) );
                if( findLabelIter == labels.end() )
                {
                    labels.push_back( thisLabel );
                    labelCount.push_back( 1 );
                }
                else
                {
                    size_t pos( findLabelIter - labels.begin() );
                    labelCount[pos] = labelCount[pos] + 1;
                }
            }
            if (labels.size() == 1)
            {
                bestLabels.push_back( labels.front() );
            }
            else
            {
                // remove labels of value 1 (no juelich label)
                for(size_t j(0); j < labels.size();)
                {
                    if ( labels[j] == 1 )
                    {
                        labels.erase( labels.begin() + j );
                        labelCount.erase( labelCount.begin() + j );
                    }
                    else
                    {
                        ++j;
                    }
                }

                // assingn each voxel to highest occurence label
                std::vector<size_t>::iterator findBestLabel( std::max_element( labelCount.begin(), labelCount.end() ) );
                size_t pos( findBestLabel - labelCount.begin() );
                bestLabels.push_back( labels[pos]  );
            }
        }

        //// JOIN BASE NODES BY LABEL AT 0.5 ////

        std::vector<bool> doneLabelFlags( bestLabels.size() ,false );

        std::vector<size_t> labelNodes;

        m_tree.m_nodes.erase(m_tree.m_nodes.end()-1);

        for(size_t i(0); i < bestLabels.size(); ++i)
        {
            // get nodes for that label
            if(doneLabelFlags[i])
            {
                continue;
            }
            size_t thisLabel( bestLabels[i] );

            std::vector<size_t> nodes4label;
            nodes4label.reserve( bestLabels.size() );
            nodes4label.push_back( bases[i] );

            for(size_t j(i+1); j < bestLabels.size(); ++j)
            {
                if( bestLabels[j] == thisLabel )
                {
                    nodes4label.push_back( bases[j] );
                    doneLabelFlags[j] = true;
                }
                else
                {
                    continue;
                }
            }
            doneLabelFlags[i] = true;

            // create new node for the label
            size_t labelNodeID( m_tree.m_nodes.size() );
            dist_t labelNodeDist( 0.5 );
            size_t labelNodeHlevel( 2 );
            size_t labelNodeSize(0);
            std::vector<nodeID_t> labelNodeChildren;

            for(size_t j(0); j < nodes4label.size(); ++j)
            {
                labelNodeChildren.push_back( std::make_pair( true, nodes4label[j] ) );
                labelNodeSize += m_tree.m_nodes[ nodes4label[j] ].getSize();
                m_tree.m_nodes[ nodes4label[j] ].setParent( std::make_pair( true, labelNodeID ) );
            }

            WHnode labelNode( std::make_pair( true, labelNodeID ), labelNodeChildren, labelNodeSize, labelNodeDist, labelNodeHlevel );
            m_tree.m_nodes.push_back( labelNode );

            labelNodes.push_back( labelNodeID );
        }

        //// JOIN LABEL NODES AT 1 ////

        // create new root
        size_t rootNodeID( m_tree.m_nodes.size() );
        dist_t rootNodeDist( 1 );
        size_t rootNodeHlevel( 3 );
        size_t rootNodeSize( m_tree.m_leaves.size() );
        std::vector<nodeID_t> rootNodeChildren;

        for(size_t i(0); i < labelNodes.size(); ++i)
        {
            rootNodeChildren.push_back( std::make_pair( true, labelNodes[i] ) );
            m_tree.m_nodes[ labelNodes[i] ].setParent( std::make_pair( true, rootNodeID ) );
        }

        WHnode rootNode( std::make_pair( true, rootNodeID ), rootNodeChildren, rootNodeSize, rootNodeDist, rootNodeHlevel );
        m_tree.m_nodes.push_back( rootNode );

        m_tree.check();
    }
    else
    {

        //erase root node
        m_tree.m_nodes.clear();

        std::vector<size_t> labels;
        std::vector<bool> donelabels;
        labels.reserve( m_tree.m_leaves.size());
        donelabels.reserve( m_tree.m_leaves.size());

        for(size_t i(0); i < m_tree.m_leaves.size(); ++i)
        {

            size_t thisLabel( m_partImage[m_tree.m_coordinates[i].m_x][m_tree.m_coordinates[i].m_y][m_tree.m_coordinates[i].m_z] );
            if( thisLabel == 0 )
            {
                throw std::runtime_error( "ERROR: Label is 0 for this voxel" );
            }
            labels.push_back( thisLabel );
            donelabels.push_back( false );
        }






    }

    ////// WRITE TREE //////////

    m_partTreeReady = true;

    if( m_verbose )
    {
        std::cout << m_tree.getReport() << std::endl;
    }
    if( m_logfile != 0 )
    {
        ( *m_logfile ) << m_tree.getReport() << std::endl;
    }

    m_tree.m_treeName = ( "partitionTree" );
    writeTree();
}


void image2treeBuilder::writeTree() const
{
    if( ( !m_partTreeReady ) || m_outputFolder.empty() )
    {
        std::cerr << "ERROR @ treeBuilder::writeTree(): Tree is not ready, or outputfolder is not set" << std::endl;
        return;
    }

    m_tree.writeTree( m_outputFolder + "/" + m_tree.m_treeName + ".txt", m_niftiMode );
        //    m_tree.writeTreeDebug( m_outputFolder + "/" + m_tree.m_treeName + "_debug.txt" );
    //    m_tree.writeTreeOldWalnut(m_outputFolder+"/"+m_tree.m_treeName+"_4ow.txt");

    if( m_verbose )
    {
        std::cout << "Written standard tree file in: " << m_outputFolder << "/" << m_tree.m_treeName << ".txt" << std::endl;
        //        std::cout << "Written debug tree file in: " << m_outputFolder << "/" << m_tree.m_treeName << "_debug.txt" << std::endl;
        //        std::cout<<"Written walnut tree file in: "<< m_outputFolder <<"/"<<m_tree.m_treeName<<"_4ow.txt"<<std::endl;
    }
    if( m_logfile != 0 )
    {
        ( *m_logfile ) << "Standard tree file in:\t" << m_outputFolder << "/" << m_tree.m_treeName << ".txt" << std::endl;
        //        ( *m_logfile ) << "Debug tree file in:\t" << m_outputFolder << "/" << m_tree.m_treeName << "_debug.txt" << std::endl;
        //        (*m_logfile)<<"Walnut tree file in:\t"<< m_outputFolder <<"/"<<m_tree.m_treeName<<"_4ow.txt"<<std::endl;
    }

    return;
} // end treeBuilder::writeTree() -------------------------------------------------------------------------------------
