#ifndef GRAPHTREEBUILDER_H
#define GRAPHTREEBUILDER_H

        // parallel execution
        #include <omp.h>

        // std library
        #include <vector>
        #include <list>
        #include <string>
        #include <map>
        #include <fstream>
        #include <ctime>
        #include <climits>
        #include <sstream>

        // boost library
        #include <boost/filesystem.hpp>

        // hClustering
        #include "WHcoord.h"
        #include "distBlock.h"
        #include "WFileParser.h"
        #include "WHtree.h"
        #include "vistaManager.h"


        typedef enum
        {
            TG_SINGLE, TG_COMPLETE, TG_AVERAGE, TG_WEIGHTED
        } tgraph_t;

        class graphTreeBuilder
        {
        public:

            explicit graphTreeBuilder( std::string roiFilename, bool verbose = true):
                m_roiLoaded( false ), m_treeReady( false ), m_logfile( 0 ), m_verbose( verbose )
            {
                readRoi( roiFilename );
            }

            ~graphTreeBuilder()
            {
            }

            void log( std::ofstream* const logfile )
            {
                m_logfile = logfile;
            }
            bool ready() const
            {
                return m_roiLoaded;
            }
            void setInputFolder( std::string inputFolder )
            {
                m_inputFolder = inputFolder;
            }
            void setOutputFolder( std::string outputFolder )
            {
                m_outputFolder = outputFolder;
            }
            size_t roiSize() const
            {
                return m_roi.size();
            }

            bool readRoi( const std::string roiFilename );
            void buildGraph( const tgraph_t graphMethod );
            void writeTree() const;

        private:
            bool m_verbose;
            bool m_roiLoaded;
            bool m_treeReady;
            WHtree m_tree;
            std::vector< WHcoord > m_roi;

            WHcoord m_datasetSize;
            HC_GRID m_datasetGrid;
            std::ofstream *m_logfile;
            std::string m_inputFolder;
            std::string m_outputFolder;

            WHnode* fetchNode( const nodeID_t &thisNode, std::vector< WHnode > &leaves, std::vector< WHnode > &nodes ) const;

            void loadDistMatrix( std::vector< std::vector< float > >* const distMatrix ) const;

            dist_t newGraphDist( const dist_t &distance1, const dist_t &distance2, const size_t &size1, const size_t &size2,
                            const tgraph_t &graphMethod ) const;
        };

#endif // GRAPHTREEBUILDER_H
