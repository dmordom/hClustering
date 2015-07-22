#include "output.h"

// "write_output():" Writes the output files (roi_block_index)
void write_output(std::string path,
                  std::vector<std::pair<size_t,size_t> > &roi_block_index) {



    // ========== Write roi_block_index file ==========

    // Write new seed voxel mask that corresponds to the tree file
    std::string roi_block_index_filename(path+"/roi_index.txt");
    std::cout<< "Writing roi_block_index in \""<< roi_block_index_filename <<"\""<< std::endl;
    FILE * roi_block_index_file;
    roi_block_index_file = fopen (roi_block_index_filename.c_str(),"w");
    ;

    fprintf(roi_block_index_file,"#distindex\n");

    for( size_t i = 0; i < roi_block_index.size(); ++i )
    {
        fprintf(roi_block_index_file,"%04d b %03d i %04d\n",i, roi_block_index[i].first,roi_block_index[i].second);
    }

    fprintf(roi_block_index_file,"#enddistindex");


    fclose(roi_block_index_file);



}
