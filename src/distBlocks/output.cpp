#include "output.h"

// "write_output():" Writes the output files (roi_block_index)
void write_output(std::string path,
                  std::vector<coordinate> &roivect,
                  std::vector<std::pair<size_t,size_t> > &roi_block_index) {



    // ========== Write roi_block_index file ==========

    // Write new seed voxel mask that corresponds to the tree file
    std::string roi_block_index_filename(path+"/roi_index.txt");
    std::cout<< "Writing roi_block_index in \""<< roi_block_index_filename <<"\""<< std::endl;
    FILE * roi_block_index_file;
    roi_block_index_file = fopen (roi_block_index_filename.c_str(),"w");
    std::vector<std::pair<size_t,size_t> >::const_iterator index_iter=roi_block_index.begin();

    fprintf(roi_block_index_file,"#distindex\n");

    for( std::vector<coordinate>::const_iterator coord_iter=roivect.begin(); coord_iter != roivect.end() ; ++coord_iter, ++index_iter)
        fprintf(roi_block_index_file,"%03d %03d %03d b %03d i %04d\n",coord_iter->x,coord_iter->y,coord_iter->z, index_iter->first, index_iter->second);

    fprintf(roi_block_index_file,"#enddistindex");


    fclose(roi_block_index_file);



}
