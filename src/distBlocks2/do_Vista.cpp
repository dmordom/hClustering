#include "do_Vista.h"

// "get_Vtract_th()": extracts tractogram data from a vista file and returns it thresholded in float array form
void get_Vtract_th (std::string filename, unsigned char *tractogram, size_t tract_length, float tractThreshold) {

    // ========== Initialize variables ==========

    VAttrList list;
    VAttrListPosn position;
    VImage tractImage(NULL);
    size_t n_elements(0), count(0);
    unsigned char tract_value_char;
    bool is_char(true);


    // ========== Read Image ==========


        //#pragma omp critical (file)  // this piece of code is never executed simultaneously by several threads to prevent errors
        {
            FILE *cor_file(NULL);
            cor_file = VOpenInputFile ((char *)filename.c_str(), TRUE);
            list = VReadFile (cor_file, NULL);
            if (!cor_file) {
                fclose (cor_file);
                VError ("get_Vtract_th(): Failed to open input file '%s'", (char *)filename.c_str());
            }
            fclose(cor_file);
        }


        if (!list) {
            VError ("get_Vtract_th(): Failed to read input file '%s'", (char *)filename.c_str());
        }


    // Verify attributes
    for (VFirstAttr (list, & position); VAttrExists (& position); VDeleteAttr  (& position)) {
        if (VGetAttrRepn (& position) != VImageRepn) continue;
        VGetAttrValue (& position, NULL, VImageRepn, & tractImage);
        ++count;
        if (VPixelRepn(tractImage) == VUByteRepn)
            is_char = true;
        else if (VPixelRepn(tractImage) == VFloatRepn)
            is_char = false;
        else
            VError("get_Vtract_th(): Error: tractogram image must be of type char or float");
        if (VImageNBands(tractImage)!=1 && VImageNRows(tractImage)!= 1)
            VError("get_Vtract_th(): Error: tractogram image must have 1 row and 1 band only");
        n_elements = VImageNColumns(tractImage);
    }

    if (count>1)
        VError("get_Vtract_th(): Error: tractogram file has more than one image");

    // ========== Write vector tract ==========
    if (n_elements!=tract_length)
        VError("get_Vtract_th(): Error: actual tractogram size is different from input value");


    if (is_char){  //it is a seed tractogram in char format

        unsigned char int_thres = (unsigned char) (tractThreshold*255.);

        // Copy tractogram data
        for (int i=0; i<n_elements; ++i) {
            tract_value_char = VPixel(tractImage,0,0,i,VUByte);

            if ((tractThreshold != 0) && (tract_value_char < int_thres) )
                tractogram[i]=0;
            else
                tractogram[i]=(tract_value_char);
        }
    }
    else { //it is a node tractogram in float format
        VError("get_Vtract_th(): Error: tractogram is float");
    }

    // clean up
    VDestroyAttrList (list);
    VDestroyImage(tractImage);

    return;

} // end "get_Vtract_th()" -----------------------------------------------------------------



// "get_Vtract()": extracts tractogram data from a vista file and returns it in std vector form
std::vector<char> get_Vtract (std::string filename) {

    // ========== Initialize variables ==========

    std::vector<char> tractOut;
    VAttrList list;
    VAttrListPosn position;
    VImage tractImage(NULL);
    int n_elements(0), count(0);
    unsigned char tract_value_char;
    bool is_char(true);


    // ========== Read Image ==========


        //#pragma omp critical (file)  // this piece of code is never executed simultaneously by several threads to prevent errors
        {
            FILE *cor_file(NULL);
            cor_file = VOpenInputFile ((char *)filename.c_str(), TRUE);
            list = VReadFile (cor_file, NULL);
            if (!cor_file) {
                fclose (cor_file);
                VError ("get_Vtract(): Failed to open input file '%s'", (char *)filename.c_str());
            }
            fclose(cor_file);
        }


        if (!list) {
            VError ("get_Vtract(): Failed to read input file '%s'", (char *)filename.c_str());
        }


    // Verify attributes
    for (VFirstAttr (list, & position); VAttrExists (& position); VDeleteAttr  (& position)) {
        if (VGetAttrRepn (& position) != VImageRepn) continue;
        VGetAttrValue (& position, NULL, VImageRepn, & tractImage);
        ++count;
        if (VPixelRepn(tractImage) == VUByteRepn)
            is_char = true;
        else if (VPixelRepn(tractImage) == VFloatRepn)
            is_char = false;
        else
            VError("get_Vtract(): Error: tractogram image must be of type char or float");
        if (VImageNBands(tractImage)!=1 && VImageNRows(tractImage)!= 1)
            VError("get_Vtract(): Error: tractogram image must have 1 row and 1 band only");
        n_elements = VImageNColumns(tractImage);
    }

    if (count>1)
        VError("get_Vtract(): Error: tractogram file has more than one image");

    // ========== Write vector tract ==========

    if (is_char){  //it is a seed tractogram in char format

        // Copy tractogram data
        tractOut.reserve(n_elements); // Reserve memory in vector to increase speed
        for (int i=0; i<n_elements; ++i) {
            tract_value_char = VPixel(tractImage,0,0,i,VUByte);
            tractOut.push_back(tract_value_char);
        }
    }
    else { //it is a node tractogram in float format

        VError("get_Vtract_th(): Error: tractogram is float");

    }

    // clean up
    VDestroyAttrList (list);
    VDestroyImage(tractImage);

    return tractOut;

} // end "get_Vtract()" -----------------------------------------------------------------



// "write_dist_block()": writes down distance block to file
void write_dist_block (std::string filename, std::vector<std::vector<float> > &dist_block) {

    const int rows(dist_block.size());
    const int columns(dist_block[0].size());

    VImage v_dist_block = VCreateImage(1,rows,columns,VFloatRepn);

    for (std::vector<float>::size_type row_count(0); row_count<rows; ++row_count) {

        for (std::vector<float>::size_type col_count(0); col_count<columns; ++col_count) {

            VPixel(v_dist_block,0,row_count,col_count,VFloat) = dist_block[row_count][col_count];

        }
    }

    // write tractogram to file
    if ( ! WriteVImage( (char *)filename.c_str() , v_dist_block ) )
        VError( "write_Vtract(): Failed to open output tractogram file " );

    // clean up
    VDestroyImage( v_dist_block );

    return;

} // end "write_dist_block()" -----------------------------------------------------------------

// "write_Vtract()": write compact tractogram data into vista format file
void write_Vtract (std::string filename, std::vector<float> &tractogram) {


    const int columns(tractogram.size());

    VImage vtract = VCreateImage(1,1,columns,VFloatRepn);

    std::vector<float>::size_type v_count(0);
    for (std::vector<float>::const_iterator iter=tractogram.begin(); iter != tractogram.end() ; ++iter, ++v_count)
        VPixel(vtract,0,0,v_count,VFloat) = *iter;

    // write tractogram to file
    if ( ! WriteVImage( (char *)filename.c_str() , vtract ) )
        VError( "write_Vtract(): Failed to open output tractogram file " );

    // clean up
    VDestroyImage( vtract );

    return;

} // end "write_Vtract()" -----------------------------------------------------------------



// "WriteVImage()": write Vista image
int WriteVImage( const VString Name, VImage &Image )
{

    VAttrListPosn pos;       /* position in list */
    VBoolean      success;   /* success flag     */
    VAttrList list;

    //#pragma omp critical (file)  // this piece of code is never executed simultaneously by several threads to prevent errors
   {
       FILE*  file;      /* output file      */
       list = VCreateAttrList();
       VAppendAttr( list, "image", NULL, VImageRepn, Image );
       /* open file */
       file = fopen (Name, "w");
       if (!file) {
           VError ("WriteVImage(): Failed to open output vista file ");
       }
       /* write file */
       success = VWriteFile (file, list);
       fclose (file);
   }

   if (!success)
   { VError ("WriteVImage(): Failed to write output file '%s'", Name);
   }

   /* remove images */
   for (VFirstAttr (list, &pos); VAttrExists (&pos); VNextAttr (&pos))
      if (VGetAttrRepn (&pos) == VImageRepn)
         VSetAttrValue (&pos, NULL, VImageRepn, NULL);

   /* clean-up*/
   VDestroyAttrList (list);

    return 1;

} // end "WriteVImage()" -----------------------------------------------------------------


