#include "get_names.h"

// "get_Vname()": returns the filename of the tractogram corresponding to the given seed voxel coordinates
std::string get_Vname (std::string path, size_t seed) {


    std::stringstream stream;
    std::string number;
    char buffer [7];
    sprintf (buffer, "%06d", seed);
    stream << buffer;
    stream >> number;

    std::string filename(path + "/compact_" + number + ".v");

    return filename;
}// end "get_Vname()" -----------------------------------------------------------------





