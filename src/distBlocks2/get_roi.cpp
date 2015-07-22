
#include "get_roi.h"

void read_tree(std::string roi_filename,
               std::vector<coordinate> &roivect,
               std::map<coordinate,size_t> &roimap) {



    // Read ROI file
    std::ifstream roi_file(roi_filename.c_str(), std::ifstream::in);
    std::string roi_line;
    size_t coord_count(0);
    coord_t temp, x_temp, y_temp, z_temp;

    if (!roi_file)
        throw std::runtime_error ("ERROR [read_tree()]: unable to open ROI file");

    while (getline(roi_file, roi_line)) {
        int dim_count(1);
        std::istringstream stream(roi_line);

        while(stream >> temp) {
            switch (dim_count) {
                case 1:
                    x_temp = temp;
                    break;
                case 2:
                    y_temp = temp;
                    break;
                case 3:
                    z_temp = temp;
                    break;
                default:
                    throw std::runtime_error ("ERROR [read_tree()]: too many elements per row in roi file.\nformat of ROI file must be x y z numeric coordinates per line separated by whitespaces");
            }
            ++dim_count;
        }

        if(stream.bad())
            throw std::runtime_error ("ERROR [read_tree()]: roi stream corrupted");
        if(stream.fail() && !stream.eof())
            throw std::runtime_error ("ERROR [read_tree()]: roi data may not be correctly formatted\nformat of ROI file must be x y z numeric coordinates per line separated by whitespaces");
        if (dim_count != 4)
            throw std::runtime_error ("ERROR [read_tree()]: too few elements per row in roi file.\nformat of ROI file must be x y z numeric coordinates per line separated by whitespaces");
        coordinate coord_temp(x_temp,y_temp,z_temp, 0);
        roimap.insert(std::make_pair(coord_temp,coord_count++)); //Insert new coordinate
        roivect.push_back(coord_temp);
    }
    roi_file.close();



    return;

}
