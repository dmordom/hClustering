#ifndef COORDINATE_H
#define COORDINATE_H

// std library
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

// typedefs
typedef short coord_t;

//Class to contain a seed voxel coordinate. consisting on x,y,z values and a subject number
class coordinate
{
public:

    // constructors
    coordinate(): x(0), y(0), z(0), subject(0) {}
    coordinate(coord_t x_init, coord_t y_init, coord_t z_init, unsigned short subject_init): x(x_init), y(y_init), z(z_init), subject(subject_init) {}
    coordinate(const coordinate &object): x(object.x), y(object.y), z(object.z), subject(object.subject) {}
    ~coordinate() {}

    // member operators
    coordinate& operator =(const coordinate &rhs);

    // member functions

    // "phys_dist()": returns the euclidean distance between this voxel and the input voxel
    float phys_dist(const coordinate voxel) const;

    // "get_phys_neighbours()": returns a vector containing the coordinates of the physical neighbours adjacent to the voxel, the neighbourhood level is defined by nb_level
    std::vector<coordinate> get_phys_neighbours(const coordinate maxdim, const int nb_level) const;

    // "get_name_string()": returns a string with the coordinates of the voxel in the form "xxx_yyy_zzz"
    std::string get_name_string() const;

    // data members
    coord_t x;
    coord_t y;
    coord_t z;
    unsigned short subject;
};

// non-member operators
std::ostream&   operator <<(std::ostream& os, const coordinate& object);
bool            operator < (const coordinate &lhs, const coordinate &rhs);
bool            operator ==(const coordinate &lhs, const coordinate &rhs);


#endif // COORDINATE_H
