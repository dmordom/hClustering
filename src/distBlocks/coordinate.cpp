#include "coordinate.h"

// member functions

// "phys_dist()": returns the euclidean distance between this voxel and the input voxel
float coordinate::phys_dist(const coordinate voxel) const{
    double temp ( sqrt( pow((x-voxel.x),2.) + pow((y-voxel.y),2.) + pow((z-voxel.z),2.) ) );
    float phys_dist(temp);
    return phys_dist;
}// end "phys_dist()" -----------------------------------------------------------------


// "get_phys_neighbours()": returns a vector containing the coordinates of the physical neighbours adjacent to the voxel, the neighbourhood level is defined by nb_level
std::vector<coordinate> coordinate::get_phys_neighbours(const coordinate maxdim, const int nb_level) const {

    std::vector<coordinate> phys_neighbours;
    coordinate temp;

    // NBhood level 6
    if (nb_level == 6) {
        for (coord_t i=-1; i<=1 ; i+=2) {
            temp=*this;
            temp.z = this->z + i;
            if ( !((temp.z<0)||(temp.y<0)||(temp.x<0)||(temp.z>maxdim.z)||(temp.y>maxdim.y)||(temp.x>maxdim.x))) // voxel is not out of mask boundaries
                phys_neighbours.push_back(temp);

            temp=*this;
            temp.y = this->y + i;
            if ( !((temp.z<0)||(temp.y<0)||(temp.x<0)||(temp.z>maxdim.z)||(temp.y>maxdim.y)||(temp.x>maxdim.x))) // voxel is not out of mask boundaries
                phys_neighbours.push_back(temp);

            temp=*this;
            temp.x = this->x + i;
            if ( !((temp.z<0)||(temp.y<0)||(temp.x<0)||(temp.z>maxdim.z)||(temp.y>maxdim.y)||(temp.x>maxdim.x))) // voxel is not out of mask boundaries
                phys_neighbours.push_back(temp);
        }
    }

    // NBhood level 18
    if (nb_level == 18) {
        for (coord_t i=-1; i<=1 ; ++i) { // loop through basic cube
            for (coord_t j=-1; j<=1 ; ++j) {
                for (coord_t k=-1; k<=1 ; ++k) {
                    if ( ((i==-1)||(i==1)) && ((j==-1)||(j==1)) && ((k==-1)||(k==1)) )  // dont use the cube corners as neighbours
                        continue;
                    if (i==0 && j==0 && k==0) // dont use the voxel itself as neighbour
                        continue;
                    temp.z = this->z + i;
                    temp.y = this->y + j;
                    temp.x = this->x + k;
                    if ( (temp.z<0)||(temp.y<0)||(temp.x<0)||(temp.z>maxdim.z)||(temp.y>maxdim.y)||(temp.x>maxdim.x)) // voxel is out of mask boundaries
                        continue;
                    temp.subject = this->subject;
                    phys_neighbours.push_back(temp);
                }
            }
        }
    }

    // NBhood level 26
    if (nb_level == 26) {
        for (coord_t i=-1; i<=1 ; ++i) {
            for (coord_t j=-1; j<=1 ; ++j) {
                for (coord_t k=-1; k<=1 ; ++k) {
                    if (i==0 && j==0 && k==0) // dont use the voxel itself as neighbour
                        continue;
                    temp.z = this->z + i;
                    temp.y = this->y + j;
                    temp.x = this->x + k;
                    if ( (temp.z<0)||(temp.y<0)||(temp.x<0)||(temp.z>maxdim.z)||(temp.y>maxdim.y)||(temp.x>maxdim.x)) // voxel is out of mask boundaries
                        continue;
                    temp.subject = this->subject;
                    phys_neighbours.push_back(temp);
                }
            }
        }
    }

    // NBhood level 32
    if (nb_level == 32) {
        for (coord_t i=-1; i<=1 ; ++i) {
            for (coord_t j=-1; j<=1 ; ++j) {
                for (coord_t k=-1; k<=1 ; ++k) {
                    if (i==0 && j==0 && k==0) // dont use the voxel itself as neighbour
                        continue;
                    temp.z = this->z + i;
                    temp.y = this->y + j;
                    temp.x = this->x + k;
                    if ( (temp.z<0)||(temp.y<0)||(temp.x<0)||(temp.z>maxdim.z)||(temp.y>maxdim.y)||(temp.x>maxdim.x)) // voxel is out of mask boundaries
                        continue;
                    temp.subject = this->subject;
                    phys_neighbours.push_back(temp);
                }
            }
        }

        for (coord_t i=-2; i<=2 ; i+=2) {
            for (coord_t j=-2; j<=2 ; j+=2) {
                for (coord_t k=-2; k<=2 ; k+=2) {
                    if (i==0 && j==0 && k==0) // dont use the voxel itself as neighbour
                        continue;
                    if ( ((i!=0)&&(j!=0)) || ((j!=0)&&(k!=0)) || ((k!=0)&&(i!=0)) )
                        continue;
                    temp.z = this->z + i;
                    temp.y = this->y + j;
                    temp.x = this->x + k;
                    if ( (temp.z<0)||(temp.y<0)||(temp.x<0)||(temp.z>maxdim.z)||(temp.y>maxdim.y)||(temp.x>maxdim.x)) // voxel is out of mask boundaries
                        continue;
                    temp.subject = this->subject;
                    phys_neighbours.push_back(temp);
                }
            }
        }
    }



    return phys_neighbours;
}// end "phys_dist()" -----------------------------------------------------------------


// "get_name_string()": returns a string with the coordinates of the voxel in the form "xxx_yyy_zzz"
std::string coordinate::get_name_string() const {

    std::stringstream stream;
    std::string name;
    stream << x <<"_"<< y <<"_"<< z;
    stream >> name;
    return name;

}// end "phys_dist()" -----------------------------------------------------------------


// member operators
coordinate& coordinate::operator =(const coordinate &rhs) {

    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    subject = rhs.subject;
    return *this;
}// end "operator =" -----------------------------------------------------------------


// non-member operators

std::ostream& operator <<(std::ostream& os, const coordinate& object) {

    std::string x_fill;
    if(object.x < 100)
        x_fill += "0";
    if(object.x < 10)
        x_fill += "0";
    std::string y_fill(" ");
    if(object.y < 100)
        y_fill += "0";
    if(object.y < 10)
        y_fill += "0";
    std::string z_fill(" ");
    if(object.z < 100)
        z_fill += "0";
    if(object.z < 10)
        z_fill += "0";
    std::string s_fill(" ");
    if(object.subject < 100)
        s_fill += "0";
    if(object.subject < 10)
        s_fill += "0";

    os << x_fill << object.x << y_fill << object.y << z_fill << object.z;
//    os << x_fill << object.x << y_fill << object.y << z_fill << object.z << s_fill << object.subject;
    return os;
}// end "operator <<" -----------------------------------------------------------------


bool operator <(const coordinate &lhs, const coordinate &rhs) {
    if (lhs.z < rhs.z) // z1<z2
        return true;
    else if (lhs.z > rhs.z) // z1>z2
        return false;
    else if (lhs.y < rhs.y) // z1=z2, y1<y2
        return true;
    else if (lhs.y > rhs.y) // z1=z2, y1>y2
        return false;
    else if (lhs.x < rhs.x) // z1=z2, y1=y2, x1<x2
        return true;
    else if (lhs.x > rhs.x) // z1=z2, y1=y2, x1>x2
        return false;
    else if (lhs.subject < rhs.subject) // z1=z2, y1=y2, x1=x2 sub1<sub2
        return true;
    else if (lhs.subject > rhs.subject) // z1=z2, y1=y2, x1=x2 sub1>sub2
        return false;
    else  // coord1 = coord2
        return false;
}// end "operator <" -----------------------------------------------------------------


bool operator ==(const coordinate &lhs, const coordinate &rhs) {
    if (lhs.z == rhs.z)
        if (lhs.y == rhs.y)
            if (lhs.x == rhs.x)
                return true;

    return false;
}// end "operator ==" -----------------------------------------------------------------
