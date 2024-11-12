#include "Atom.h"
#include <ciso646>

Atom::Atom(const unsigned long id,
           const double rx,
           const double ry,
           const double rz,
           const double vx,
           const double vy,
           const double vz,
           const double charge)
{
  _id = id;
  _r[0] = rx;
  _r[1] = ry;
  _r[2] = rz;
  _v[0] = vx;
  _v[1] = vy;
  _v[2] = vz;
  _charge = charge;
}

unsigned long Atom::GetID() const { return _id; }

void Atom::SetID(unsigned long id) { _id = id; }

Real_mpi Atom::r(unsigned short dimension) const { return _r[dimension]; }

vtkm::Vec3f Atom::Position3D() const
{
    vtkm::Vec3f positon ={ static_cast<vtkm::FloatDefault>(_r[0]),static_cast<vtkm::FloatDefault>(_r[1]),static_cast<vtkm::FloatDefault>(_r[2]) };
    return positon;
}

vtkm::Vec3f Atom::Velocity3D() const
{
    vtkm::Vec3f velocity = { static_cast<vtkm::FloatDefault>(_v[0]), static_cast<vtkm::FloatDefault>(_v[1]), static_cast<vtkm::FloatDefault>(_v[2]) };
    return velocity;
}

void Atom::Setr(unsigned short dimension, Real_mpi r) { _r[dimension] = r; }

Real_mpi Atom::v(unsigned short dimension) const { return _v[dimension]; }

void Atom::Setv(unsigned short dimension, Real_mpi v) { _v[dimension] = v; }

bool Atom::InBox(const double *bounding_box_min, const double *bounding_box_max) const {
    bool in_box = true;
#ifdef __INTEL_COMPILER
    #pragma float_control(precise, on)
			#pragma fenv_access(on)
#endif
    for (int dimension = 0; dimension < 3; ++dimension) {
        in_box &= (r(dimension) >= bounding_box_min[dimension] and r(dimension) < bounding_box_max[dimension]);
    }

    return in_box;
}

void Atom::setF(unsigned short d, double F)  { _F[d] = F; }

double Atom::F(unsigned short d) const  {return _F[d]; }

std::array<double, 3> Atom::F_arr() {
    std::array<double, 3> f_ret{F(0), F(1), F(2)};
    return f_ret;
}



