#pragma once
#include "device_types.h"
#include "types.h"
#include "model/box.h"
namespace op {

    template <typename DEVICE>
    struct ShakeAOp {
      void operator()(const rbmd::Id num_angle, 
                      const rbmd::Real dt,
                      const rbmd::Real fmt2v, 
                      Box* box,
                      const rbmd::Id* atom_id_to_idx,
                      const rbmd::Real* mass, 
                      const rbmd::Id* atoms_type,
                      const Id3* angle_id_vec,
                      rbmd::Real* shake_px,
                      rbmd::Real* shake_py,
                      rbmd::Real* shake_pz,
                      rbmd::Real* shake_vx,
                      rbmd::Real* shake_vy,
                      rbmd::Real* shake_vz,
                      const rbmd::Real* fx,
                      const rbmd::Real* fy, 
                      const rbmd::Real* fz, 
                      rbmd::Id* flag_px,
                      rbmd::Id* flag_py, 
                      rbmd::Id* flag_pz);
    };

    template <typename DEVICE>
    struct ShakeBOp {
        void operator()(rbmd::Id num_angle,
                        rbmd::Real dt,
                        rbmd::Real fmt2v,
                        Box* box,
                        const rbmd::Id* atom_id_to_idx,
                        const rbmd::Real* mass,
                        const rbmd::Id* atoms_type,
                        const Id3* angle_id_vec,
                        rbmd::Real* px,
                        rbmd::Real* py,
                        rbmd::Real* pz,
                        rbmd::Real* shake_vx,
                        rbmd::Real* shake_vy,
                        rbmd::Real* shake_vz,
                        const rbmd::Real* fx,
                        const rbmd::Real* fy,
                        const rbmd::Real* fz);
    };

    template <>
    struct ShakeAOp<device::DEVICE_GPU> {
      void operator()(rbmd::Id num_angle,
                      rbmd::Real dt,
                      rbmd::Real fmt2v,
                      Box* box,
                      const rbmd::Id* atom_id_to_idx,
                      const rbmd::Real* mass,
                      const rbmd::Id* atoms_type,
                      const Id3* angle_id_vec,
                      rbmd::Real* shake_px,
                      rbmd::Real* shake_py,
                      rbmd::Real* shake_pz,
                      rbmd::Real* shake_vx,
                      rbmd::Real* shake_vy,
                      rbmd::Real* shake_vz,
                      const rbmd::Real* fx,
                      const rbmd::Real* fy,
                      const rbmd::Real* fz,
                      rbmd::Id* flag_px,
                      rbmd::Id* flag_py,
                      rbmd::Id* flag_pz);
    };

    template <>
    struct ShakeBOp<device::DEVICE_GPU>  {
        void operator()(const rbmd::Id num_angle,
                        const rbmd::Real dt,
                        const rbmd::Real fmt2v,
                        Box* box,
                        const rbmd::Id* atom_id_to_idx,
                        const rbmd::Real* mass,
                        const rbmd::Id* atoms_type,
                        const Id3* angle_id_vec,
                        rbmd::Real* px,
                        rbmd::Real* py,
                        rbmd::Real* pz,
                        rbmd::Real* shake_vx,
                        rbmd::Real* shake_vy,
                        rbmd::Real* shake_vz,
                        const rbmd::Real* fx,
                        const rbmd::Real* fy,
                        const rbmd::Real* fz);
    };
}  // namespace op