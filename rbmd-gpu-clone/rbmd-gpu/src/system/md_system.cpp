#include "md_system.h"

#include "base/locator.h"
#include "base/memory/memory_op.h"
#include "common/device_types.h"
#include "near_force/direct_truncation/direct_truncation_op.h"

int MDSystem::Evolve() {
  auto& structure_info = _md_data._structure_info;
  auto& potential_data = _md_data._potential_data;
  auto& structure_data = _md_data._structure_data;
  auto& nAtoms = structure_info._num_atoms;
  rbmd::Real *d_dt = nullptr, *d_fmt2v = nullptr, *d_mass = nullptr;
  rbmd::Real3 *d_v = nullptr, *d_force = nullptr, *d_position;
  rbmd::Id3* d_cellid;
  Locator* locator = new Locator();
  Locator* d_locator;

  rbmd::Real dt = 0.5, fmt2v = 1.0;
  rbmd::Real3* force = new rbmd::Real3[nAtoms];
  op::resize_memory_op<rbmd::Real, device::DEVICE_GPU>()(d_dt, 1);
  op::resize_memory_op<rbmd::Real, device::DEVICE_GPU>()(d_fmt2v, 1);
  op::resize_memory_op<rbmd::Real, device::DEVICE_GPU>()(
      d_mass, potential_data._mass.size());
  op::resize_memory_op<rbmd::Real3, device::DEVICE_GPU>()(d_v, nAtoms);
  op::resize_memory_op<rbmd::Real3, device::DEVICE_GPU>()(d_force, nAtoms);
  op::resize_memory_op<rbmd::Real3, device::DEVICE_GPU>()(d_position, nAtoms);
  op::resize_memory_op<rbmd::Id3, device::DEVICE_GPU>()(d_cellid, nAtoms);
  op::resize_memory_op<Locator, device::DEVICE_GPU>()(d_locator, 1);

  op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>()(d_dt, &dt, 1);
  op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>()(d_fmt2v, &fmt2v, 1);
  op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>()(
      d_mass, potential_data._mass.data(), 1);
  op::sync_memory_h2d_op<rbmd::Real3, device::DEVICE_GPU>()(
      d_v, structure_data._velocities.data(), nAtoms);
  op::sync_memory_h2d_op<rbmd::Real3, device::DEVICE_GPU>()(d_force, force,
                                                            nAtoms);
  op::sync_memory_h2d_op<rbmd::Real3, device::DEVICE_GPU>()(
      d_position, structure_data._positions.data(), nAtoms);
  op::sync_memory_h2d_op<Locator, device::DEVICE_GPU>()(d_locator, locator, 1);

  op::direct_truncation_op<rbmd::Real, device::DEVICE_GPU>()(
      1, structure_info._num_atoms, rbmd::Real3(0, 0, 0),
      rbmd::Real3(10, 10, 10), rbmd::Id3(10, 10, 10), d_cellid, d_dt, d_fmt2v,
      d_mass, d_locator, d_position, d_v, d_force);

  // op::delete_memory_op<rbmd::Real, device::DEVICE_GPU>(d_dt);
  // op::delete_memory_op<rbmd::Real, device::DEVICE_GPU>(d_fmt2v);
  // op::delete_memory_op<rbmd::Real, device::DEVICE_GPU>(d_mass);
  // op::delete_memory_op<rbmd::Real3, device::DEVICE_GPU>(d_v);
  // op::delete_memory_op<rbmd::Real3, device::DEVICE_GPU>(d_force);

  return 0;
}

int MDSystem::PreSolve() { return 0; }

int MDSystem::Solve() { return 0; }

int MDSystem::PostSolve() { return 0; }
