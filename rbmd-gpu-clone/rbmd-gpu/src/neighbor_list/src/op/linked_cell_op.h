#pragma once
#include "common/device_types.h"
#include "common/types.h"
#include "linked_cell/linked_cell.h"

namespace op {
template <typename DEVICE>
struct InitializeCellOp {
  void operator()(LinkedCellDeviceDataPtr* linked_cell, Box* box, Cell* cells,
                  rbmd::Id total_cells);
};

template <typename DEVICE>
struct AssignAtomsToCellOp {
  void operator()(rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, Box* d_box,
                  LinkedCellDeviceDataPtr* linked_cell, Cell* cells,
                  rbmd::Id* per_atom_cell_id, rbmd::Id total_atoms_num);
};

template <typename DEVICE>
struct ComputeCellRangesIndicesOp {
  void operator()(rbmd::Id* sorted_cell_index,
                  rbmd::Id* d_in_atom_list_start_index,
                  rbmd::Id* d_in_atom_list_end_index, rbmd::Id num_atoms);
};

// template <typename DEVICE>
// struct ComputeCellAtomsCountOp {
//   void operator()(Cell* cells,
//                   rbmd::Id* d_in_atom_list_start_index,
//                   rbmd::Id* d_in_atom_list_end_index, rbmd::Id num_cells);;
// };

template <>
struct InitializeCellOp<device::DEVICE_GPU> {
  void operator()(LinkedCellDeviceDataPtr* linked_cell, Box* box, Cell* cells,
                  rbmd::Id total_cells);
};

template <>
struct AssignAtomsToCellOp<device::DEVICE_GPU> {
  void operator()(rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, Box* d_box,
                  LinkedCellDeviceDataPtr* linked_cell, Cell* cells,
                  rbmd::Id* per_atom_cell_id, rbmd::Id total_atoms_num);
};

template <>
struct ComputeCellRangesIndicesOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* sorted_cell_index,
                  rbmd::Id* d_in_atom_list_start_index,
                  rbmd::Id* d_in_atom_list_end_index, rbmd::Id num_atoms);
};

// template <>
// struct ComputeCellAtomsCountOp<device::DEVICE_GPU> {
//   void operator()(Cell* cells,
//                   rbmd::Id* d_in_atom_list_start_index,
//                   rbmd::Id* d_in_atom_list_end_index, rbmd::Id num_cells);
// };

template <typename DEVICE>
struct MapAtomidToIdxOp {
  void operator()(rbmd::Id* d_atomid2idx, rbmd::Id* d_sorted_atom_idx,
                  rbmd::Id num_atoms);
};

template <>
struct MapAtomidToIdxOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* d_atomid2idx, rbmd::Id* d_sorted_atom_idx,
                  rbmd::Id num_atoms);
};

}  // namespace op
