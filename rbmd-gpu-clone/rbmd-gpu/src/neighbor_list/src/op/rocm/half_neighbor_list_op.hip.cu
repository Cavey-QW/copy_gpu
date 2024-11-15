#include <hip/hip_runtime.h>

#include <hipcub/hipcub.hpp>

#include "common/device_types.h"
#include "common/rbmd_define.h"
#include "common/types.h"
#include "half_neighbor_list_op.h"
#include "model/box.h"
#include "src/op/full_neighbor_list_op.h"

namespace op {
__global__ void ComputeHalfNeighbors(rbmd::Id* per_dimension_cells,
                                     rbmd::Id* neighbor_cell,
                                     rbmd::Id neighbor_num, rbmd::Id total_cell,
                                     rbmd::Id cell_count_within_cutoff) {
  const unsigned int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell_idx < total_cell) {
    rbmd::Id idx_z =
        cell_idx / (per_dimension_cells[0] * per_dimension_cells[1]);
    rbmd::Id idx_y =
        (cell_idx - idx_z * per_dimension_cells[0] * per_dimension_cells[1]) /
        per_dimension_cells[0];
    rbmd::Id idx_x = cell_idx - per_dimension_cells[0] *
                                    (idx_y + per_dimension_cells[1] * idx_z);
    int neighbor_count = 0;
    for (int dz = -cell_count_within_cutoff; dz <= cell_count_within_cutoff;
         ++dz) {
      // TODO cell in cutoff
      // RBL 2！
      for (int dy = -cell_count_within_cutoff; dy <= cell_count_within_cutoff;
           ++dy) {
        for (int dx = -cell_count_within_cutoff; dx <= cell_count_within_cutoff;
             ++dx) {
          const int offset =
              (dz * per_dimension_cells[1] + dy) * per_dimension_cells[0] + dx;
          if (offset >= 0) {
            int neighbor_z =
                (idx_z + dz + per_dimension_cells[2]) % per_dimension_cells[2];
            int neighbor_y =
                (idx_y + dy + per_dimension_cells[1]) % per_dimension_cells[1];
            int neighbor_x =
                (idx_x + dx + per_dimension_cells[0]) % per_dimension_cells[0];

            rbmd::Id neighbour_cell_idx =
                (neighbor_z * per_dimension_cells[1] + neighbor_y) *
                    per_dimension_cells[0] +
                neighbor_x;

            neighbor_cell[cell_idx * neighbor_num + neighbor_count] =
                neighbour_cell_idx;
            ++neighbor_count;
          }
        }
      }
    }
  }
}

/// 一个原子对应一个线程束   建议调用的时候 -EPSILON
__global__ void EstimateHalfNeighborList(
    rbmd::Id* __restrict__ per_atom_cell_id,
    rbmd::Id* __restrict__ in_atom_list_start_index,
    rbmd::Id* __restrict__ in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* __restrict__ px,
    rbmd::Real* __restrict__ py, rbmd::Real* __restrict__ pz,
    rbmd::Id* __restrict__ neighbour_num,
    rbmd::Id* __restrict__ max_neighbour_num, Box* __restrict__ box,
    rbmd::Id* __restrict__ neighbor_cell, rbmd::Id neighbor_cell_num,
    bool without_pbc_or_rbl) {
  // cutoff2是平方
  extern __shared__ hipcub::WarpReduce<rbmd::Id>::TempStorage
      reduce_temp_storage[];
  rbmd::Id atom_neighbor_num = MIN_NBNUM;
  __shared__ rbmd::Real shared_px[BLOCK_SIZE];
  __shared__ rbmd::Real shared_py[BLOCK_SIZE];
  __shared__ rbmd::Real shared_pz[BLOCK_SIZE];
  // 计算当前线程处理的原子的索引    ---- 当前线程原子就是atom_idx
  unsigned int atom_idx = (blockIdx.x * blockDim.x + threadIdx.x) / warpSize;
  if (atom_idx < total_atom_num) {
    rbmd::Real distance = 0;
    // 计算当前线程在warp中的位置
    unsigned int lane_id = (blockIdx.x * blockDim.x + threadIdx.x) % warpSize;
    atom_neighbor_num = 0;
    // 当前所属的cell
    rbmd::Id cell_idx = __ldg(&per_atom_cell_id[atom_idx]);
    if (threadIdx.x < BLOCK_SIZE) {
      shared_px[threadIdx.x] = __ldg(&px[atom_idx]);
      shared_py[threadIdx.x] = __ldg(&py[atom_idx]);
      shared_pz[threadIdx.x] = __ldg(&pz[atom_idx]);
    }
    __syncthreads();
    for (int i = 0; i < neighbor_cell_num; ++i) {
      rbmd::Id neighbour_cell_idx =
          __ldg(&neighbor_cell[cell_idx * neighbor_cell_num + i]);
      for (rbmd::Id neighbor_atom_idx =
               __ldg(&in_atom_list_start_index[neighbour_cell_idx]) + lane_id;
           neighbor_atom_idx < __ldg(&in_atom_list_end_index[neighbour_cell_idx]);
           neighbor_atom_idx += warpSize) {
        bool valid_neighbor;
        if (cell_idx == neighbour_cell_idx || without_pbc_or_rbl == true) {
          // 如果是同一个cell，只考虑索引大于当前原子的邻居
          valid_neighbor = (neighbor_atom_idx > atom_idx);
        } else {
          // 如果是不同的cell，考虑所有不等于当前原子的邻居
          valid_neighbor = (neighbor_atom_idx != atom_idx);
        }
        if (valid_neighbor) {
          distance = CaculateDistance(
              box, shared_px[threadIdx.x], shared_py[threadIdx.x],
              shared_pz[threadIdx.x], __ldg(&px[neighbor_atom_idx]),
              __ldg(&py[neighbor_atom_idx]), __ldg(&pz[neighbor_atom_idx]));
          if (distance < cutoff_2) {
            atom_neighbor_num++;
          }
        }
      }
    }
    atom_neighbor_num = hipcub::WarpReduce<rbmd::Id>(
                            reduce_temp_storage[threadIdx.x / warpSize])
                            .Sum(atom_neighbor_num);
    if (lane_id == 0) {
      neighbour_num[atom_idx] = atom_neighbor_num;
      // 线程束对齐
      max_neighbour_num[atom_idx] =
          MAX(((rbmd::Id)CEIL(atom_neighbor_num * SAFE_ZONE) + warpSize - 1) /
                  warpSize * warpSize,
              MIN_NBNUM);
    }
  }
}

//! Note : 算子有问题，会多次导致 *should_realloc = true 暂时先用NoWarp
// TODO ！ 考虑可能是 hipcub::WarpScan<int>(temp_storage[threadIdx.x /
// warpSize]).ExclusiveSum(is_neighbor, offset);导致的
// __global__ void GenerateHalfNeighborList(
//     rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
//     rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
//     rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
//     rbmd::Id* max_neighbor_num, rbmd::Id* neighbor_start,
//     rbmd::Id* neighbor_end, rbmd::Id* neighbors, Box* d_box,
//     rbmd::Id* should_realloc, rbmd::Id* neighbor_cell, rbmd::Id
//     neighbor_cell_num, bool without_pbc) {
//   const unsigned int atom_idx =
//       (blockIdx.x * blockDim.x + threadIdx.x) / warpSize;
//   const unsigned int lane_id =
//       (blockIdx.x * blockDim.x + threadIdx.x) % warpSize;
//
//   __shared__ hipcub::WarpScan<int>::TempStorage temp_storage[WARP_SIZE];
//
//   rbmd::Id neighbor_num = 0;
//
//   if (atom_idx < total_atom_num) {
//     neighbor_num = neighbor_start[atom_idx];
//     rbmd::Id cell_idx = per_atom_cell_id[atom_idx];
//
//     for (int i = 0; i < neighbor_cell_num; ++i) {
//       rbmd::Id neighbour_cell_idx =
//           neighbor_cell[cell_idx * neighbor_cell_num + i];
//       rbmd::Id start = in_atom_list_start_index[neighbour_cell_idx];
//       rbmd::Id end = in_atom_list_end_index[neighbour_cell_idx];
//
//       for (rbmd::Id base_idx = start; base_idx < end; base_idx += warpSize) {
//         rbmd::Id neighbor_atom_idx = base_idx + lane_id;
//         int is_neighbor = 0;
//         bool valid_neighbor = false;
//         if (cell_idx == neighbour_cell_idx ||
//             without_pbc == true) {  // TODO 暂未测试
//           // 如果是同一个cell，只考虑索引大于当前原子的邻居
//           valid_neighbor = (neighbor_atom_idx > atom_idx);
//         } else {
//           valid_neighbor = (neighbor_atom_idx != atom_idx);
//         }
//         if (valid_neighbor) {
//           rbmd::Real distance =
//               CaculateDistance(d_box, px[atom_idx], py[atom_idx],
//               pz[atom_idx],
//                                px[neighbor_atom_idx], py[neighbor_atom_idx],
//                                pz[neighbor_atom_idx]);
//
//           if (distance < cutoff_2) {
//             is_neighbor = 1;
//           }
//         }
//
//         int offset;
//         hipcub::WarpScan<int>(temp_storage[threadIdx.x / warpSize])
//             .ExclusiveSum(is_neighbor, offset);
//         if (is_neighbor) {
//           neighbors[neighbor_num + offset] = neighbor_atom_idx;
//         }
//
//         neighbor_num += hipcub::ShuffleIndex<WARP_SIZE, rbmd::Id>(
//             offset + is_neighbor, warpSize - 1,
//             0xffffffff);  // 广播线程束最后一个  half错误会是广播导致的吗
//         // __shfl_sync(0xffffffff, offset + is_neighbor, warpSize - 1);
//         // 兼容dcu
//       }
//     }
//
//     if (lane_id == 0) {
//       neighbor_end[atom_idx] = neighbor_num;
//       rbmd::Id my_total_neighbor_num = neighbor_num -
//       neighbor_start[atom_idx]; if (my_total_neighbor_num >
//       max_neighbor_num[atom_idx] &&
//           atom_idx == 816) {
//         atomicOr(should_realloc,RBMD_TRUE);  // TODO atomicOr
//       }
//     }
//   }
// }

__global__ void GenerateHalfNeighborListNoWarp(
    rbmd::Id* __restrict__ per_atom_cell_id,
    rbmd::Id* __restrict__ in_atom_list_start_index,
    rbmd::Id* __restrict__ in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* __restrict__ px,
    rbmd::Real* __restrict__ py, rbmd::Real* __restrict__ pz,
    rbmd::Id* __restrict__ max_neighbor_num,
    rbmd::Id* __restrict__ neighbor_start, rbmd::Id* __restrict__ neighbor_end,
    rbmd::Id* __restrict__ neighbors, Box* __restrict__ d_box,
    rbmd::Id* __restrict__ should_realloc, rbmd::Id* __restrict__ neighbor_cell,
    rbmd::Id neighbor_cell_num, bool without_pbc) {
  // cutoff2是平方
  __shared__ rbmd::Real shared_px[BLOCK_SIZE];
  __shared__ rbmd::Real shared_py[BLOCK_SIZE];
  __shared__ rbmd::Real shared_pz[BLOCK_SIZE];
  // 计算当前线程处理的原子的索引    ---- 当前线程原子就是atom_idx
  const unsigned int atom_idx = (blockIdx.x * blockDim.x + threadIdx.x);
  rbmd::Real distance = 0;
  rbmd::Id neighbor_num = 0;
  if (atom_idx < total_atom_num) {
    neighbor_num = neighbor_start[atom_idx];  // 变量名不太好，这里没全部变为0
    // 当前所属的cell
    rbmd::Id cell_idx = __ldg(&per_atom_cell_id[atom_idx]);
    if (threadIdx.x < BLOCK_SIZE) {
      shared_px[threadIdx.x] = __ldg(&px[atom_idx]);
      shared_py[threadIdx.x] = __ldg(&py[atom_idx]);
      shared_pz[threadIdx.x] = __ldg(&pz[atom_idx]);
    }
    __syncthreads();
    // 遍历当前原子的所有邻居cell
    for (int i = 0; i < neighbor_cell_num; ++i) {
      rbmd::Id neighbour_cell_idx =
            __ldg(&neighbor_cell[cell_idx * neighbor_cell_num + i]);
      rbmd::Id start = __ldg(&in_atom_list_start_index[neighbour_cell_idx]);
      rbmd::Id end = __ldg(&in_atom_list_end_index[neighbour_cell_idx]);
      for (rbmd::Id neighbor_atom_idx = start; neighbor_atom_idx < end;
           ++neighbor_atom_idx) {
        bool valid_neighbor;
        if (cell_idx == neighbour_cell_idx || without_pbc == true) {
          // 如果是同一个cell，只考虑索引大于当前原子的邻居
          valid_neighbor = (neighbor_atom_idx > atom_idx);
        } else {
          // 如果是不同的cell，考虑所有不等于当前原子的邻居
          valid_neighbor = (neighbor_atom_idx != atom_idx);
        }
        if (valid_neighbor) {
          distance = CaculateDistance(
              d_box, shared_px[threadIdx.x], shared_py[threadIdx.x],
              shared_pz[threadIdx.x], __ldg(&px[neighbor_atom_idx]),
              __ldg(&py[neighbor_atom_idx]), __ldg(&pz[neighbor_atom_idx]));
          if (distance < cutoff_2) {
            neighbors[neighbor_num] = neighbor_atom_idx;
            ++neighbor_num;
          }
        }
      }
    }
    neighbor_end[atom_idx] = neighbor_num;
    rbmd::Id my_total_neighbor_num = neighbor_num - neighbor_start[atom_idx];
    if (my_total_neighbor_num > max_neighbor_num[atom_idx]) {
      atomicOr(should_realloc, RBMD_TRUE);
    }
  }
}

void ComputeHalfNeighborsOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
    rbmd::Id neighbor_num, rbmd::Id total_cell,
    rbmd::Id cell_count_within_cutoff) {
  unsigned int blocks_per_grid = (total_cell + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(ComputeHalfNeighbors<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      per_dimension_cells, neighbor_cell, neighbor_num, total_cell,
      cell_count_within_cutoff));
}

void EstimateHalfNeighborListOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
    rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* neighbour_num, rbmd::Id* max_neighbour_num, Box* box,
    rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num,
    bool without_pbc_or_rbl) {
  unsigned int threads_per_block = BLOCK_SIZE;
  unsigned int warps_per_block = threads_per_block / WARP_SIZE;
  unsigned int blocks_per_grid =
      (total_atom_num + warps_per_block - 1) / warps_per_block;
  CHECK_KERNEL(
      EstimateHalfNeighborList<<<blocks_per_grid, BLOCK_SIZE,
                                 sizeof(hipcub::WarpReduce<int>::TempStorage) *
                                     (BLOCK_SIZE / WARP_SIZE),
                                 0>>>(
          per_atom_cell_id, in_atom_list_start_index, in_atom_list_end_index,
          cutoff_2, total_atom_num, px, py, pz, neighbour_num,
          max_neighbour_num, box, neighbor_cell, neighbor_cell_num,
          without_pbc_or_rbl));
}

void GenerateHalfNeighborListOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
    rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* max_neighbor_num, rbmd::Id* neighbor_start,
    rbmd::Id* neighbor_end, rbmd::Id* neighbors, Box* d_box,
    rbmd::Id* should_realloc, rbmd::Id* neighbor_cell,
    rbmd::Id neighbor_cell_num, bool without_pbc) {
  // unsigned int threads_per_block = BLOCK_SIZE;
  // unsigned int warps_per_block = threads_per_block / WARP_SIZE;
  // unsigned int blocks_per_grid =
  //     (total_atom_num + warps_per_block - 1) / warps_per_block;
  unsigned int blocks_per_grid = (total_atom_num + BLOCK_SIZE - 1) / BLOCK_SIZE;

  CHECK_KERNEL(
      GenerateHalfNeighborListNoWarp<<<blocks_per_grid, BLOCK_SIZE,
                                0,
                                 0>>>(
          per_atom_cell_id, in_atom_list_start_index, in_atom_list_end_index,
          cutoff_2, total_atom_num, px, py, pz, max_neighbor_num,
          neighbor_start, neighbor_end, neighbors, d_box, should_realloc,
          neighbor_cell, neighbor_cell_num, without_pbc));
}

}  // namespace op