#include "CellManager.h"
#include "LogTool.h"
#include "IdexMappingTool.hpp"
#include <ciso646>
CellManager &CellManager::Instance() {
    static CellManager instance;
    return instance;
}

std::array<unsigned long, 3> CellManager::MapLinearTo3DIndex(
        unsigned long cell_index) const {
//    assert(_initlized);
    if (!_initlized){
        GlobalLogger::error("Error in CellManager,Not Init!");
    }
    return OneToThreeD(cell_index, _num_cells_per_dimension);
}



// TODO 需要注意
bool CellManager::IsCellOfType(CellManager::CellType cell_type,
                               unsigned long cell_index) const {
    if (!_initlized){
        GlobalLogger::error("Error in CellManager,Not Init!");
    }
    bool result = false;
    auto index_3d = MapLinearTo3DIndex(cell_index);
    switch (cell_type) {
        case CellType::HALO:  // 这里真的会被调用到吗?根据halo_cell_falg已经返回了
            result = false;
            for (int dimension = 0; dimension < 3; ++dimension) {
                result |= (index_3d[dimension] == 0) or
                          (index_3d[dimension] == _num_cells_per_dimension[dimension] - 1);  // 1D
            }
            break;
        case CellType::BOUNDARY:
            result = false;
            for (int dimension = 0; dimension < 3; ++dimension) {
                result |= (index_3d[dimension] == 1) or
                          (index_3d[dimension] == _num_cells_per_dimension[dimension] - 2);  // 1D
            }
            break;
        case CellType::INNER:
            result = true;
            for (int dimension = 0; dimension < 3; ++dimension) {
                result &= (index_3d[dimension] >= 2) and
                          (index_3d[dimension] <= _num_cells_per_dimension[dimension] - 3);
            }
            break;
        case CellType::INNERMOST:
            result = true;
            for (int dimension = 0; dimension < 3; ++dimension) {
                result &= (index_3d[dimension] >= 3) and
                          (index_3d[dimension] <= _num_cells_per_dimension[dimension] - 4);
            }
            break;
    }
    return result;
}

bool CellManager::IsHaloCell(unsigned long cell_index) const {
    return _halo_cell_flags[cell_index];
}

bool CellManager::IsBoundaryCell(unsigned long cell_index) const {
    return IsCellOfType(CellType::BOUNDARY, cell_index);
}

bool CellManager::IsInnerCell(unsigned long cell_index) const {
    return IsCellOfType(CellType::INNER, cell_index);
}

bool CellManager::IsInnerMostCell(unsigned long cell_index) const {
    return IsCellOfType(CellType::INNERMOST, cell_index);
}

double CellManager::GetCellBorder(unsigned long cell_index, int dimension) const {
    if (!_initlized){
        GlobalLogger::error("Error in CellManager,Not Init!");
    }
//    assert(cell_index <= _num_cells_per_dimension[dimension]);
    if (cell_index > _num_cells_per_dimension[dimension]){
        GlobalLogger::error("cell_index > _num_cells_per_dimension[dimension]");
    }
    double result{};
    if (cell_index == 0) {
        result = _halo_bounding_box_min[dimension];
    } else if (cell_index == _halo_width_in_num_cells[dimension]) {
        result = _bounding_box_min[dimension];
    } else if (cell_index == _num_cells_per_dimension[dimension] - _halo_width_in_num_cells[dimension]) {
        result = _bounding_box_max[dimension];
    } else if (cell_index == _num_cells_per_dimension[dimension]) {
        result = _halo_bounding_box_max[dimension];
    } else {
        result = cell_index * _cell_length[dimension] + _halo_bounding_box_min[dimension];  // TODO ! narrow cell的box
    }
    return result;
}

double CellManager::GetCellBoundingBoxMin(unsigned long cell_index, int dimension) const {
    auto dimension_index = MapLinearTo3DIndex(cell_index)[dimension];
    return GetCellBorder(dimension_index, dimension);
}

double CellManager::GetCellBoundingBoxMax(unsigned long cell_index, int dimension) const {
    auto dimension_index = MapLinearTo3DIndex(cell_index)[dimension];
    return GetCellBorder(dimension_index + 1, dimension);
}

void CellManager::init(int cellsPerDim[3], double haloBoxMin[3], double haloBoxMax[3], double boxMin[3],
double boxMax[3], double cellLength[3], int haloWidthInNumCells[3]) {
    int totalNumCells = 1;
    for (int d = 0; d < 3; ++d) {
        _num_cells_per_dimension[d] = cellsPerDim[d];
        _halo_bounding_box_min[d] = haloBoxMin[d];
        _halo_bounding_box_max[d] = haloBoxMax[d];
        _bounding_box_min[d] = boxMin[d];
        _bounding_box_max[d] = boxMax[d];
        _cell_length[d] = cellLength[d];
        _halo_width_in_num_cells[d] = haloWidthInNumCells[d];

        totalNumCells *= cellsPerDim[d];
    }

    _halo_cell_flags.resize(totalNumCells);

    int runningIndex = 0;

    for (unsigned z = 0; z < _num_cells_per_dimension[2]; ++z) {
        bool isHaloZ = (z < _halo_width_in_num_cells[2] or z >= _num_cells_per_dimension[2] - _halo_width_in_num_cells[2]);

        for (unsigned y = 0; y < _num_cells_per_dimension[1]; ++y) {
            bool isHaloY = (y < _halo_width_in_num_cells[1] or y >= _num_cells_per_dimension[1] - _halo_width_in_num_cells[1]);

            for (unsigned x = 0; x < _num_cells_per_dimension[0]; ++x) {
                bool isHaloX =
                    (x < _halo_width_in_num_cells[0] or x >= _num_cells_per_dimension[0] - _halo_width_in_num_cells[0]);

                if (isHaloZ or isHaloY or isHaloX) {
                    _halo_cell_flags.at(runningIndex) = true;
                } else {
                    _halo_cell_flags.at(runningIndex) = false;
                }
                runningIndex++;
            }
        }
    }
    _initlized = true;
}

