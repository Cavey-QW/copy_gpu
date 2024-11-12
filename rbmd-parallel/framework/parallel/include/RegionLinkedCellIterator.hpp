#pragma once
#include "LinkedCellIterator.hpp"
#include <ciso646>

class RegionLinkedCellIterator : public LinkedCellIterator{
public:
    RegionLinkedCellIterator ();
    RegionLinkedCellIterator (Type t, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const int startCellIndex_arg, const int regionDimensions_arg[3], const int globalDimensions_arg[3], const double startRegion_arg[3], const double endRegion_arg[3]);
    RegionLinkedCellIterator& operator=(const RegionLinkedCellIterator& other);

    ~RegionLinkedCellIterator()  override = default;

    constexpr RegionLinkedCellIterator(const RegionLinkedCellIterator&) = default;

    void operator++() override;
private:

    CellIndex_T getGlobalCellIndex();
    void next_non_empty_cell() override;

    CellIndex_T _baseX;
    CellIndex_T _baseY;
    CellIndex_T _baseZ;
    CellIndex_T _localCellIndex;
    CellIndex_T _regionDimensions[3];
    CellIndex_T _globalDimensions[3];

    double _startRegion[3];
    double _endRegion[3];
};


inline RegionLinkedCellIterator :: RegionLinkedCellIterator () : LinkedCellIterator(), _localCellIndex(0) {
}

inline RegionLinkedCellIterator :: RegionLinkedCellIterator (Type t, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const int startCellIndex_arg, const int regionDimensions_arg[3], const int globalDimensions_arg[3], const double startRegion_arg[3], const double endRegion_arg[3]) :
        LinkedCellIterator(t, cells_arg, offset_arg, stride_arg, false), _localCellIndex (offset_arg) {
    
    for (int d = 0; d < 3; d++) {
        _regionDimensions[d] = regionDimensions_arg[d];
        _globalDimensions[d] = globalDimensions_arg[d];
        _startRegion[d] = startRegion_arg[d];
        _endRegion[d] = endRegion_arg[d];
    }

    _baseZ = startCellIndex_arg / (_globalDimensions[0] * _globalDimensions[1]);
    _baseY = (startCellIndex_arg % (_globalDimensions[0] * _globalDimensions[1])) / _globalDimensions[0];
    _baseX = (startCellIndex_arg % (_globalDimensions[0] * _globalDimensions[1])) % _globalDimensions[0];

    _cell_index = getGlobalCellIndex();
    updateCellIteratorCell();

    rbmd_assert(_cells != nullptr);

    const CellContainer_T& cells = *_cells;
    const CellIndex_T numCellsInRegion = _regionDimensions[2] * _regionDimensions[1] * _regionDimensions[0];

    // if _cell_index is out of the region => invalid <=> (_localCellIndex >= numCellsInRegion)
    if (_localCellIndex < numCellsInRegion) {
        if (cells[_cell_index].IsNotEmpty() and (this->operator*()).InBox(_startRegion, _endRegion)) {
            // if the particle is good, leave it
        } else {
            // else, find a particle in the box
            this->operator++();
        }
    } else {
        // set to invalid
        _cells = nullptr;
    }
}

inline RegionLinkedCellIterator& RegionLinkedCellIterator::operator=(const RegionLinkedCellIterator& other) {
    rbmd_assert(_stride == other._stride);
    _cells = other._cells;
    _cell_index = other._cell_index;
    _cell_iterator = other._cell_iterator;
    _baseX = other._baseX;
    _baseY = other._baseY;
    _baseZ = other._baseZ;
    _localCellIndex = other._localCellIndex;
    for(int d = 0; d < 3; d++){
        _regionDimensions[d] = other._regionDimensions[d];
        _globalDimensions[d] = other._globalDimensions[d];
        _startRegion[d] = other._startRegion[d];
        _endRegion[d] = other._endRegion[d];
    }
    return *this;
}

inline void RegionLinkedCellIterator :: operator ++() {
    do{
        LinkedCellIterator :: operator++();
    } while (isValid() and !(this->operator*()).InBox(_startRegion, _endRegion));
}

inline void RegionLinkedCellIterator :: next_non_empty_cell() {
    //cellIndex should always be the index in the cell array (_cells member variable in LinkedCells)
    rbmd_assert(_cells != nullptr);

    const CellContainer_T& cells = *_cells;
    const CellIndex_T numCellsInRegion = _regionDimensions[2] * _regionDimensions[1] * _regionDimensions[0];

    // find the next non-empty cell
    for (_localCellIndex += _stride; _localCellIndex < numCellsInRegion; _localCellIndex += _stride) {

        Cell c = cells.at(getGlobalCellIndex());

        if (c.IsNotEmpty() and (_type == ALL_CELLS or not c.IsHaloCell())) {
            _cell_index = getGlobalCellIndex();
            updateCellIteratorCell();
            break;
        }
    }
}

inline LinkedCellIterator::CellIndex_T RegionLinkedCellIterator :: getGlobalCellIndex() {
    CellIndex_T dz = _localCellIndex / (_regionDimensions[0] * _regionDimensions[1]);
    CellIndex_T dy = (_localCellIndex % (_regionDimensions[0] * _regionDimensions[1])) / _regionDimensions[0];
    CellIndex_T dx = (_localCellIndex % (_regionDimensions[0] * _regionDimensions[1])) % _regionDimensions[0];
    return (_baseX + dx) + (_baseY + dy) * _globalDimensions[0] + (_baseZ + dz) * _globalDimensions[0] * _globalDimensions[1];
}