#pragma once


#include <vector>
#include "LogTool.h"
#include "Cell.h"
#include "SingleCellIterator.hpp"
#include <ciso646>

class LinkedCellIterator {   // TODO 到CellManager
public:
    enum Type {
        ALL_CELLS=0, /* iterates every cell */
        ONLY_INNER_AND_BOUNDARY=1, /* iterates only inner and boundary cells, i.e. no halo cells */
    };

    typedef std::vector<Cell> CellContainer_T;
    typedef CellContainer_T* CellContainer_T_ptr;
    typedef size_t CellIndex_T;
    typedef size_t MolIndex_T;

    LinkedCellIterator ();
    LinkedCellIterator (Type t_arg, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize=true);
    LinkedCellIterator& operator=(const LinkedCellIterator& other);
    virtual ~LinkedCellIterator() = default;
    constexpr LinkedCellIterator(const LinkedCellIterator&) = default;

    Atom& operator *  () const;
    Atom* operator -> () const;

    void deleteCurrentParticle();

    CellIndex_T getCellIndex(){return _cell_index;}

    bool isValid() const {
        return _cells != nullptr and _cell_index < _cells->size() and _cell_iterator.isValid();
    }

    virtual void operator ++ ();
protected:
    SingleCellIterator<Cell> _cell_iterator;
    Type _type;
    CellContainer_T_ptr _cells;

    CellIndex_T _cell_index;

    const CellIndex_T _stride;

    virtual void next_non_empty_cell();
    virtual void updateCellIteratorCell();
};


inline LinkedCellIterator :: LinkedCellIterator () : _cell_iterator(), _type(ALL_CELLS), _cells (nullptr), _cell_index (0), _stride (1) {
}

inline LinkedCellIterator :: LinkedCellIterator (Type t_arg, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize) :
        _cell_iterator(&(cells_arg->front())), _type(t_arg), _cells (cells_arg), _cell_index (offset_arg), _stride (stride_arg) {
    if ( !initialize ) {
        return;
    }

    updateCellIteratorCell();

//    assert(_cells != nullptr);
    if (_cells == nullptr){
        GlobalLogger::error("cells is empty");
    }

    const CellContainer_T& cells = *_cells;

    if(_cell_index < cells.size()) {
        if(cells.at(_cell_index).IsEmpty() or (_type == ONLY_INNER_AND_BOUNDARY and cells.at(_cell_index).IsHaloCell())) {
            next_non_empty_cell();
        }
        /*
        else {
            // leave object as is
        }
        */
    }
}

inline LinkedCellIterator& LinkedCellIterator::operator=(const LinkedCellIterator& other) {
//    assert(_stride == other._stride);
    if ((_stride != other._stride)){
        GlobalLogger::error("_stride != other._stride");
    }
    _cell_iterator = other._cell_iterator;
    _type = other._type;
    _cells = other._cells;
    _cell_index = other._cell_index;
    return *this;
}

inline void LinkedCellIterator :: next_non_empty_cell() {
    //assert(_cells != nullptr);
    if (_cells == nullptr){
        GlobalLogger::error("cells is empty");
    }
    const CellContainer_T& cells = *_cells;
    const CellIndex_T numCells = cells.size();

    // find the next non-empty cell

    for (_cell_index += _stride; _cell_index < numCells; _cell_index += _stride) {
        Cell c = cells.at(_cell_index);

        // if we want only inner/boundary cells: check if it is a halo cell
        if(_type == ONLY_INNER_AND_BOUNDARY and c.IsHaloCell()){
            continue;
        }

        // only use this cell if it is not empty
        if(c.IsNotEmpty()) {
            updateCellIteratorCell();
            break;
        }
    }
}

inline void LinkedCellIterator :: operator ++ () {

    if (_cell_iterator.isValid()) {
        ++_cell_iterator;
    }

    // don't merge into if-else, _cell_iterator may become invalid after ++

    if (not _cell_iterator.isValid()) {
        next_non_empty_cell();
    }
}

inline Atom& LinkedCellIterator :: operator * () const {
    // .at method performs automatically an out-of-bounds check
//    assert(&_cells->at(_cell_index) == dynamic_cast<const Cell*>(_cell_iterator.getCell()));
    if (&_cells->at(_cell_index) != dynamic_cast<const Cell*>(_cell_iterator.getCell()))
    {
        GlobalLogger::error("LinkedCellIterator :: operator * () ERROR");
    }
    return _cell_iterator.operator *();
}

// no clue why this returns a pointer
inline Atom* LinkedCellIterator:: operator -> () const {
    return &(this->operator*());
}

inline void LinkedCellIterator :: deleteCurrentParticle () {
    _cell_iterator.deleteCurrentParticle();
}

inline void LinkedCellIterator :: updateCellIteratorCell() {
    if(_cell_index < _cells->size()) {
        _cell_iterator = SingleCellIterator<Cell>(&_cells->at(_cell_index));
    }
}




