#pragma once
#include <cstddef>
#include "Cell.h"

template<class Cell>
class SingleCellIterator {
public:
    SingleCellIterator();
    SingleCellIterator(Cell *cell_arg,
                       size_t index_arg = 0) : _cell(cell_arg), _mol_index(index_arg), _currentParticleDeleted(false) {
    }
    SingleCellIterator& operator=(const SingleCellIterator& other);
    ~SingleCellIterator(){}

    constexpr SingleCellIterator(const SingleCellIterator&) = default;

    Atom& operator *  () const {
        // .at method performs automatically an out-of-bounds check
        Atom *moleculePtr = nullptr;


        _cell->moleculesAtNew(_mol_index, moleculePtr);

        return *moleculePtr;
    }
    Atom* operator -> () const;

    void deleteCurrentParticle() {
        _cell->deleteMoleculeByIndex(_mol_index);
        _currentParticleDeleted = true;
    }

    size_t getIndex() const {
        return _mol_index;
    }

    Cell * getCell() const { return _cell; }

    bool isValid() const {
        return _cell != nullptr and _mol_index < static_cast<size_t>(_cell->GetAtomsCount());
    }

    void operator ++() {
        // if the "current" particle was deleted, then there is a new particle at _mol_index
        // and the _mol_index value should not be incremented.
        _mol_index += _currentParticleDeleted ? 0 : 1;

        _currentParticleDeleted = false;
    }

private:

    Cell * _cell;
    size_t _mol_index;
    bool _currentParticleDeleted;

};

template<class Cell>
inline SingleCellIterator<Cell>::SingleCellIterator() : _cell(nullptr), _mol_index(0), _currentParticleDeleted(false) {
}

template<class Cell>
inline SingleCellIterator<Cell>& SingleCellIterator<Cell>::operator=(const SingleCellIterator& other) {
    _cell = other._cell;
    _mol_index = other._mol_index;
    _currentParticleDeleted = other._currentParticleDeleted;
    // no operator= for Molecule and it should not be needed.
    return *this;
}

// no clue why this returns a pointer
template<class Cell>
inline Atom* SingleCellIterator<Cell>:: operator -> () const {
    return &(this->operator*());
}

