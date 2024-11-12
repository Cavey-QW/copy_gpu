#include "Cell.h"
#include "CellManager.h"
#include "SingleCellIterator.hpp"
#include "LogTool.h"
#include <ciso646>
Cell::Cell() = default;

Cell::~Cell() = default;

unsigned long Cell::GetID() const { return _id; }

void Cell::SetID(unsigned long id) { _id = id;}   //2024年4月1日 这里的问题。要检查还有没有没有实现的函数

void Cell::EmptyAtoms() {
    _atoms.clear();
}

bool Cell::AddAtomToCell(Atom &atom) {
    bool inserted = false;
    _atoms.emplace_back(atom);
    inserted = true;
    return inserted;
}

unsigned long Cell::GetAtomsCount() const {
    return _atoms.size();
}

bool Cell::IsEmpty() const {
    return _atoms.empty();
}

double Cell::GetBoundingBoxMin(int dimension) const {
    return CellManager::Instance().GetCellBoundingBoxMin(_id, dimension);
}

double Cell::GetBoundingBoxMax(int dimension) const {
    return CellManager::Instance().GetCellBoundingBoxMax(_id, dimension);
}

bool Cell::IsHaloCell() const {
    return CellManager::Instance().IsHaloCell(_id);
}

bool Cell::IsBoundaryCell() const {
    return CellManager::Instance().IsBoundaryCell(_id);
}

bool Cell::IsInnerCell() const {
    return CellManager::Instance().IsInnerCell(_id);
}

bool Cell::IsInnerMostCell() const {
    return CellManager::Instance().IsInnerMostCell(_id);
}

void Cell::AssignCellToHaloRegion() {
    //assert(IsHaloCell());
    if (!IsHaloCell()){
        GlobalLogger::error("!IsHaloCell()!!");
    }
}

bool Cell::InBox(const Atom &atom) const {
    double bounding_box_min[3] = {GetBoundingBoxMin(0),GetBoundingBoxMin(1),GetBoundingBoxMin(2)};
    double bounding_box_max[3] = {GetBoundingBoxMax(0),GetBoundingBoxMax(1),GetBoundingBoxMax(2)};
    return atom.InBox(bounding_box_min,bounding_box_max);
}

SingleCellIterator<Cell> Cell::iterator() {
    return SingleCellIterator<Cell>(this);
}

bool Cell::IsNotEmpty() const {
    return not IsEmpty();
}

template<typename T, typename A, typename UIntType>
void fastRemove(std::vector<T, A>& v, UIntType index) {
    using std::swap;

    rbmd_assert(std::is_integral<UIntType>::value);

    // assumption: if the vector is empty, then the method is not called at all
//    // i.e. v is not empty()
//    assert(not v.empty());
//
//    // assert that the iterator points inside this vector
//    assert(index >= 0);
//    assert(index < v.size());
    if (v.empty() || index<0 || index> v.size()){
        GlobalLogger::error("error in fastRemove");
    }
    // either std::swap, or a user-written function in the namespace
    swap(v[index], v.back());
    // note: swapping is necessary (as opposed to only copying )

    v.pop_back(); // end() will change and may become equal to pos, which is intended
}

bool Cell::deleteMoleculeByIndex(size_t index) {
    bool found = true;
    fastRemove(_atoms, index);
    return found;}

bool Cell::TestPointInCell(const double *point) const {
    double boxMin[3] =  {GetBoundingBoxMin(0),GetBoundingBoxMin(1),GetBoundingBoxMin(2)};
    double boxMax[3] = {GetBoundingBoxMax(0),GetBoundingBoxMax(1),GetBoundingBoxMax(2)};
    return boxMin[0] <= point[0] && boxMin[1] <= point[1] && boxMin[2] <= point[2] &&
           point[0] < boxMax[0] && point[1] < boxMax[1] && point[2] < boxMax[2];
}

void Cell::increaseMoleculeStorage(size_t numExtraMols) {
    _atoms.reserve(_atoms.size() + numExtraMols);
}

void Cell::PreUpdateLeavingMolecules() {
    _leaveing_atoms.clear();

#ifndef NDEBUG
    const size_t size_total = _atoms.size(); // for debugging, see below
#endif

    for (auto it = iterator(); it.isValid(); ++it) {
        //it->setSoA(nullptr);

        const bool isStaying = InBox(*it);

        if (isStaying) {
            // don't do anything, just advance iterator
        } else {
            _leaveing_atoms.push_back(*it);
            it.deleteCurrentParticle();
        }
    }

//    assert(_atoms.size() + _leaveing_atoms.size() == size_total); // any molecules lost?
    if ((_atoms.size() + _leaveing_atoms.size() != size_total)){
        GlobalLogger::error("atom lost in PreUpdateLeavingMolecules");
    }
}


void Cell::UpdateLeavingMolecules(Cell &otherCell) {
    for (auto m = otherCell._leaveing_atoms.begin(); m != otherCell._leaveing_atoms.end(); ++m) { // loop over all indices
        if (InBox(*m)) { // if molecule moves in this cell
            AddAtomToCell(*m);
        }
    }
}

void Cell::PostUpdateLeavingMolecules() {
    _leaveing_atoms.clear();
}

std::vector<Atom> &Cell::GetAtoms() {
    return this->_atoms;
}

