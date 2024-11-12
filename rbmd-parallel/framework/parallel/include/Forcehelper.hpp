
#include <variant>
#include <ctgmath>
#include <iostream>
#include "LinkedCell.h"
#include "LinkedCellIterator.hpp"

class Forcehelper{
public:
    static std::variant<LinkedCellIterator, SingleCellIterator<Cell>> AddValuesAndGetIterator(
            LinkedCell* moleculeContainer, const double* position,
            std::variant<LinkedCellIterator, SingleCellIterator<Cell>> previousIterator,
            Atom& haloMolecule) {
// Reuse the previous iterator in case the particles match.
        bool usePreviousIterator = false;
        std::visit(
                [&](auto& previousIter) {
                    if (previousIter.isValid()) {
// Increment pointer!
                        ++previousIter;

                        if (previousIter.isValid()) {
                            const double epsi = moleculeContainer->GetCutoff() * 1e-6;
// If the particle pointed to by previousIter is close to the given position, use it.
                            if (fabs(previousIter->r(0) - position[0]) <= epsi and
                                fabs(previousIter->r(1) - position[1]) <= epsi and
                                fabs(previousIter->r(2) - position[2]) <= epsi) {
                                usePreviousIterator = true;
                            }
                        }
                    }
                },
                previousIterator);

// If the previousIterator matches, use it!
        auto originalIter = usePreviousIterator ? previousIterator : moleculeContainer->GetMoleculeAtPosition(position);

        std::visit(
                [&](auto originalIter) {
                    if (not originalIter.isValid()) {
// This should not happen
                        std::cout << "Original molecule not usePreviousIterator";
                        exit(1);
                    }

                    rbmd_assert(originalIter->GetID() == haloMolecule.GetID());

// Add the force values.
                    originalIter->Fadd(haloMolecule.F_arr().data());
                    originalIter->Madd(haloMolecule.M_arr().data());
                    originalIter->Viadd(haloMolecule.Vi_arr().data());
                },
                originalIter);
        return originalIter;
    }


};


