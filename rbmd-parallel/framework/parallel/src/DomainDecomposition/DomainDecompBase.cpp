#include "DomainDecompBase.h"

#include <cmath>
#include <csignal>

#include "Forcehelper.hpp"




void DomainDecompBase::PopulateHaloLayerWithCopiesDirect(const HaloRegion& haloRegion,
                                                         LinkedCell* moleculeContainer, bool positionCheck) const {
    double shift[3];
    for (int dim = 0; dim < 3; dim++) {
        shift[dim] = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);
        shift[dim] *= haloRegion.offset[dim] * (-1);  // for halo regions the shift has to be multiplied with -1; as the
        // offset is in negative direction of the shift
    }

#if defined (_OPENMP)
#pragma omp parallel
#endif
    {
        // traverse and gather all boundary particles in the cells
        for (auto i = moleculeContainer->RegionIterator(haloRegion.rmin, haloRegion.rmax,
                                                        LinkedCellIterator::ONLY_INNER_AND_BOUNDARY);
             i.isValid(); ++i) {
            Atom m = *i;
            for (int dim = 0; dim < 3; dim++) {
                if (shift[dim] != 0) {
                    m.Setr(dim, m.r(dim) + shift[dim]);
                    if (positionCheck) {
                        // checks if the molecule has been shifted to inside the domain due to rounding errors.
                        if (shift[dim] < 0) {
                            // if the shift was negative, it is now in the lower part of the domain -> min
                            if (m.r(dim) >= moleculeContainer->getBoundingBoxMin(dim)) {
                                // in the lower part it was wrongly shifted if it is at least boxMin
                                Real_mpi r = moleculeContainer->getBoundingBoxMin(dim);
                                m.Setr(dim, std::nexttoward(
                                        r, r - 1.f));  // ensures that r is smaller than the boundingboxmin
                            }
                        } else {  // shift > 0
                            if (m.r(dim) < moleculeContainer->getBoundingBoxMax(dim)) {
                                // in the lower part it was wrongly shifted if
                                // std::nextafter: returns the next bigger value of _boundingBoxMax
                                Real_mpi r = moleculeContainer->getBoundingBoxMax(dim);
                                m.Setr(dim, r);  // ensures that r is at least boundingboxmax
                            }
                        }
                    }
                }
            }
            //GlobalLogger::info("Halo calling in PopulateHaloLayerWithCopiesDirect");
            moleculeContainer->AddHaloParticle(m);
        }
    }
}


void DomainDecompBase::PopulateHaloLayerWithCopies(unsigned int dim, LinkedCell *moleculeContainer) const {

    double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

    // molecules that have crossed the lower boundary need a positive shift
    // molecules that have crossed the higher boundary need a negative shift
    // loop over -+1 for dim=0, -+2 for dim=1, -+3 for dim=2
    const int sDim = dim+1;
    double interactionLength = moleculeContainer->GetCutoff();

    for(int direction = -sDim; direction < 2*sDim; direction += 2*sDim) {
        double shift = copysign(shiftMagnitude, static_cast<double>(-direction));

        double startRegion[3]{moleculeContainer->getBoundingBoxMin(0) - interactionLength,
                              moleculeContainer->getBoundingBoxMin(1) - interactionLength,
                              moleculeContainer->getBoundingBoxMin(2) - interactionLength};
        double endRegion[3]{moleculeContainer->getBoundingBoxMax(0) + interactionLength,
                            moleculeContainer->getBoundingBoxMax(1) + interactionLength,
                            moleculeContainer->getBoundingBoxMax(2) + interactionLength};

        if (direction < 0) {
            startRegion[dim] = moleculeContainer->getBoundingBoxMin(dim);
            endRegion[dim] = moleculeContainer->getBoundingBoxMin(dim) + interactionLength;
        } else {
            startRegion[dim] = moleculeContainer->getBoundingBoxMax(dim) - interactionLength;
            endRegion[dim] = moleculeContainer->getBoundingBoxMax(dim);
        }


#if defined (_OPENMP)
#pragma omp parallel shared(startRegion, endRegion)
#endif
        {
            // traverse and gather all boundary particles in the cells
            for (auto i = moleculeContainer->RegionIterator(startRegion, endRegion, LinkedCellIterator::ALL_CELLS);
                 i.isValid(); ++i) {
                Atom m = *i;
                m.Setr(dim, m.r(dim) + shift);
                // checks if the molecule has been shifted to inside the domain due to rounding errors.
                if (shift < 0) {  // if the shift was negative, it is now in the lower part of the domain -> min
                    if (m.r(dim) >= moleculeContainer->getBoundingBoxMin(dim)) { // in the lower part it was wrongly shifted if
                        Real_mpi r = moleculeContainer->getBoundingBoxMin(dim);
                        m.Setr(dim, std::nexttoward(r, r - 1.f));  // ensures that r is smaller than the boundingboxmin
                    }
                } else {  // shift > 0
                    if (m.r(dim) < moleculeContainer->getBoundingBoxMax(dim)) { // in the lower part it was wrongly shifted if
                        // std::nextafter: returns the next bigger value of _boundingBoxMax
                        Real_mpi r = moleculeContainer->getBoundingBoxMax(dim);
                        m.Setr(dim, std::nexttoward(r, r + 1.f));  // ensures that r is bigger than the boundingboxmax
                    }
                }
                //GlobalLogger::info("Halo calling in PopulateHaloLayerWithCopies");
                moleculeContainer->AddHaloParticle(m);
            }
        }
    }
}

void DomainDecompBase::addLeavingMolecules(std::vector<Atom> &invalidMolecules, LinkedCell *moleculeContainer) {
        for (auto& molecule : invalidMolecules) {
            for (auto dim : {0, 1, 2}) {
                auto shift = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);
                auto r = molecule.r(dim);
                if (r < moleculeContainer->getBoundingBoxMin(dim)) {
                    r = r + shift;
                    if (r >= moleculeContainer->getBoundingBoxMax(dim)) {
                        r = std::nextafter(moleculeContainer->getBoundingBoxMax(dim), -1);
                    }
                } else if (r >= moleculeContainer->getBoundingBoxMax(dim)) {
                    r = r - shift;
                    if (r < moleculeContainer->getBoundingBoxMin(dim)) {
                        r = moleculeContainer->getBoundingBoxMin(dim);
                    }
                }
                molecule.Setr(dim, r);
            }
        }
        moleculeContainer->addParticles(invalidMolecules);
        invalidMolecules.clear();

}

void DomainDecompBase::HandleDomainLeavingParticles(unsigned int dim, LinkedCell *moleculeContainer) const  {
    const double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

    // molecules that have crossed the lower boundary need a positive shift
    // molecules that have crossed the higher boundary need a negative shift
    // loop over -+1 for dim=0, -+2 for dim=1, -+3 for dim=2
    const int sDim = dim+1;
    for(int direction = -sDim; direction < 2*sDim; direction += 2*sDim) {
        double shift = copysign(shiftMagnitude, static_cast<double>(-direction));

        double cutoff = moleculeContainer->GetCutoff();
        double startRegion[3]{moleculeContainer->getBoundingBoxMin(0) - cutoff,
                              moleculeContainer->getBoundingBoxMin(1) - cutoff,
                              moleculeContainer->getBoundingBoxMin(2) - cutoff};
        double endRegion[3]{moleculeContainer->getBoundingBoxMax(0) + cutoff,
                            moleculeContainer->getBoundingBoxMax(1) + cutoff,
                            moleculeContainer->getBoundingBoxMax(2) + cutoff};

        if (direction < 0) {
            endRegion[dim] = moleculeContainer->getBoundingBoxMin(dim);
        } else {
            startRegion[dim] = moleculeContainer->getBoundingBoxMax(dim);
        }

#if defined (_OPENMP)
#pragma omp parallel shared(startRegion, endRegion)
#endif
        {
            // traverse and gather all halo particles in the cells
            for (auto i = moleculeContainer->RegionIterator(startRegion, endRegion, LinkedCellIterator::ALL_CELLS);
                 i.isValid(); ++i) {
                Atom m = *i;
                m.Setr(dim, m.r(dim) + shift);
                // some additional shifting to ensure that rounding errors do not hinder the correct placement
                if (shift < 0) {  // if the shift was negative, it is now in the lower part of the domain -> min
                    if (m.r(dim) <= moleculeContainer->getBoundingBoxMin(dim)) { // in the lower part it was wrongly shifted if
                        m.Setr(dim, moleculeContainer->getBoundingBoxMin(dim));  // ensures that r is at least the boundingboxmin
                    }
                } else {  // shift > 0
                    if (m.r(dim) >= moleculeContainer->getBoundingBoxMax(dim)) { // in the lower part it was wrongly shifted if
                        // std::nexttoward: returns the next bigger value of _boundingBoxMax
                        Real_mpi r = moleculeContainer->getBoundingBoxMax(dim);
                        m.Setr(dim, std::nexttoward(r, r - 1.f));  // ensures that r is smaller than the boundingboxmax
                    }
                }
                moleculeContainer->AddParticle(m);
                moleculeContainer->deleteMolecule(i, false); // removeFromContainer = true;
            }
        }
    }
}

void DomainDecompBase::handleForceExchange(unsigned int dim, LinkedCell *moleculeContainer) const  {
    const double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

    // direction +1/-1 for dim 0, +2/-2 for dim 1, +3/-3 for dim 2
    //const int direction = dim+1;
    const int sDim = dim + 1;
    for (int direction = -sDim; direction < 2 * sDim; direction += 2 * sDim) {
        double shift = copysign(shiftMagnitude, static_cast<double>(-direction));

        // Loop over all halo particles in the positive direction
        double cutoff = moleculeContainer->GetCutoff();
        double startRegion[3]{moleculeContainer->getBoundingBoxMin(0) - cutoff,
                              moleculeContainer->getBoundingBoxMin(1) - cutoff,
                              moleculeContainer->getBoundingBoxMin(2) - cutoff};
        double endRegion[3]{moleculeContainer->getBoundingBoxMax(0) + cutoff,
                            moleculeContainer->getBoundingBoxMax(1) + cutoff,
                            moleculeContainer->getBoundingBoxMax(2) + cutoff};

        if (direction < 0) {
            endRegion[dim] = moleculeContainer->getBoundingBoxMin(dim);
        } else {
            startRegion[dim] = moleculeContainer->getBoundingBoxMax(dim);
        }

#if defined (_OPENMP)
#pragma omp parallel shared(startRegion, endRegion)
#endif
        {
            double shiftedPosition[3];

            decltype(moleculeContainer->GetMoleculeAtPosition(shiftedPosition)) originalPreviousIter{};

            for (auto haloIter = moleculeContainer->RegionIterator(startRegion, endRegion, LinkedCellIterator::ALL_CELLS);
                 haloIter.isValid(); ++haloIter) {
                // Add force of halo particle to original particle (or other duplicates)
                // that have a distance of -'shiftMagnitude' in the current direction
                shiftedPosition[0] = haloIter->r(0);
                shiftedPosition[1] = haloIter->r(1);
                shiftedPosition[2] = haloIter->r(2);
                shiftedPosition[dim] += shift;
                originalPreviousIter = Forcehelper::AddValuesAndGetIterator(moleculeContainer, shiftedPosition, originalPreviousIter, *haloIter);
            }
        }
    }
}