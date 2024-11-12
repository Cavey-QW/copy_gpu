#pragma once
#include "Atom.h"


typedef enum {
    MOLECULE_MOLECULE = 0,      /**< molecule molecule */
    MOLECULE_HALOMOLECULE = 1,  /**< molecule - halo molecule */
    MOLECULE_MOLECULE_FLUID = 2 /**< molecule - molecule (fluid) */
} PairType;


class ParticlePairsHandler {
public:
	//! Constructor
	ParticlePairsHandler() {}

	//! Destructor
	virtual ~ParticlePairsHandler() {
	}

	//! @brief things to be done before particle pairs are processed
	virtual void init() = 0;

	//! @brief things to be done after particle pairs are processed
	virtual void finish() = 0;

	//! @brief things to be done for each particle pair
	//!
	//! @param particle1 first particle
	//! @param particle2 second particle
	//! @param distanceVector[3] distance between the two particles
	//! @param pairType describes whether the pair is a original pair(0) or a duplicated pair(1)
	//!                 for details about pair types see comments on traversePairs() in ParticleContainer
	virtual double processPair(Atom& particle1, Atom& particle2, double distanceVector[3], PairType pairType, double dd, bool calculateLJ) = 0;

protected:
	// RDF* _rdf; calculation moved to RDFCellProcessor! Legacy implementation will be slower, but we want to always use the vectorized one, which is now faster.
};






