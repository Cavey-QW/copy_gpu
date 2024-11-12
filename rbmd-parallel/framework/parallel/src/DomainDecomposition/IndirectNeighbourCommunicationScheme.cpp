#include <csignal>
#include "IndirectNeighbourCommunicationScheme.h"

#include "DomainServiceLocator.h"
#include "NeighborAcquirer.h"

void IndirectNeighbourCommunicationScheme::initExchangeMoleculesMPI1D(LinkedCell *moleculeContainer,
MessageType msgType, bool removeRecvDuplicates, unsigned short d, DomainDecompBase *domainDecomp) {
    if (_coversWholeDomain[d]) {
        // use the sequential version

        switch (msgType) {
            case LEAVING_AND_HALO_COPIES:
                domainDecomp->HandleDomainLeavingParticles(d, moleculeContainer);           // TODO 放到基类是不是有问题呢？
            domainDecomp->PopulateHaloLayerWithCopies(d, moleculeContainer);
            break;
            case LEAVING_ONLY:
                domainDecomp->HandleDomainLeavingParticles(d, moleculeContainer);
            break;
            case HALO_COPIES:
                domainDecomp->PopulateHaloLayerWithCopies(d, moleculeContainer);
            break;
            case FORCES:
                domainDecomp->handleForceExchange(d, moleculeContainer);
            break;
        }

    } else {

        const int numNeighbours = (*_neighbours)[d].size();
        std::vector<Atom> dummy;
        for (int i = 0; i < numNeighbours; ++i) {
            //Log::global_log->debug() << "Rank " << domainDecomp->getRank() << " is initiating communication to" << std::endl;
            (*_neighbours)[d][i].initSend(moleculeContainer, domainDecomp->GetMPIComm(),
                    domainDecomp->getMPIParticleType(), msgType, dummy, false, true/*do halo position change*/);
        }

    }
}

void IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1D(LinkedCell *moleculeContainer,
    MessageType msgType, bool removeRecvDuplicates, unsigned short d, DomainDecompBase *domainDecomp) {
	if (_coversWholeDomain[d]) {
		return;
	}
	const int numNeighbours = (*_neighbours)[d].size();
	// the following implements a non-blocking recv scheme, which overlaps unpacking of
	// messages with waiting for other messages to arrive
	bool allDone = false;
	double startTime = MPI_Wtime();

	double waitCounter = 50.0;
	double deadlockTimeOut = 360.0;
	// Log::global_log->set_mpi_output_all();
	for (int i = 0; i < numNeighbours; ++i) { // reset receive status
		if (domainDecomp->GetCurrentRank() != (*_neighbours)[d][i].getRank()) {
			(*_neighbours)[d][i].resetReceive();
		}
	}

	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= (*_neighbours)[d][i].testSend();
		}

		// get the counts and issue the Irecv-s
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= (*_neighbours)[d][i].iprobeCount(domainDecomp->GetMPIComm(),
					domainDecomp->getMPIParticleType());
		}

		// testRecv
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= (*_neighbours)[d][i].testRecv(moleculeContainer, removeRecvDuplicates, msgType==FORCES);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			// Log::global_log->warning()
			// 		<< "IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
			// 		<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds"
			// 		<< std::endl;
			waitCounter += 5.0;
			for (int i = 0; i < numNeighbours; ++i) {
				(*_neighbours)[d][i].deadlockDiagnosticSendRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			// Log::global_log->error()
			// 		<< "IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock error: Rank "
			// 		<< domainDecomp->getRank() << " is waiting for more than " << deadlockTimeOut << " seconds"
			// 		<< std::endl;
			for (int i = 0; i < numNeighbours; ++i) {
				(*_neighbours)[d][i].deadlockDiagnosticSendRecv();
			}
			//Simulation::exit(457);
		}

	} // while not allDone
	//Log::global_log->set_mpi_output_root(0);
}

void IndirectNeighbourCommunicationScheme::exchangeMoleculesMPI1D(LinkedCell *moleculeContainer, MessageType msgType,
bool removeRecvDuplicates, unsigned short d, DomainDecompBase *domainDecomp) {
//    MPI_Barrier(MPI_COMM_WORLD); // wait for all processes to reach this point (to avoid race conditions in output)
//    sleep(15); // wait for all processes to reach this point (to avoid race conditions in output)
	initExchangeMoleculesMPI1D(moleculeContainer,  msgType, removeRecvDuplicates, d, domainDecomp);
	finalizeExchangeMoleculesMPI1D(moleculeContainer,  msgType, removeRecvDuplicates, d, domainDecomp);

}

void IndirectNeighbourCommunicationScheme::convert1StageTo3StageNeighbours(
	const std::vector<CommunicationPartner> &commPartners, std::vector<std::vector<CommunicationPartner>> &neighbours,
	HaloRegion &ownRegion, double cutoffRadius) {
	//TODO: extend for anything else than full shell
	//TODO: implement conversion of 1StageTo3StageNeighbours

	for (const CommunicationPartner& commPartner : commPartners) {
		if (!commPartner.isFaceCommunicator()) {
			continue;  // if commPartner is not a face sharing communicator, we can ignore it!
		}
		unsigned int d = commPartner.getFaceCommunicationDirection();
		neighbours[d].push_back(commPartner);
		neighbours[d].back().enlargeInOtherDirections(d, cutoffRadius); // do this more wisely if multiple neighbours exist in that direction.
	}
}

void IndirectNeighbourCommunicationScheme::exchangeMoleculesMPI(LinkedCell *moleculeContainer, MessageType msgType,
                                                                bool removeRecvDuplicates, DomainDecompBase *domainDecomp, bool doHaloPositionCheck) {
    for (unsigned int d = 0; d < GetCommDims(); d++) {
        exchangeMoleculesMPI1D(moleculeContainer, msgType, removeRecvDuplicates, d, domainDecomp);
    }
}

void IndirectNeighbourCommunicationScheme::initCommunicationPartners(double cutoffRadius,
DomainDecompBase *domainDecomp, LinkedCell *moleculeContainer) { // if this one is used, push pull should not (at least for now) be set

    // corners of the process-specific domain
    double rmin[DIM]; // lower corner
    double rmax[DIM]; // higher corner

    for (int d = 0; d < DIM; d++) {
        rmin[d] = domainDecomp->GetBoundingBoxMin(d);
        rmax[d] = domainDecomp->GetBoundingBoxMax(d);

        // TODO: this should be safe, as long as molecules don't start flying around
        // at the speed of one cutoffRadius per time step
    }

    for (unsigned int d = 0; d < _commDimms; d++) {
        (*_neighbours)[d].clear();
    }
    HaloRegion ownRegion = { rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0, cutoffRadius}; // region of the box
    // ---
    std::vector<HaloRegion> haloRegions = _zonalMethod->getLeavingExportRegions(ownRegion, cutoffRadius, _coversWholeDomain); // halo regions (outside of the box)
    // ---
    std::vector<CommunicationPartner> commPartners;
    for (HaloRegion haloRegion : haloRegions) { // determine who to communicate with - who's region is in the haloRegions vector?
        auto newCommPartners = domainDecomp->GetNeighboursFromHaloRegion(haloRegion, cutoffRadius);
        commPartners.insert(commPartners.end(), newCommPartners.begin(), newCommPartners.end());
    }
    _fullShellNeighbours = commPartners;   // 2024年4月2日 这里确定没有问题
    convert1StageTo3StageNeighbours(commPartners, (*_neighbours), ownRegion, cutoffRadius);
    //squeeze neighbours -> only a single send, if rightneighbour == leftneighbour
    for (unsigned int d = 0; d < _commDimms; d++) {
        (*_neighbours)[d]= NeighborAcquirer::squeezePartners((*_neighbours)[d]);
    }
}

void IndirectNeighbourCommunicationScheme::prepareNonBlockingStageImpl(LinkedCell *moleculeContainer,
unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates, DomainDecompBase *domainDecomp) {
    rbmd_assert(stageNumber < GetCommDims());
	initExchangeMoleculesMPI1D(moleculeContainer, msgType, removeRecvDuplicates, stageNumber, domainDecomp);
}

void IndirectNeighbourCommunicationScheme::finishNonBlockingStageImpl(LinkedCell *moleculeContainer,
unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates, DomainDecompBase *domainDecomp) {
    rbmd_assert(stageNumber < GetCommDims());
	finalizeExchangeMoleculesMPI1D(moleculeContainer, msgType, removeRecvDuplicates, stageNumber, domainDecomp);
}

