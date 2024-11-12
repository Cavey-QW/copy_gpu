#pragma  once
#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include "CellProcessor.h"
#include "OriginalCellPairTraversal.hpp"


template<class CellTemplate>
class TraversalTuner {
	friend class LinkedCellsTest;

public:
	// Probably remove this once autotuning is implemented
	enum traversalNames {
		ORIGINAL = 0,
		// C08      = 1,
		// C04      = 2,
		// SLICED   = 3,
		// HS       = 4,
		// MP       = 5,
		// C08ES    = 6,
		// NT       = 7,
		// // quicksched has to be the last traversal!
		// QSCHED   = 8,
	};

	TraversalTuner();

	~TraversalTuner();

	void findOptimalTraversal();


	/**
	 * Rebuild the traversals.
	 * @param cells The vector of cells.
	 * @param dims The dimensions (cells per dimension, including halo!)
	 * @param cellLength The length of the cells.
	 * @param cutoff The cutoff radius.
	 */
	void rebuild(std::vector<CellTemplate> &cells,
				 const std::array<unsigned long, 3> &dims, double cellLength[3], double cutoff);

	void traverseCellPairs(CellProcessor &cellProcessor);

	void traverseCellPairs(traversalNames name, CellProcessor &cellProcessor);

	void traverseCellPairsOuter(CellProcessor &cellProcessor);

	void traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage, unsigned stageCount);


	bool isTraversalApplicable(traversalNames name, const std::array<unsigned long, 3> &dims) const; // new


	traversalNames getSelectedTraversal() const {
		return selectedTraversal;
	}

	CellPairTraversals<Cell> *getCurrentOptimalTraversal() { return _optimalTraversal; }

private:
	std::vector<CellTemplate>* _cells;
	std::array<unsigned long, 3> _dims;

	traversalNames selectedTraversal;

	std::vector<std::pair<CellPairTraversals<CellTemplate> *, CellPairTraversalData *> > _traversals;

	CellPairTraversals<CellTemplate> *_optimalTraversal;

	unsigned _cellsInCutoff = 1;
};

template<class CellTemplate>
TraversalTuner<CellTemplate>::TraversalTuner() : _cells(nullptr), _dims(), _optimalTraversal(nullptr) {
	// defaults:
	selectedTraversal =  ORIGINAL;

	auto *origData = new OriginalCellPairTraversalData;


	_traversals = {
			std::make_pair(nullptr, origData),
	};
#ifdef QUICKSCHED
	struct QuickschedTraversalData *quiData = new QuickschedTraversalData;
	quiData->taskBlockSize = {{2, 2, 2}};
	if (std::is_base_of<ParticleCellBase, CellTemplate>::value) {
		_traversals.push_back(std::make_pair(nullptr, quiData));
	}
#endif
}

template<class CellTemplate>
TraversalTuner<CellTemplate>::~TraversalTuner() {
	for (auto t : _traversals) {
		if (t.first != nullptr)
			delete (t.first);
		if (t.second != nullptr)
			delete (t.second);
	}
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::findOptimalTraversal() {
	// TODO implement autotuning here! At the moment the traversal is chosen via readXML!

	_optimalTraversal = _traversals[selectedTraversal].first;
	if (dynamic_cast<OriginalCellPairTraversal<CellTemplate> *>(_optimalTraversal))
			std::cout << "Using OriginalCellPairTraversal." << std::endl;

	// // log traversal
	// if (dynamic_cast<HalfShellTraversal<CellTemplate> *>(_optimalTraversal))
	// 	std::cout << "Using HalfShellTraversal." << std::endl;
	// else if (dynamic_cast<OriginalCellPairTraversal<CellTemplate> *>(_optimalTraversal))
	// 	std::cout << "Using OriginalCellPairTraversal." << std::endl;
// 	else if (dynamic_cast<C08CellPairTraversal<CellTemplate> *>(_optimalTraversal))
// 		Log::global_log->info() << "Using C08CellPairTraversal without eighthShell." << std::endl;
// 	else if (dynamic_cast<C08CellPairTraversal<CellTemplate, true> *>(_optimalTraversal))
// 		Log::global_log->info() << "Using C08CellPairTraversal with eighthShell." << std::endl;
// 	else if (dynamic_cast<C04CellPairTraversal<CellTemplate> *>(_optimalTraversal))
// 		Log::global_log->info() << "Using C04CellPairTraversal." << std::endl;
// 	else if (dynamic_cast<MidpointTraversal<CellTemplate> *>(_optimalTraversal))
// 		Log::global_log->info() << "Using MidpointTraversal." << std::endl;
// 	else if (dynamic_cast<NeutralTerritoryTraversal<CellTemplate> *>(_optimalTraversal))
// 		Log::global_log->info() << "Using NeutralTerritoryTraversal." << std::endl;
// 	else if (dynamic_cast<QuickschedTraversal<CellTemplate> *>(_optimalTraversal)) {
// 		Log::global_log->info() << "Using QuickschedTraversal." << std::endl;
// #ifndef QUICKSCHED
// 		Log::global_log->error() << "MarDyn was compiled without Quicksched Support. Aborting!" << std::endl;
// 		Simulation::exit(1);
// #endif
// 	} else if (dynamic_cast<SlicedCellPairTraversal<CellTemplate> *>(_optimalTraversal))
// 		Log::global_log->info() << "Using SlicedCellPairTraversal." << std::endl;
// 	else
// 		Log::global_log->warning() << "Using unknown traversal." << std::endl;
//
// 	if (_cellsInCutoff > _optimalTraversal->maxCellsInCutoff()) {
// 		Log::global_log->error() << "Traversal supports up to " << _optimalTraversal->maxCellsInCutoff()
// 							<< " cells in cutoff, but value is chosen as " << _cellsInCutoff << std::endl;
		throw std::runtime_error("traversaltuner error");
	}



template<class CellTemplate>
void TraversalTuner<CellTemplate>::rebuild(std::vector<CellTemplate> &cells, const std::array<unsigned long, 3> &dims,
										   double cellLength[3], double cutoff) {
	_cells = &cells; // new - what for?
	_dims = dims; // new - what for?

	for (size_t i = 0ul; i < _traversals.size(); ++i) {
		auto& [traversalPointerReference, traversalData] = _traversals[i];
		// decide whether to initialize or rebuild
		if (traversalPointerReference == nullptr) {
			switch (i) {
				case traversalNames ::ORIGINAL:
					traversalPointerReference = new OriginalCellPairTraversal<CellTemplate>(cells, dims);
					break;
				// case traversalNames::C08:
				// 	traversalPointerReference = new C08CellPairTraversal<CellTemplate>(cells, dims);
				// 	break;
				default:
					throw std::runtime_error("Unknown traversal data");
					// Log::global_log->error() << "Unknown traversal data found in TraversalTuner._traversals!" << std::endl;
					// Simulation::exit(1);
			}
		}
		traversalPointerReference->rebuild(cells, dims, cellLength, cutoff, traversalData);
	}
	_optimalTraversal = nullptr;
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairs(CellProcessor &cellProcessor) {
	if (not _optimalTraversal) {
		findOptimalTraversal();
	}
	_optimalTraversal->traverseCellPairs(cellProcessor);
}

template<class CellTemplate>
inline void TraversalTuner<CellTemplate>::traverseCellPairs(traversalNames name,
		CellProcessor& cellProcessor) {
	//if (name == getSelectedTraversal()) {
		traverseCellPairs(cellProcessor);
	// } else {
	// 	SlicedCellPairTraversal<CellTemplate> slicedTraversal(*_cells, _dims);
	// 	switch(name) {
	// 	case SLICED:
	// 		slicedTraversal.traverseCellPairs(cellProcessor);
	// 		break;
	// 	default:
	// 		Log::global_log->error()<< "Calling traverseCellPairs(traversalName, CellProcessor&) for something else than the Sliced Traversal is disabled for now. Aborting." << std::endl;
	// 		mardyn_exit(1);
	//
	// 		break;
	// 	}
	// }
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsOuter(CellProcessor &cellProcessor) {
	if (not _optimalTraversal) {
		findOptimalTraversal();
	}
	_optimalTraversal->traverseCellPairsOuter(cellProcessor);
}

template<class CellTemplate>
void TraversalTuner<CellTemplate>::traverseCellPairsInner(CellProcessor &cellProcessor, unsigned stage,
														  unsigned stageCount) {
	if (not _optimalTraversal) {
		findOptimalTraversal();
	}
	_optimalTraversal->traverseCellPairsInner(cellProcessor, stage, stageCount);
}

template<class CellTemplate>
inline bool TraversalTuner<CellTemplate>::isTraversalApplicable(
		traversalNames name, const std::array<unsigned long, 3> &dims) const {
	bool ret = true;
	switch(name) {
// 	case SLICED:
// 		ret = SlicedCellPairTraversal<CellTemplate>::isApplicable(dims);
// 		break;
// 	case QSCHED:
// #ifdef QUICKSCHED
// 		ret = true;
// #else
// 		ret = false;
// #endif
// 		break;
	// case C08:
	// 	ret = true;
	// 	break;
	// case C04:
	// 	ret = true;
	// 	break;
	case ORIGINAL:
		ret = true;
		break;
		default:
			throw std::runtime_error("unknown traversal ");
		//Log::global_log->warning() << "unknown traversal given in TraversalTuner::isTraversalApplicable, assuming that is applicable" << std::endl;
	}
	return ret;
}

