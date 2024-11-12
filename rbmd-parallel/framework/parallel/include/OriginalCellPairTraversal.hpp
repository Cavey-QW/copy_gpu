#pragma  once

#include "LogTool.h"
#include <stdexcept>

#include "CellProcessor.h"
#include "IdexMappingTool.hpp"
#include "CellPairTraversalData.hpp"

struct OriginalCellPairTraversalData : CellPairTraversalData {
};

template <class CellTemplate>
class OriginalCellPairTraversal: public CellPairTraversals<CellTemplate> {
public:
	OriginalCellPairTraversal(std::vector<CellTemplate> &cells, const std::array<unsigned long, 3> &dims) :
			CellPairTraversals<CellTemplate>(cells, dims) {
		computeNeighbourOffsets();
	}

	virtual ~OriginalCellPairTraversal() {}

	virtual void rebuild(std::vector<CellTemplate>& cells,
						 const std::array<unsigned long, 3>& dims, double cellLength[3], double cutoff,
						 struct CellPairTraversalData *data);

	void traverseCellPairs(CellProcessor& cellProcessor);
	void traverseCellPairsOuter(CellProcessor& cellProcessor);
	void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount);

	virtual void processBaseCell(CellProcessor& cellProcessor, unsigned long cellIndex) const;

protected:
	// couldn't get it to compile when the enum is part of the class, so making it global...
	enum TraverseType {
		ALL_CELLS,
		INNER_CELLS,
		OUTER_CELLS
	};

	void computeNeighbourOffsets();
	void traverseCellPairsBackend(CellProcessor& cellProcessor, unsigned loIndex, unsigned hiIndex, TraverseType type) const;


	std::array<long, 13> _forwardNeighbourOffsets; //!< Neighbours that come in the total ordering after a cell
	std::array<long, 13> _backwardNeighbourOffsets; //!< Neighbours that come in the total ordering before a cell
	std::vector<unsigned long> _innerMostCellIndices; //!< Vector containing the indices (for the cells vector) of all inner cells (without boundary)
};

template<class CellTemplate>
void OriginalCellPairTraversal<CellTemplate>::rebuild(std::vector<CellTemplate> &cells,
													  const std::array<unsigned long, 3> &dims, double cellLength[3], double cutoff,
													  struct CellPairTraversalData *data) {
	if (dynamic_cast<OriginalCellPairTraversalData *>(data)) {
		CellPairTraversals<CellTemplate>::rebuild(cells, dims, cellLength, cutoff, data);
		computeNeighbourOffsets();

		_innerMostCellIndices.clear();

		auto maxIndex = 1;
		for (auto d : dims)
			maxIndex *= d;

		for (auto i = 0; i < maxIndex; ++i) {
			if (this->_cells->at(i).IsInnerMostCell()){
				_innerMostCellIndices.push_back(i);
			}
		}
	} else {
		// Log::global_log->error() << "OriginalCellPairTraversalDat::rebuild was called with incompatible Traversal data!" << std::endl;
		// Simulation::exit(-1);
		throw std::runtime_error("ebuild was called with incompatible Traversal data!");
	}
}

template<class CellTemplate>
inline void OriginalCellPairTraversal<CellTemplate>::computeNeighbourOffsets() {
	// Log::global_log->debug() << "Setting up cell neighbour index lists." << std::endl;

	std::fill(_forwardNeighbourOffsets.begin(), _forwardNeighbourOffsets.end(), 0);
	std::fill(_backwardNeighbourOffsets.begin(), _backwardNeighbourOffsets.end(), 0);
	int forwardNeighbourIndex = 0, backwardNeighbourIndex = 0;

	long maxNeighbourOffset = 0;
	long minNeighbourOffset = 0;

	std::array<long, 3> dims;
	for (int d = 0; d < 3; ++d) {
		dims[d] = static_cast<long>(this->_dims[d]);
	}

	std::array<long, 3> r;
	for (r[2] = -1; r[2] <= 1; r[2]++) {
		for (r[1] = -1; r[1] <= 1; r[1]++) {
			for (r[0] = -1; r[0] <= 1; r[0]++) {

				long offset = ThreeToOneD(r, dims);

				if (offset > 0) {
					_forwardNeighbourOffsets[forwardNeighbourIndex] = offset;
					++forwardNeighbourIndex;
					if (offset > maxNeighbourOffset) {
						maxNeighbourOffset = offset;
					}
				}
				if (offset < 0) {
					_backwardNeighbourOffsets[backwardNeighbourIndex] = abs(offset);
					++backwardNeighbourIndex;
					if (abs(offset) > minNeighbourOffset) {
						minNeighbourOffset = abs(offset);
					}
				}
			}
		}
	}

//	assert(forwardNeighbourIndex == 13);
//	assert(backwardNeighbourIndex == 13);
    if (forwardNeighbourIndex!=13 || backwardNeighbourIndex != 13){
        GlobalLogger::error("NeighbourIndex Error In computeNeighbourOffsets");
    }
	// Log::global_log->info() << "Neighbour offsets are bounded by "
	// 		<< minNeighbourOffset << ", " << maxNeighbourOffset << std::endl;

}

template<class CellTemplate>
inline void OriginalCellPairTraversal<CellTemplate>::traverseCellPairsBackend (
		CellProcessor& cellProcessor, unsigned loIndex, unsigned hiIndex,
		TraverseType type) const {
	switch(type) {
	case ALL_CELLS:
		for (unsigned i = loIndex; i < hiIndex; ++i) {
			processBaseCell(cellProcessor, i);
		}
		break;
	case INNER_CELLS:
		for (unsigned i = loIndex; i < hiIndex; ++i) {
			processBaseCell(cellProcessor, _innerMostCellIndices.at(i));
		}
		break;
	case OUTER_CELLS:
		for (unsigned i = loIndex; i < hiIndex; ++i) {
			CellTemplate& baseCell = this->_cells->at(i);
			if (!baseCell.IsInnerMostCell()) {
				processBaseCell(cellProcessor, i);
			}
		}
		break;
	}
}

template<class CellTemplate>
inline void OriginalCellPairTraversal<CellTemplate>::traverseCellPairs(CellProcessor& cellProcessor) {
	unsigned long start = 0ul;
	unsigned long end = this->_cells->size();
	traverseCellPairsBackend(cellProcessor, start, end, ALL_CELLS);
}

template<class CellTemplate>
inline void OriginalCellPairTraversal<CellTemplate>::traverseCellPairsOuter(CellProcessor& cellProcessor) {
	unsigned long start = 0ul;
	unsigned long end = this->_cells->size();
	traverseCellPairsBackend(cellProcessor, start, end, OUTER_CELLS);
}

template<class CellTemplate>
inline void OriginalCellPairTraversal<CellTemplate>::traverseCellPairsInner(CellProcessor& cellProcessor,
																			unsigned stage,
																			unsigned stageCount) {
	unsigned long start =  _innerMostCellIndices.size() * stage / stageCount;
	unsigned long end =  _innerMostCellIndices.size() * (stage+1) / stageCount;
	traverseCellPairsBackend(cellProcessor, start, end, INNER_CELLS);
}

template<class CellTemplate>
inline void OriginalCellPairTraversal<CellTemplate>::processBaseCell(CellProcessor& cellProcessor,
																	 unsigned long cellIndex) const {

	CellTemplate& currentCell = this->_cells->at(cellIndex);

	if (currentCell.IsInnerCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto& neighbourOffset : this->_forwardNeighbourOffsets) {
			CellTemplate& neighbourCell = this->_cells->at(cellIndex + neighbourOffset);
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	} else
	// loop over all boundary cells and calculate forces to forward and backward neighbours
	if (currentCell.IsBoundaryCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto& neighbourOffset : this->_forwardNeighbourOffsets) {
			CellTemplate& neighbourCell = this->_cells->at(cellIndex + neighbourOffset);
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
		// loop over all backward neighbours. calculate only forces
		// to neighbour cells in the halo region, all others already have been calculated
		for (auto& neighbourOffset : this->_backwardNeighbourOffsets) {
			CellTemplate& neighbourCell = this->_cells->at(cellIndex - neighbourOffset);
			if (neighbourCell.IsHaloCell()) {
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}
	} // if ( isBoundaryCell() )
}






