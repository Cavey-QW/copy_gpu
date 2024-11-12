
#pragma once

#include <cstddef>
#include <cmath>

#include "Cell.h"

class CellProcessor {
protected:
	double _cutoffRadiusSquare;
	double _LJCutoffRadiusSquare;

public:
	CellProcessor(const double cutoffRadius, const double LJCutoffRadius) :
		_cutoffRadiusSquare(cutoffRadius * cutoffRadius),
		_LJCutoffRadiusSquare(LJCutoffRadius * LJCutoffRadius) {}
    /**
     * virtual destructor
     */
	virtual ~CellProcessor() {}

	double getCutoffRadius() const {return sqrt(_cutoffRadiusSquare);}
	double getLJCutoffRadius() const {return sqrt(_LJCutoffRadiusSquare);}
	void setCutoffRadius(const double c) {_cutoffRadiusSquare = c * c;}
	void setLJCutoffRadius(const double ljc) {_LJCutoffRadiusSquare = ljc * ljc;}


	double getCutoffRadiusSquare() const {return _cutoffRadiusSquare;}
	double getLJCutoffRadiusSquare() const {return _LJCutoffRadiusSquare;}
	void setCutoffRadiusSquare(const double c) {_cutoffRadiusSquare = c;}
	void setLJCutoffRadiusSquare(const double ljc) {_LJCutoffRadiusSquare = ljc;}


	/**
	 * called before the traversal starts.
	 *
	 * @param numCells number of cells in window
	 */
	virtual void initTraversal() = 0;

	/**
	 * Called before a cell is touched for the first time during an interation.
	 */
	virtual void preprocessCell(Cell& cell) = 0;

	/**
	 * Called for each cell pair within the cutoff radius. Called exactly once per
	 * pair (i.e. pairs are not ordered).
	 *
	 * @note will not be called for empty cells.
	 * Sum up all macroscopic values (e.g. for hs) or only half of them (e.g. for fs)
         */
	virtual void processCellPair(Cell& cell1, Cell& cell2, bool sumAll = false) = 0;

	/**
	 * Called when this cell is the current cell.
	 *
	 * @note will not be called for empty cells.
	 */
	virtual void processCell(Cell& cell) = 0;

	virtual double processSingleMolecule(Atom* m1, Cell& cell2) = 0;

	/**
	 * Called after the cell has been considered for the last time during the traversal.
	 */
	virtual void postprocessCell(Cell& cell) = 0;

	/**
	 * Called after the traversal finished.
	 */
	virtual void endTraversal() = 0;
};

