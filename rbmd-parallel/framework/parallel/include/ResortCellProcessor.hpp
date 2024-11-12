#pragma once

#include "CellProcessor.h"

class ResortCellProcessor : public CellProcessor {
public:
    ResortCellProcessor() : CellProcessor(0.0, 0.0) {}
    void initTraversal() override {}
    void preprocessCell(Cell& ) override {}

    void processCellPair(Cell& cell1, Cell& cell2, bool sumAll = false) override;

    void processCell(Cell& cell) override {}
    double processSingleMolecule(Atom*, Cell& ) override { return 0.0;}
    void postprocessCell(Cell& ) override {}
    void endTraversal() override {}

};

inline void ResortCellProcessor::processCellPair(Cell &cell1, Cell &cell2, bool sumAll) { // does this need a bool?
    cell1.UpdateLeavingMolecules(cell2);
}
