#pragma  once

#include <vector>
#include <array>

class CellProcessor;

struct CellPairTraversalData {
    virtual ~CellPairTraversalData() {}
};

template <class CellTemplate>
class CellPairTraversals {
public:
    CellPairTraversals(
        std::vector<CellTemplate>& cells,
        const std::array<unsigned long, 3>& dims): _cells(&cells), _dims(dims) {}

    virtual ~CellPairTraversals() {}

    /**
     * Reset all necessary data without reallocation.
     */
    virtual void rebuild(std::vector<CellTemplate>& cells, const std::array<unsigned long, 3>& dims,
                         double cellLength[3], double cutoff, CellPairTraversalData* data) {
        _cells = &cells;
        _dims = dims;
    };

    virtual void traverseCellPairs(CellProcessor& cellProcessor) = 0;
    virtual void traverseCellPairsOuter(CellProcessor& cellProcessor) = 0;
    virtual void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount) = 0;

    // @brief Should the domain decomposition exchange calculated forces at the boundaries,
    // or does this traversal calculate all forces.
    virtual bool requiresForceExchange() const {return false;}

    // @brief Returns the maximum number of cells per cutoff this traversal supports.
    virtual unsigned maxCellsInCutoff() const { return 1; }

protected:
    //TODO:
    //void traverseCellPairsNoDep(CellProcessor& cellProcessor);
    std::vector<CellTemplate> * _cells;
    std::array<unsigned long, 3> _dims;
};
