#pragma once

#include <vector>
#include <memory>
#include <variant>
#include "Atom.h"
#include "LinkedCellIterator.hpp"
#include "OpenMPUtil.hpp"
#include "RegionLinkedCellIterator.hpp"
#include "TraversalTuner.hpp"
static int update_step = 0;
class LinkedCell {
public:
    explicit LinkedCell(

    );
    ~LinkedCell();
    bool IsInvalidParticleReturner();
    RegionLinkedCellIterator RegionIterator(const double startRegion[3], const double endRegion[3],
                                            LinkedCellIterator::Type type);

    unsigned long int GetCellIndexOfPoint(const double point[3]) const;

    bool AddParticle(Atom& particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false, const bool& rebuildCaches = false);


    bool AddHaloParticle(Atom& particle, bool inBoxCheckedAlready = false, bool checkWhetherDuplicate = false,
                         const bool& rebuildCaches = false);

    std::variant<LinkedCellIterator, SingleCellIterator<Cell>> GetMoleculeAtPosition(const double pos[3]) ;

    void deleteMolecule(LinkedCellIterator &moleculeIter, const bool& rebuildCaches);

    double GetCutoff() const ;

    std::vector<Atom>& GetInvalidParticlesRef() ;

    double getBoundingBoxMin(int dimension) const;

    double getBoundingBoxMax(int dimension) const ;

    bool isInBoundingBox(double r[3]) const;

    unsigned long getNumberOfParticles(LinkedCellIterator::Type t = LinkedCellIterator::ALL_CELLS);

    void addParticles(std::vector<Atom>& particles, bool checkWhetherDuplicate=false) ;

    void DeleteOuterParticles();




    bool rebuild(double bBoxMin[3], double bBoxMax[3]);

    void update();

    void SetCutoff(double rc)  { _cutoffRadius = rc; }

    LinkedCellIterator iterator (LinkedCellIterator::Type t) {
        const LinkedCellIterator::CellIndex_T offset = GetThreadNum();
        const LinkedCellIterator::CellIndex_T stride = GetNumThreads();

        return LinkedCellIterator(t, &_cells, offset, stride);
    }



    bool IsInBoundingBox(double r[3]) const;

    std::vector<Cell> &GetCells();
    int GetCellPerDim(int d) const;
    
    double _haloBoundingBoxMin[3]{}; //!< low corner of the bounding box around the linked cells (including halo)
    double _haloBoundingBoxMax[3]{}; //!< high corner of the bounding box around the linked cells (including halo)
private:
    std::vector<Cell> _cells{}; //!< Vector containing all cells (including halo)

    std::vector<unsigned long> _haloCellIndices{}; //!< Vector containing the indices of all halo cells in the _cells vector

    // std::unique_ptr<TraversalTuner<Atom>> _traversalTuner; // new TODO?

    int _cellsPerDimension[3]{}; //!< Number of Cells in each spatial dimension (including halo)
    int _haloWidthInNumCells[3]{}; //!< Halo width (in cells) in each dimension
    int _boxWidthInNumCells[3]{}; //!< Box width (in cells) in each dimension
    double _haloLength[3]{}; //!< width of the halo strip (in size units)
    double _cellLength[3]{}; //!< length of the cell (for each dimension)
    double _cellLengthReciprocal[3]{}; //!< 1.0 / _cellLength, to speed-up particle sorting
    double _cutoffRadius{}; //!< RDF/electrostatics cutoff radius
    unsigned _cellsInCutoff = 1; //!< Cells in cutoff radius -> cells with size cutoff / cellsInCutoff
    double _boundingBoxMin[3]{};
    //! coordinates of the right, upper, back corner of the bounding box
    double _boundingBoxMax[3]{};
    std::vector<Atom> _invalidParticles{};
    long int CellIndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const ;


    void InitializeCells();

    void deleteParticlesOutsideBox(double boxMin[3], double boxMax[3]) ;
    //! @brief True if all Particles are in the right cell
    //!
    //! The particles themselves are not stored in cells, but in one large
    //! list. The cells only contain pointers to particles. But at some
    //! points during the simulation, Those pointers are invalid. This happens
    //! e.g. after the integrator has changed the positions of particles.
    //! If this happens, no method must be able to access the particles via
    //! the cells. Therefore, whenever a piece of code causes the cells to
    //! become possibly invalid, _cellsValid has to be set to false. Methods
    //! accessing cells have to check whether _cellsValid is true (and e.g.
    //! abort the program if not). After the cells are updated, _cellsValid
    //! should be set to true.
    bool _cellsValid{};
    void CalculateNeighbourIndices(std::vector<long>& forwardNeighbourOffsets, std::vector<long>& backwardNeighbourOffsets) const;

    void GetCellIndicesOfRegion(const double startRegion[3], const double endRegion[3], unsigned int &startRegionCellIndex, unsigned int &endRegionCellIndex);

    void ThreeDIndexOfCellIndex(int ind, int r[3], const int dim[3]) const;

    unsigned long int GetCellIndexOfMolecule(Atom* molecule) const;
    void CheckAtomsInBox();
    // TODO 3.21 暂时不做这个  update_via_traversal目前来看好像可以替代
    void InitializeTraversal();

    std::unique_ptr<TraversalTuner<Cell>> _traversalTuner;
    void UpdateViaTraversal();
    std::array<std::pair<unsigned long, unsigned long>, 14> calculateCellPairOffsets() const;
    void UpdateViaCopies();

    void UpdateViaCopiesLongRange();

};
