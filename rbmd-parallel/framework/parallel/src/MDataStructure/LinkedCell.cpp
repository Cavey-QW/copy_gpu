#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>


#include "LinkedCell.h"

//#include <unistd.h>
#include <mpi.h>

#include "CellManager.h"
#include "ResortCellProcessor.hpp"

#include <RBMDParallelUtil.h>
#include <algorithm>
#include <csignal>

bool LinkedCell::IsInvalidParticleReturner()
{
  return false;
}

unsigned long int LinkedCell::GetCellIndexOfPoint(const double* point) const
{
  int cellIndex[3]; // 3D Cell index
  double localPoint[] = { point[0], point[1], point[2] };

  for (int dim = 0; dim < 3; dim++)
  {
    // different than getCellIndexOfMolecule!!!

    // ignore a bit of rounding, if the point is outside of the box.
    if (localPoint[dim] <= _haloBoundingBoxMin[dim])
    {
      localPoint[dim] += _cellLength[dim] * 0.5;
    }
    else if (localPoint[dim] >= _haloBoundingBoxMax[dim])
    {
      localPoint[dim] -= _cellLength[dim] * 0.5;
    }

#ifndef NDEBUG
    //this should never ever happen!
    if (localPoint[dim] < _haloBoundingBoxMin[dim] || localPoint[dim] >= _haloBoundingBoxMax[dim])
    {
      //            Log::global_log->error() << "Point is outside of halo bounding box" << std::endl;
      //            Log::global_log->error() << "Point p = (" << localPoint[0] << ", " << localPoint[1] << ", " << localPoint[2] << ")" << std::endl;
      //            Log::global_log->error() << "_haloBoundingBoxMin = (" << _haloBoundingBoxMin[0] << ", " << _haloBoundingBoxMin[1] << ", " << _haloBoundingBoxMin[2] << ")" << std::endl;
      //            Log::global_log->error() << "_haloBoundingBoxMax = (" << _haloBoundingBoxMax[0] << ", " << _haloBoundingBoxMax[1] << ", " << _haloBoundingBoxMax[2] << ")" << std::endl;
      GlobalLogger::error("Point is outside of halo bounding box");
    }
#endif

    //this version is sensitive to roundoffs, if we have molecules (initialized) precisely at position 0.0:
    //cellIndex[dim] = (int) floor(point[dim] - _haloBoundingBoxMin[dim]) * _cellLengthReciprocal[dim]);
    cellIndex[dim] = std::min(
      std::max(((int)floor((localPoint[dim] - _boundingBoxMin[dim]) * _cellLengthReciprocal[dim])) +
                 _haloWidthInNumCells[dim],
               0),
      _cellsPerDimension[dim] - 1);
  }

  int cellIndex1d = this->CellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
  // in very rare cases rounding is stupid, thus we need a check...
  //TODO: check if this can in any way be done better...
  if (_cells[cellIndex1d].TestPointInCell(localPoint))
  {
    return cellIndex1d;
  }
  else
  {
    for (int dim = 0; dim < 3; dim++)
    {
      if (localPoint[dim] < _cells[cellIndex1d].GetBoundingBoxMin(dim))
      {
        cellIndex[dim]--;
      }
      else
      {
        if (localPoint[dim] >= _cells[cellIndex1d].GetBoundingBoxMax(dim))
        {
          cellIndex[dim]++;
        }
      }
    }
    cellIndex1d = this->CellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
    if (!_cells[cellIndex1d].TestPointInCell(localPoint))
    {
      throw std::runtime_error("Test Point Failed");
    }
    return cellIndex1d;
  }
}


void LinkedCell::GetCellIndicesOfRegion(const double* startRegion,
                                        const double* endRegion,
                                        unsigned int& startRegionCellIndex,
                                        unsigned int& endRegionCellIndex)
{
  startRegionCellIndex = GetCellIndexOfPoint(startRegion);
  endRegionCellIndex = GetCellIndexOfPoint(endRegion);
}

long int LinkedCell::CellIndexOf3DIndex(long xIndex, long yIndex, long zIndex) const
{
  return (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
}


void LinkedCell::ThreeDIndexOfCellIndex(int ind, int* r, const int* dim) const
{
  r[2] = ind / (dim[0] * dim[1]);
  r[1] = (ind - r[2] * dim[0] * dim[1]) / dim[0];
  r[0] = ind - dim[0] * (r[1] + dim[1] * r[2]);
}

RegionLinkedCellIterator LinkedCell::RegionIterator(const double* startRegion,
                                                    const double* endRegion,
                                                    LinkedCellIterator::Type type)
{
  // parameter "type" not yet used
  // add functionality in a future version...
  unsigned int startRegionCellIndex;
  unsigned int endRegionCellIndex;

  // only include halo cells if iterator explicitly asks for all cells
  const auto& localBoxOfInterestMin =
    type == LinkedCellIterator::ALL_CELLS ? _haloBoundingBoxMin : _boundingBoxMin;
  const auto& localBoxOfInterestMax =
    type == LinkedCellIterator::ALL_CELLS ? _haloBoundingBoxMax : _boundingBoxMax;

  // clamp iterated region to local MPI subdomain
  const std::array<double, 3> startRegionClamped = {
    std::clamp(startRegion[0], localBoxOfInterestMin[0], localBoxOfInterestMax[0]),
    std::clamp(startRegion[1], localBoxOfInterestMin[1], localBoxOfInterestMax[1]),
    std::clamp(startRegion[2], localBoxOfInterestMin[2], localBoxOfInterestMax[2]),
  };
  const std::array<double, 3> endRegionClamped = {
    std::clamp(endRegion[0], localBoxOfInterestMin[0], localBoxOfInterestMax[0]),
    std::clamp(endRegion[1], localBoxOfInterestMin[1], localBoxOfInterestMax[1]),
    std::clamp(endRegion[2], localBoxOfInterestMin[2], localBoxOfInterestMax[2]),
  };

  // check if the clipping resulted in the region of interest box collapsing as it is not part of the local domain.
  const auto localVolumeIsZero = localBoxOfInterestMin[0] == localBoxOfInterestMax[0] or
    localBoxOfInterestMin[1] == localBoxOfInterestMax[1] or
    localBoxOfInterestMin[2] == localBoxOfInterestMax[2];
  if (localVolumeIsZero)
  {
    // return invalid iterator (_cells == nullptr)
    return RegionLinkedCellIterator{};
  }

  GetCellIndicesOfRegion(
    startRegionClamped.data(), endRegionClamped.data(), startRegionCellIndex, endRegionCellIndex);

  std::array<int, 3> start3DIndices{}, end3DIndices{};
  ThreeDIndexOfCellIndex(
    static_cast<int>(startRegionCellIndex), start3DIndices.data(), _cellsPerDimension);
  ThreeDIndexOfCellIndex(
    static_cast<int>(endRegionCellIndex), end3DIndices.data(), _cellsPerDimension);
  const std::array<int, 3> regionDimensions = {
    end3DIndices[0] - start3DIndices[0] + 1,
    end3DIndices[1] - start3DIndices[1] + 1,
    end3DIndices[2] - start3DIndices[2] + 1,
  };
  // if the iterator on this rank has nothing to iterate over invalidate it by pushing its cell offset out of range.
  const LinkedCellIterator::CellIndex_T offset =
    GetThreadNum(); // starting position         //TODO 检查是否错误
  const LinkedCellIterator::CellIndex_T stride = GetNumThreads(); // stride

  return { type,
           &_cells,
           offset,
           stride,
           static_cast<int>(startRegionCellIndex),
           regionDimensions.data(),
           _cellsPerDimension,
           startRegion,
           endRegion };
}


unsigned long int LinkedCell::GetCellIndexOfMolecule(Atom* molecule) const
{
  double r[3] = { molecule->r(0), molecule->r(1), molecule->r(2) };
  return GetCellIndexOfPoint(r);
}


bool LinkedCell::AddParticle(Atom& particle,
                             bool inBoxCheckedAlready,
                             bool checkWhetherDuplicate,
                             const bool& rebuildCaches)
{
  bool wasInserted = false;
  const bool inBox =
    inBoxCheckedAlready or particle.InBox(_haloBoundingBoxMin, _haloBoundingBoxMax);
  if (inBox)
  {
    int cellIndex = GetCellIndexOfMolecule(&particle);
    wasInserted =
      _cells[cellIndex].AddAtomToCell(particle); //TODO check Temporarily Not Implemented
    //        if (rebuildCaches) {  TODO do not hava soa current
    //            _cells[cellIndex].buildSoACaches();
    //        }
  }
  return wasInserted;
}


bool LinkedCell::AddHaloParticle(Atom& particle,
                                 bool inBoxCheckedAlready,
                                 bool checkWhetherDuplicate,
                                 const bool& rebuildCaches)
{
  //    if(particle.InBox(_boundingBoxMin, _boundingBoxMax)){
  //        std::cout << "Bounding box in domian: "<< " : " << "[" << _boundingBoxMin[0] << ", " << _boundingBoxMax[0] << "]" << " x " << "["
  //                  << _boundingBoxMin[1] << ", " << _boundingBoxMax[1] << "]" << " x " << "[" << _boundingBoxMin[2] << ", " << _boundingBoxMax[2] << "]"
  //                  << std::endl;
  //        std::cout << "x:" << particle.r(0) <<",y:" << particle.r(1) << ",z:" << particle.r(2) << std::endl;
  //
  //    }
  //    ASSERT(not particle.InBox(_boundingBoxMin,_boundingBoxMax));
  if (particle.r(0) == 0 && particle.r(1) == 0 && particle.r(2) == 0)
  {
    Logger::info("error zero xyz halo get");
    return false;
  }
  return AddParticle(particle, false, checkWhetherDuplicate, rebuildCaches);
}

std::variant<LinkedCellIterator, SingleCellIterator<Cell>> LinkedCell::GetMoleculeAtPosition(
  const double* pos)
{
  const double epsi = this->_cutoffRadius * 1e-6;
  auto index = GetCellIndexOfPoint(pos);
  auto& cell = _cells.at(index);

  // iterate through cell and compare position of molecules with given position


  for (auto cellIterator = cell.iterator(); cellIterator.isValid(); ++cellIterator)
  {
    auto& mol = *cellIterator;

    if (fabs(cellIterator->r(0) - pos[0]) <= epsi && fabs(cellIterator->r(1) - pos[1]) <= epsi &&
        fabs(cellIterator->r(2) - pos[2]) <= epsi)
    {
      // found
      return cellIterator;
    }
  }
  // not found -> return default initialized iter.
  return {};
}


double LinkedCell::GetCutoff() const
{
  return _cutoffRadius;
}

void LinkedCell::deleteMolecule(LinkedCellIterator& moleculeIter, const bool& rebuildCaches)
{
  moleculeIter.deleteCurrentParticle();

  if (rebuildCaches)
  {
    auto cellid = GetCellIndexOfMolecule(&*moleculeIter);
    if (cellid >= _cells.size())
    {
      //            Log::global_log->error_always_output()
      //                    << "coordinates for atom deletion lie outside bounding box."
      //                    << std::endl;
      exit(1);
    }
    // _cells[cellid].buildSoACaches();  //暂未实现
  }
}

std::vector<Atom>& LinkedCell::GetInvalidParticlesRef()
{
  return _invalidParticles;
}

double LinkedCell::getBoundingBoxMin(int dimension) const
{
  return this->_boundingBoxMin[dimension];
}

double LinkedCell::getBoundingBoxMax(int dimension) const
{
  return this->_boundingBoxMax[dimension];
}

bool LinkedCell::isInBoundingBox(double* r) const
{
  if (r[0] >= _boundingBoxMin[0] && r[1] >= _boundingBoxMin[1] && r[2] >= _boundingBoxMin[2] &&
      r[0] < _boundingBoxMax[0] && r[1] < _boundingBoxMax[1] && r[2] < _boundingBoxMax[2])
  {
    return true;
  }
  else
  {
    return false;
  }
}

unsigned long LinkedCell::getNumberOfParticles(LinkedCellIterator::Type t)
{
  unsigned long N = 0;
  unsigned long numCells = _cells.size();

#if defined(_OPENMP)
#pragma omp parallel for reduction(+ : N)
#endif
  for (unsigned long i = 0; i < numCells; ++i)
  {
    if ((t == LinkedCellIterator::ALL_CELLS) or (not _cells.at(i).IsHaloCell()))
    {
      N += _cells.at(i).GetAtomsCount();
    }
  }
  return N;
}

void LinkedCell::addParticles(std::vector<Atom>& particles, bool checkWhetherDuplicate)
{
  typedef std::vector<Atom>::size_type mol_index_t;
  typedef std::vector<Cell>::size_type cell_index_t;

#ifndef NDEBUG
  int oldNumberOfParticles = getNumberOfParticles();
#endif

  const mol_index_t N = particles.size();

  std::map<cell_index_t, std::vector<mol_index_t>> newPartsPerCell;

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
    std::map<cell_index_t, std::vector<mol_index_t>> local_newPartsPerCell;

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (mol_index_t i = 0; i < N; ++i)
    {
      Atom& particle = particles[i];

#ifndef NDEBUG
      if (!particle.InBox(_haloBoundingBoxMin, _haloBoundingBoxMax))
      {
        //Log::global_log->error()<<"At particle with ID "<< particle.getID()<<" assertion failed..."<<std::endl;
      }
      rbmd_assert(particle.InBox(_haloBoundingBoxMin, _haloBoundingBoxMax));
#endif

      const unsigned long cellIndex = GetCellIndexOfMolecule(&particle);
      rbmd_assert(cellIndex < _cells.size());
      local_newPartsPerCell[cellIndex].push_back(i);
    }

#if defined(_OPENMP)
#pragma omp critical(add_particles_reduce_maps)
#endif
    {
      for (auto it = local_newPartsPerCell.begin(); it != local_newPartsPerCell.end(); ++it)
      {
        cell_index_t cellIndex = it->first;
        std::vector<mol_index_t>& global_vector = newPartsPerCell[cellIndex];
        std::vector<mol_index_t>& local_vector = it->second;
        global_vector.insert(global_vector.end(), local_vector.begin(), local_vector.end());
      }
    }

#if defined(_OPENMP)
#pragma omp barrier
#endif

    const cell_index_t numCells = newPartsPerCell.size();
    const int thread_id = GetThreadNum();
    const int num_threads = GetNumThreads();
    const cell_index_t my_cells_start = (thread_id)*numCells / num_threads;
    const cell_index_t my_cells_end = (thread_id + 1) * numCells / num_threads;

    auto thread_begin = newPartsPerCell.begin();
    advance(thread_begin, my_cells_start);
    auto thread_end = newPartsPerCell.begin();
    advance(thread_end, my_cells_end);

    for (auto it = thread_begin; it != thread_end; ++it)
    {
      cell_index_t cellIndex = it->first;
      const std::vector<mol_index_t>& global_vector = it->second;

      const size_t numMolsInCell = global_vector.size();
      _cells[cellIndex].increaseMoleculeStorage(numMolsInCell);

      for (size_t j = 0; j < numMolsInCell; ++j)
      {
        const mol_index_t molIndex = global_vector[j];
        Atom& mol = particles[molIndex];
        _cells[cellIndex].AddAtomToCell(mol);
        // AddAtomToCell(mol, checkWhetherDuplicate);  原来：addParticle(mol, checkWhetherDuplicate);
      }
    }
  } // end pragma omp parallel

#ifndef NDEBUG
  int numberOfAddedParticles = getNumberOfParticles() - oldNumberOfParticles;
  //    Log::global_log->debug()<<"In LinkedCells::addParticles :"<<std::endl;
  //    Log::global_log->debug()<<"\t#Particles to be added = "<<particles.size()<<std::endl;
  //    Log::global_log->debug()<<"\t#Particles actually added = "<<numberOfAddedParticles<<std::endl;
#endif
}

void LinkedCell::DeleteOuterParticles()
{
  const size_t numHaloCells = _haloCellIndices.size();

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < numHaloCells; i++)
  {
    Cell& currentCell = _cells[_haloCellIndices[i]];
    currentCell.EmptyAtoms();
  }
}


void LinkedCell::InitializeCells()
{
  //TODO: resize _haloCellIndices
  _haloCellIndices.clear();

  long int cellIndex = 0;

  CellManager::Instance().init(_cellsPerDimension,
                               _haloBoundingBoxMin,
                               _haloBoundingBoxMax,
                               _boundingBoxMin,
                               _boundingBoxMax,
                               _cellLength,
                               _haloWidthInNumCells);


  for (int iz = 0; iz < _cellsPerDimension[2]; ++iz)
  {
    for (int iy = 0; iy < _cellsPerDimension[1]; ++iy)
    {
      for (int ix = 0; ix < _cellsPerDimension[0]; ++ix)
      {
        cellIndex = CellIndexOf3DIndex(ix, iy, iz);
        Cell& cell = _cells[cellIndex];
        cell.SetID(cellIndex); //set the index of the cell to the index of it...

        if (ix < _haloWidthInNumCells[0] || iy < _haloWidthInNumCells[1] ||
            iz < _haloWidthInNumCells[2] || ix >= _cellsPerDimension[0] - _haloWidthInNumCells[0] ||
            iy >= _cellsPerDimension[1] - _haloWidthInNumCells[1] ||
            iz >= _cellsPerDimension[2] - _haloWidthInNumCells[2])
        {
          cell.AssignCellToHaloRegion();
          _haloCellIndices.push_back(cellIndex); //TODO 2024年4月3日 没有问题 对比调试
        }
      }
    }
  }
}


void LinkedCell::deleteParticlesOutsideBox(double boxMin[3], double boxMax[3])
{
  // This should be unimportant
#if defined(_OPENMP)
#pragma omp parallel
#endif
  for (auto it = iterator(LinkedCellIterator::ALL_CELLS); it.isValid(); ++it)
  {
    bool outside = not it->InBox(boxMin, boxMax);
    if (outside)
    {
      it.deleteCurrentParticle();
    }
  }
}

// TODO 3.21 暂时不做这个
// void LinkedCell::InitializeTraversal() {
//     std::array<long unsigned, 3> dims;
//     for (int d = 0; d < 3; ++d) {
//         dims[d] = _cellsPerDimension[d];
//     }
//     _traversalTuner->rebuild(_cells, dims, _cellLength, _cutoffRadius);
// }

void LinkedCell::CheckAtomsInBox()
{
  std::vector<Atom> badMolecules;
  unsigned numBadMolecules = 0;

#if defined(_OPENMP)
#pragma omp parallel reduction(+ : numBadMolecules)
#endif
  {
    for (LinkedCellIterator tM = iterator(LinkedCellIterator::ALL_CELLS); tM.isValid(); ++tM)
    {
      if (not tM->InBox(_haloBoundingBoxMin, _haloBoundingBoxMax))
      {
        numBadMolecules++;

#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          badMolecules.push_back(*tM);
        }
      }
    }
  }

  if (numBadMolecules > 0)
  {
    std::cout << "Found " << numBadMolecules << "Atom  outside of bounding box:" << std::endl;
    for (auto& m : badMolecules)
    {
      std::cout << "Atom (id=" << m.GetID() << "), (current position: x=" << m.r(0)
                << ", y=" << m.r(1) << ", z=" << m.r(2) << ")" << std::endl;
    }
    std::cout << "The bounding box is: [" << _haloBoundingBoxMin[0] << ", "
              << _haloBoundingBoxMax[0] << ") x [" << _haloBoundingBoxMin[1] << ", "
              << _haloBoundingBoxMax[1] << ") x [" << _haloBoundingBoxMin[2] << ", "
              << _haloBoundingBoxMax[2] << ")" << std::endl;
    std::cout << "Particles will be lost. Aborting simulation." << std::endl;
    throw std::runtime_error("error in LinkedCell updating");
  }
}

void LinkedCell::InitializeTraversal()
{
  std::array<long unsigned, 3> dims;
  for (int d = 0; d < 3; ++d)
  {
    dims[d] = _cellsPerDimension[d];
  }
  _traversalTuner->rebuild(_cells, dims, _cellLength, _cutoffRadius);
}

void LinkedCell::UpdateViaTraversal()
{
  auto resortCellProcessor = new ResortCellProcessor();
  _traversalTuner->traverseCellPairs(*resortCellProcessor);
}

void LinkedCell::CalculateNeighbourIndices(std::vector<long>& forwardNeighbourOffsets,
                                           std::vector<long>& backwardNeighbourOffsets) const
{
  //Log::global_log->debug() << "Setting up cell neighbour indice lists." << std::endl;

  // 13 neighbors for _haloWidthInNumCells = 1 or 64 for =2
  int maxNNeighbours = ((2 * _haloWidthInNumCells[0] + 1) * (2 * _haloWidthInNumCells[1] + 1) *
                          (2 * _haloWidthInNumCells[2] + 1) -
                        1) /
    2;

  // Resize offset vector to number of neighbors and fill with 0
  forwardNeighbourOffsets.reserve(maxNNeighbours);
  backwardNeighbourOffsets.reserve(maxNNeighbours);

  int forwardNeighbourIndex = 0, backwardNeighbourIndex = 0;

  double xDistanceSquare;
  double yDistanceSquare;
  double zDistanceSquare;
  double cutoffRadiusSquare = pow(_cutoffRadius, 2);
  for (int zIndex = -_haloWidthInNumCells[2]; zIndex <= _haloWidthInNumCells[2]; zIndex++)
  {
    // The distance in one dimension is the width of a cell multiplied with the number
    // of cells between the two cells (this is received by subtracting one of the
    // absolute difference of the cells, if this difference is not zero)
    if (zIndex != 0)
    {
      zDistanceSquare = pow((abs(zIndex) - 1) * _cellLength[2], 2);
    }
    else
    {
      zDistanceSquare = 0;
    }
    for (int yIndex = -_haloWidthInNumCells[1]; yIndex <= _haloWidthInNumCells[1]; yIndex++)
    {
      if (yIndex != 0)
      {
        yDistanceSquare = pow((abs(yIndex) - 1) * _cellLength[1], 2);
      }
      else
      {
        yDistanceSquare = 0;
      }
      for (int xIndex = -_haloWidthInNumCells[0]; xIndex <= _haloWidthInNumCells[0]; xIndex++)
      {
        if (xIndex != 0)
        {
          xDistanceSquare = pow((abs(xIndex) - 1) * _cellLength[0], 2);
        }
        else
        {
          xDistanceSquare = 0;
        }
        if (xDistanceSquare + yDistanceSquare + zDistanceSquare <= cutoffRadiusSquare)
        {
          long int offset = CellIndexOf3DIndex(xIndex, yIndex, zIndex);
          if (offset > 0)
          {
            forwardNeighbourOffsets.emplace_back(offset); // now vector
            ++forwardNeighbourIndex;
          }
          if (offset < 0)
          {
            backwardNeighbourOffsets.emplace_back(abs(offset)); // now vector
            ++backwardNeighbourIndex;
          }
        }
      }
    }
  }
  if (_haloWidthInNumCells[0] == 1 and _haloWidthInNumCells[1] == 1 and
      _haloWidthInNumCells[2] == 1)
  {
    rbmd_assert(forwardNeighbourIndex == maxNNeighbours);
    rbmd_assert(backwardNeighbourIndex == maxNNeighbours);
  }
  else
  {
    rbmd_assert(forwardNeighbourIndex <= maxNNeighbours);
    rbmd_assert(backwardNeighbourIndex <= maxNNeighbours);
  }
}

std::array<std::pair<unsigned long, unsigned long>, 14> LinkedCell::calculateCellPairOffsets() const
{
  long int o = CellIndexOf3DIndex(0, 0, 0); // origin
  long int x = CellIndexOf3DIndex(1, 0, 0); // displacement to the right
  long int y = CellIndexOf3DIndex(0, 1, 0); // displacement ...
  long int z = CellIndexOf3DIndex(0, 0, 1);
  long int xy = CellIndexOf3DIndex(1, 1, 0);
  long int yz = CellIndexOf3DIndex(0, 1, 1);
  long int xz = CellIndexOf3DIndex(1, 0, 1);
  long int xyz = CellIndexOf3DIndex(1, 1, 1);

  // minimize number of cells simultaneously in memory:
  std::array<std::pair<unsigned long, unsigned long>, 14> cellPairOffsets;

  cellPairOffsets[0] = std::make_pair(o, xyz);
  // evict xyz

  cellPairOffsets[1] = std::make_pair(o, yz);
  cellPairOffsets[2] = std::make_pair(x, yz);
  // evict yz

  cellPairOffsets[3] = std::make_pair(o, x);

  cellPairOffsets[4] = std::make_pair(o, xy);
  cellPairOffsets[5] = std::make_pair(xy, z);
  // evict xy

  cellPairOffsets[6] = std::make_pair(o, z);
  cellPairOffsets[7] = std::make_pair(x, z);
  cellPairOffsets[8] = std::make_pair(y, z);
  // evict z

  cellPairOffsets[9] = std::make_pair(o, y);
  cellPairOffsets[10] = std::make_pair(x, y);
  // evict x

  cellPairOffsets[11] = std::make_pair(o, xz);
  cellPairOffsets[12] = std::make_pair(y, xz);
  // evict xz

  cellPairOffsets[13] = std::make_pair(o, o);

  return cellPairOffsets;
}

constexpr int getChunkSize(size_t loop_size, size_t max_num_chunks, size_t max_chunk_size)
{
  return static_cast<int>(
    std::max(std::min(loop_size / max_num_chunks, max_chunk_size), static_cast<size_t>(1ul)));
}

void LinkedCell::UpdateViaCopies()
{
  const std::vector<Cell>::size_type numCells = _cells.size();
  std::vector<long> forwardNeighbourOffsets;  // now vector
  std::vector<long> backwardNeighbourOffsets; // now vector
  CalculateNeighbourIndices(forwardNeighbourOffsets, backwardNeighbourOffsets);

  // magic numbers: empirically determined to be somewhat efficient.
  const int chunk_size = getChunkSize(_cells.size(), 10000, 100);
#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(_OPENMP)
#pragma omp for schedule(dynamic, chunk_size)
#endif
    //if (update_step > 42)
    //{
    //  RBMDParallelUtil::printAtomInfo(51, this, "PreUpdateLeavingMolecules之前");
    //}
    for (std::vector<Cell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++)
    {
      _cells[cellIndex].PreUpdateLeavingMolecules();   /// 删除了51？
    }
    //if (update_step > 42)
    //{
    //  RBMDParallelUtil::printAtomInfo(51, this, "PreUpdateLeavingMolecules之后");
    //}

#if defined(_OPENMP)
#pragma omp for schedule(dynamic, chunk_size)
#endif
    //if (update_step > 42)
    //{
    //  RBMDParallelUtil::printAtomInfo(51, this, "UpdateLeavingMolecules之前");
    //}
    for (std::vector<Cell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++)
    {
      Cell& cell = _cells[cellIndex];

      for (unsigned long j = 0; j < backwardNeighbourOffsets.size(); j++)
      {
        const unsigned long neighbourIndex =
          cellIndex - backwardNeighbourOffsets.at(j); // now vector
        if (neighbourIndex >= _cells.size())
        {
          // handles cell_index < 0 (indices are unsigned!)
          rbmd_assert(cell.IsHaloCell());
          continue;
        }
        cell.UpdateLeavingMolecules(_cells[neighbourIndex]);
      }

      for (unsigned long j = 0; j < forwardNeighbourOffsets.size(); j++)
      {
        const unsigned long neighbourIndex =
          cellIndex + forwardNeighbourOffsets.at(j); // now vector
        if (neighbourIndex >= numCells)
        {
          rbmd_assert(cell.IsHaloCell());
          continue;
        }
        cell.UpdateLeavingMolecules(_cells[neighbourIndex]);
      }
    }
    //if (update_step > 42)
    //{
    //  RBMDParallelUtil::printAtomInfo(51, this, "UpdateLeavingMolecules之前");
    //}

#if defined(_OPENMP)
#pragma omp for schedule(dynamic, chunk_size)
#endif
    //if (update_step > 42)
    //{
    //  RBMDParallelUtil::printAtomInfo(51, this, "PostUpdateLeavingMolecules之前");
    //}
    for (std::vector<Cell>::size_type cellIndex = 0; cellIndex < _cells.size(); cellIndex++)
    {
      _cells[cellIndex].PostUpdateLeavingMolecules();
    }
  } // end pragma omp parallel
  //if (update_step > 42)
  //{
  //  RBMDParallelUtil::printAtomInfo(51, this, "PostUpdateLeavingMolecules之后");
  //}
}

bool LinkedCell::rebuild(double bBoxMin[3], double bBoxMax[3]) {
    // Log::global_log->info() << "REBUILD OF LinkedCells" << std::endl;
    for (int i = 0; i < 3; i++) {
        this->_boundingBoxMin[i] = bBoxMin[i];
        this->_boundingBoxMax[i] = bBoxMax[i];
        //		_haloWidthInNumCells[i] = ::ceil(_cellsInCutoff);
        _haloWidthInNumCells[i] = _cellsInCutoff;
    }
    // Log::global_log->info() << "Bounding box: " << "[" << bBoxMin[0] << ", " << bBoxMax[0] << "]" << " x " << "["
    // 		<< bBoxMin[1] << ", " << bBoxMax[1] << "]" << " x " << "[" << bBoxMin[2] << ", " << bBoxMax[2] << "]"
    // 		<< std::endl;

    int numberOfCells = 1;

    // Log::global_log->info() << "Using " << _cellsInCutoff << " cells in cutoff." << std::endl;
    float rc = (_cutoffRadius / _cellsInCutoff);

    for (int dim = 0; dim < 3; dim++) {
        _boxWidthInNumCells[dim] = floor((_boundingBoxMax[dim] - _boundingBoxMin[dim]) / rc);

        _cellsPerDimension[dim] = _boxWidthInNumCells[dim] + 2 * _haloWidthInNumCells[dim];

        // in each dimension at least one layer of (inner+boundary) cells necessary
        if (_cellsPerDimension[dim] == 2 * _haloWidthInNumCells[dim]) {

            // Log::global_log->error_always_output() << "LinkedCells::rebuild: region too small" << std::endl;
            // Simulation::exit(1);
            throw std::runtime_error("LinkedCells::rebuild: region too small");
        }

        numberOfCells *= _cellsPerDimension[dim];    //包含halo

        double diff = _boundingBoxMax[dim] - _boundingBoxMin[dim];
        _cellLength[dim] = diff / _boxWidthInNumCells[dim];
        _cellLengthReciprocal[dim] = _boxWidthInNumCells[dim] / diff;

        _haloLength[dim] = _haloWidthInNumCells[dim] * _cellLength[dim];

        _haloBoundingBoxMin[dim] = _boundingBoxMin[dim] - _haloLength[dim];
        _haloBoundingBoxMax[dim] = _boundingBoxMax[dim] + _haloLength[dim];
    }

    // Log::global_log->info() << "Cells per dimension (incl. halo): " << _cellsPerDimension[0] << " x "
    // 		<< _cellsPerDimension[1] << " x " << _cellsPerDimension[2] << std::endl;


    _cells.resize(numberOfCells);
    bool sendParticlesTogether = true;
    // If the width of the inner region is less than the width of the halo region
    // leaving particles and halo copy must be sent separately.
    if (_boxWidthInNumCells[0] < 2 * _haloWidthInNumCells[0]
        || _boxWidthInNumCells[1] < 2 * _haloWidthInNumCells[1]
        || _boxWidthInNumCells[2] < 2 * _haloWidthInNumCells[2]) {
        sendParticlesTogether = false;
        // throw std::runtime_error("LinkedCells (constructor): bounding box too small for calculated cell length");
    }


    InitializeCells();

    // TODO: We loose particles here as they are not communicated to the new owner
    // delete all Particles which are outside of the halo region
    deleteParticlesOutsideBox(_haloBoundingBoxMin, _haloBoundingBoxMax);   // TODO：gpu版本应该还是需要----用于动态负载均衡的
    // InitializeTraversal();  // TODO 3.21 暂时不做这个

    _cellsValid = false;

    return sendParticlesTogether;
}

void LinkedCell::UpdateViaCopiesLongRange() {
    const std::vector<Cell>::size_type numCells = _cells.size();

    for (std::vector<Cell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++) {
        _cells[cellIndex].PreUpdateLeavingMolecules(); // 不在当前cell的放到离开粒子
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic, chunk_size)
#endif
    for (std::vector<Cell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++) {
        Cell &cell = _cells[cellIndex];
        for (std::vector<Cell>::size_type neighbor_cell = 0; neighbor_cell < numCells; neighbor_cell++) {
            if (neighbor_cell >= _cells.size()) continue;
            cell.UpdateLeavingMolecules(_cells[neighbor_cell]);
        }
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic, chunk_size)
#endif
    for (std::vector<Cell>::size_type cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
        _cells[cellIndex].PostUpdateLeavingMolecules();
    }
} // end pragma omp parallel


void LinkedCell::update()
{
  //  update_step++;
  //  if (update_step > 42)
  //{
  //   RBMDParallelUtil::printAtomInfo(51, this, "进入Update时");
  //   RBMDParallelUtil::printAtomInfo(61, this, "进入Update时");
  //   RBMDParallelUtil::printAtomInfo(10, this, "进入Update时");
  //
  //}
#ifndef NDEBUG
  CheckAtomsInBox();
#endif

  // TODO: replace via a cellProcessor and a traverseCells call ?         -----------  UpdateViaTraversal 好像没用？！！！！！！！
#ifndef ENABLE_REDUCED_MEMORY_MODE
  //if (update_step > 42)
  //{
  //    RBMDParallelUtil::printAtomInfo(51, this, "UpdateViaCopies之前");
  //    //RBMDParallelUtil::printAtomInfo(61, this, "UpdateViaCopies之前");
  //    //RBMDParallelUtil::printAtomInfo(10, this, "UpdateViaCopies之前");
  //
  //}
  UpdateViaCopies(); // 普通模式只调用里这个
  //if (update_step > 42)
  //{
  //    RBMDParallelUtil::printAtomInfo(51, this, "UpdateViaCopies之后");
  //    //RBMDParallelUtil::printAtomInfo(944, this, "UpdateViaCopies之后");
  //}
#else
  //	update_via_coloring();
  std::array<long unsigned, 3> dims = { static_cast<long unsigned>(_cellsPerDimension[0]),
                                        static_cast<long unsigned>(_cellsPerDimension[1]),
                                        static_cast<long unsigned>(_cellsPerDimension[2]) };
  UpdateViaTraversal(); //没有调用这个
#endif
  _cellsValid = true;


#ifndef NDEBUG
  unsigned numBadMolecules = 0;

#if defined(_OPENMP)
#pragma omp parallel reduction(+ : numBadMolecules)
#endif
  {
    for (LinkedCellIterator tM = iterator(LinkedCellIterator::ALL_CELLS); tM.isValid(); ++tM)
    {
      if (not _cells[tM.getCellIndex()].InBox(*tM))
      {
        numBadMolecules++;
        // Log::global_log->error_always_output() << "particle " << tM->getID() << " in cell " << tM.getCellIndex()
        // 		<< ", which is" << (_cells[tM.getCellIndex()].isBoundaryCell() ? "" : " NOT")
        // 		<< " a boundarycell is outside of its cell after LinkedCells::update()." << std::endl;
        // Log::global_log->error_always_output() << "particle at (" << tM->r(0) << ", " << tM->r(1) << ", " << tM->r(2) << ")"
        // 		<< std::endl << "cell: [" << _cells[tM.getCellIndex()].getBoxMin(0) << ", "
        // 		<< _cells[tM.getCellIndex()].getBoxMax(0) << "] x [" << _cells[tM.getCellIndex()].getBoxMin(1)
        // 		<< ", " << _cells[tM.getCellIndex()].getBoxMax(1) << "] x ["
        // 		<< _cells[tM.getCellIndex()].getBoxMin(2) << ", " << _cells[tM.getCellIndex()].getBoxMax(2) << "]"
        // 		<< std::endl;
      }
    }
  }


  if (numBadMolecules > 0)
  {
    // Log::global_log->error() << "Found " << numBadMolecules << " outside of their correct cells. Aborting." << std::endl;
    // Simulation::exit(311);
    throw std::runtime_error("Found  numBadMolecules in update");
  }
#endif
}



LinkedCell::~LinkedCell() {}

LinkedCell::LinkedCell() {}

bool LinkedCell::IsInBoundingBox(double* r) const
{
  if (r[0] >= _boundingBoxMin[0] && r[1] >= _boundingBoxMin[1] && r[2] >= _boundingBoxMin[2] &&
      r[0] < _boundingBoxMax[0] && r[1] < _boundingBoxMax[1] && r[2] < _boundingBoxMax[2])
  {
    return true;
  }
  else
  {
    return false;
  }
}

std::vector<Cell>& LinkedCell::GetCells()
{
  return this->_cells;
}

int LinkedCell::GetCellPerDim(const int d) const
{
  return _cellsPerDimension[d];
}
