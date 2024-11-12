#include "DomainDecompositionMPI.h"

#include <RBMDParallelUtil.h>

#include <csignal>
#include <iostream>

#include "CommunicatorNonBlockingScheduler.h"
#include "DomainServiceLocator.h"
#include "FullShell.hpp"
#include "IndirectNeighbourCommunicationScheme.h"
#include "MPIChecker.hpp"

void DomainDecompositionMPI::InitMPIGridDims() {
  rbmd_assert(DIM == 3);
  int period[DIM] = {1, 1, 1};  // 1(true) when using periodic boundary
                                // conditions in the corresponding dimension
  int reorder =
      1;  // 1(true) if the ranking may be reordered by MPI_Cart_create
  {
    auto numProcsGridSize =
        _gridSize[0] * _gridSize[1] * _gridSize[2];  // TODO What doing?
    if (numProcsGridSize != _total_ranks and numProcsGridSize != 0) {
      // Log::global_log->error() << "DomainDecomposition: Wrong grid size
      // given!" << std::endl; Log::global_log->error() << "\tnumProcs is " <<
      // _total_ranks << "," << std::endl; Log::global_log->error() << "\tbut
      // grid is " << _gridSize[0] << " x " << _gridSize[1] << " x " <<
      // _gridSize[2] << std::endl; Log::global_log->error() << "\tresulting in
      // " << numProcsGridSize << " subdomains!" << std::endl;
      // Log::global_log->error() << "\tplease check your input file!" <<
      // std::endl; Simulation::exit(2134);
      throw std::runtime_error("error mpi grdims!");
    }
  }

  MPI_CHECK(MPI_Dims_create(_total_ranks, DIM, _gridSize.data()));
  MPI_CHECK(MPI_Cart_create(_mpi_comm, DIM, _gridSize.data(), period, reorder,
                            &_mpi_comm));

  std::cout << "MPI grid dimensions: " << _gridSize[0] << ", " << _gridSize[1]
            << ", " << _gridSize[2] << std::endl;
  MPI_CHECK(MPI_Comm_rank(_mpi_comm, &_current_rank));
  MPI_CHECK(MPI_Cart_coords(_mpi_comm, _current_rank, DIM, _coords));
  std::cout << "MPI coordinate of current process: " << _coords[0] << ", "
            << _coords[1] << ", " << _coords[2] << std::endl;
}

void DomainDecompositionMPI::prepareNonBlockingStageImpl(
    LinkedCell *moleculeContainer, unsigned int stageNumber,
    MessageType msgType, bool removeRecvDuplicates) {
  rbmd_assert(stageNumber < _neighbourCommunicationScheme->GetCommDims());
  _neighbourCommunicationScheme->prepareNonBlockingStageImpl(
      moleculeContainer, stageNumber, msgType, removeRecvDuplicates, this);
}

void DomainDecompositionMPI::finishNonBlockingStageImpl(
    LinkedCell *moleculeContainer, unsigned int stageNumber,
    MessageType msgType, bool removeRecvDuplicates) {
  _neighbourCommunicationScheme->finishNonBlockingStageImpl(
      moleculeContainer, stageNumber, msgType, removeRecvDuplicates, this);
}

std::vector<CommunicationPartner>
DomainDecompositionMPI::GetNeighboursFromHaloRegion(
    const HaloRegion &haloRegion, double cutoff) {
  // TODO: change this method for support of midpoint rule, half shell, eighth
  // shell, Neutral Territory
  //  currently only one process per region is possible.
  int rank;
  int regionCoords[DIM];
  for (unsigned int d = 0; d < DIM; d++) {
    regionCoords[d] = _coords[d] + haloRegion.offset[d];
  }
  // TODO: only full shell! (otherwise more neighbours possible)
  MPI_CHECK(
      MPI_Cart_rank(GetMPIComm(), regionCoords,
                    &rank));  // does automatic shift for periodic boundaries
  double haloLow[3];
  double haloHigh[3];
  double boundaryLow[3];
  double boundaryHigh[3];
  double shift[3];
  bool enlarged[3][2];

  for (unsigned int d = 0; d < DIM; d++) {
    haloLow[d] = haloRegion.rmin[d];
    haloHigh[d] = haloRegion.rmax[d];
    // TODO: ONLY FULL SHELL!!!

    boundaryLow[d] =
        haloRegion.rmin[d] -
        haloRegion.offset[d] * haloRegion.width;  // rmin[d] if offset[d]==0
    boundaryHigh[d] =
        haloRegion.rmax[d] -
        haloRegion.offset[d] *
            haloRegion.width;  // if offset[d]!=0 : shift by cutoff in negative
                               // offset direction
    if (_coords[d] == 0 and haloRegion.offset[d] == -1) {
      shift[d] =
          DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(d);
    } else if (_coords[d] == _gridSize[d] - 1 and haloRegion.offset[d] == 1) {
      shift[d] =
          -DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(d);
    } else {
      shift[d] = 0.;
    }
    enlarged[d][0] = false;
    enlarged[d][1] = false;
  }
  // initialize using initializer list - here a vector with one element is
  // created
  std::vector<CommunicationPartner> temp;
  temp.emplace_back(rank, haloLow, haloHigh, boundaryLow, boundaryHigh, shift,
                    haloRegion.offset, enlarged);
  return temp;
}

void DomainDecompositionMPI::ExchangeMoleculesMPI(LinkedCell *moleculeContainer,
                                                  MessageType msgType,
                                                  bool doHaloPositionCheck,
                                                  bool removeRecvDuplicates) {
  _neighbourCommunicationScheme->exchangeMoleculesMPI(
      moleculeContainer, msgType, removeRecvDuplicates, this,
      doHaloPositionCheck);
}

DomainDecompositionMPI::DomainDecompositionMPI(int current_rank, int total_rank)
    : _mpiParticleType(NULL), _mpiParticleForceType(NULL) {
  _mpi_comm = MPI_COMM_WORLD;
  _coll_comm_mpi = std::make_unique<CommunicatorNonBlockingScheduler>();
  _neighbourCommunicationScheme =
      std::make_unique<IndirectNeighbourCommunicationScheme>(new FullShell());
  _gridSize = {0, 0, 0};
  this->_current_rank = current_rank;
  this->_total_ranks = total_rank;
  InitMPIGridDims();
  this->CommunicatorInit(_mpi_comm, 1024);  // TODO 2024年4月2日 暂时这样
}

DomainDecompositionMPI::DomainDecompositionMPI(
    MPI_Comm mpi_comm, const std::array<int, 3> &gridSize, int current_rank,
    int total_rank)
    : _mpiParticleType(NULL), _mpiParticleForceType(NULL) {
  _mpi_comm = mpi_comm;
  _gridSize = gridSize;
  _coll_comm_mpi = std::make_unique<CommunicatorNonBlockingScheduler>();
  _neighbourCommunicationScheme =
      std::make_unique<IndirectNeighbourCommunicationScheme>(new FullShell());
  this->_current_rank = current_rank;
  this->_total_ranks = total_rank;
  InitMPIGridDims();
  this->CommunicatorInit(_mpi_comm, 1024);
}

DomainDecompositionMPI::~DomainDecompositionMPI() {
  free(_mpiParticleType);
  free(_mpiParticleForceType);
  MPI_Comm_free(&_mpi_comm);
  //    _gridSize.;  //TODO
}

void DomainDecompositionMPI::CommunicatorInitializeMemory(const int numValues,
                                                          int key) {
  _coll_comm_mpi->InitializeMemory(numValues);
}

void DomainDecompositionMPI::CommunicatorFinalize() {
  _coll_comm_mpi->Finalize();
}

void DomainDecompositionMPI::CommunicatorAppendInt(int intValue) {
  _coll_comm_mpi->AppendInt(intValue);
}

void DomainDecompositionMPI::CommunicatorAppendUnsignedLong(
    unsigned long unsigned_long_value) {
  _coll_comm_mpi->AppendUnsignedLong(unsigned_long_value);
}

void DomainDecompositionMPI::CommunicatorAppendDouble(double doubleValue) {
  _coll_comm_mpi->AppendDouble(doubleValue);
}

void DomainDecompositionMPI::CommunicatorAppendLongDouble(
    long double longDoubleValue) {
  _coll_comm_mpi->AppendLongDouble(longDoubleValue);
}

void DomainDecompositionMPI::CommunicatorAppendFloat(float floatValue) {
  _coll_comm_mpi->AppendFloat(floatValue);
}

int DomainDecompositionMPI::CommunicatorGetInt() {
  return _coll_comm_mpi->GetInt();
}

unsigned long DomainDecompositionMPI::CommunicatorGetUnsignedLong() {
  return _coll_comm_mpi->GetUnsignedLong();
}

float DomainDecompositionMPI::CommunicatorGetFloat() {
  return _coll_comm_mpi->GetFloat();
}

double DomainDecompositionMPI::CommunicatorGetDouble() {
  return _coll_comm_mpi->GetDouble();
}

long double DomainDecompositionMPI::CommunicatorGetLongDouble() {
  return _coll_comm_mpi->GetLongDouble();
}

void DomainDecompositionMPI::CommunicatorAllreduceSum() {
  _coll_comm_mpi->AllReduceSum();
}

void DomainDecompositionMPI::CommunicatorAllreduceCustom(ReduceType type) {
  _coll_comm_mpi->AllReduceCustom(type);
}

void DomainDecompositionMPI::CommunicatorScanSum() {
  _coll_comm_mpi->ScanSum();
}

void DomainDecompositionMPI::CommunicatorBroadcast(int root) {
  _coll_comm_mpi->Broadcast(root);
}

void DomainDecompositionMPI::SetSendLeavingWithCopiesStrategy(
    bool send_together) {  // TODO 可能错误
  CommunicatorInitializeMemory(1, 0);
  CommunicatorAppendInt(!send_together);
  CommunicatorAllreduceSum();
  _sendLeavingAndCopiesSeparately = CommunicatorGetInt();
  CommunicatorFinalize();
}

bool DomainDecompositionMPI::IsSendLeavingWithCopies() {
  return _sendLeavingAndCopiesSeparately == 0;
}

void DomainDecompositionMPI::PrepareNonBlockingStage(
    bool forceRebalancing, LinkedCell *moleculeContainer,
    unsigned int stageNumber) {
  if (this->IsSendLeavingWithCopies()) {
    prepareNonBlockingStageImpl(moleculeContainer, stageNumber,
                                LEAVING_AND_HALO_COPIES);
  } else {
    // Would first need to send leaving, then halo -> not good for overlapping!
    // Log::global_log->error() << "nonblocking P2P using separate messages for
    // leaving and halo is currently not "
    //                        "supported. Please use the indirect neighbor
    //                        communication scheme!"
    //                     << std::endl;
    // Simulation::exit(235861);
    throw std::runtime_error("can not use");
  }
}

void DomainDecompositionMPI::FinishNonBlockingStage(
    bool forceRebalancing, LinkedCell *moleculeContainer,
    unsigned int stageNumber) {
  if (this->IsSendLeavingWithCopies()) {
    finishNonBlockingStageImpl(moleculeContainer, stageNumber,
                               LEAVING_AND_HALO_COPIES);
  } else {
    // Would first need to send leaving, then halo -> not good for overlapping!
    // Log::global_log->error()
    //     << "nonblocking P2P using separate messages for leaving and halo is
    //     currently not supported." << std::endl;
    // Simulation::exit(235861);
    throw std::runtime_error("can not use");
  }
}

bool DomainDecompositionMPI::QueryBalanceAndExchangeNonBlocking(
    bool force_re_balancing, LinkedCell *linked_Cell, double etime) {
  return true;
}

void DomainDecompositionMPI::BalanceAndExchange(double last_traversal_time,
                                                bool force_re_balancing,
                                                LinkedCell* linked_Cell)
{  
  //    MPI_Barrier(MPI_COMM_WORLD); // wait for all processes to reach this
  //    point (to avoid race conditions in output) sleep(15); // wait for all
  //    processes to reach this point (to avoid race conditions in output)
  if (this->IsSendLeavingWithCopies())
  {
    //if (step > 20)
    //{
    //  RBMDParallelUtil::printAtomInfo(500, linked_Cell, "ExchangeMoleculesMPI之前");
    //  RBMDParallelUtil::printAtomInfo(944, linked_Cell, "ExchangeMoleculesMPI之前");
    //}
    ExchangeMoleculesMPI(linked_Cell, LEAVING_AND_HALO_COPIES);
    //if (step > 20)
    //{
    //  RBMDParallelUtil::printAtomInfo(500, linked_Cell, "ExchangeMoleculesMPI之后");
    //  RBMDParallelUtil::printAtomInfo(944, linked_Cell, "ExchangeMoleculesMPI之后");
    //}
  }
  else
  {
    if (this->IsSendLeavingWithCopies())
    {
      //if (step > 20)
      //{
      //  RBMDParallelUtil::printAtomInfo(500, linked_Cell, "ExchangeMoleculesMPI之前");
      //  RBMDParallelUtil::printAtomInfo(944, linked_Cell, "ExchangeMoleculesMPI之前");
      //}
      // Log::global_log->debug() << "DD: Sending Leaving." << std::endl;
      ExchangeMoleculesMPI(linked_Cell, LEAVING_ONLY);
      linked_Cell->DeleteOuterParticles();
      // Log::global_log->debug() << "DD: Sending Halos." << std::endl;
      ExchangeMoleculesMPI(linked_Cell, HALO_COPIES);
      //if (step > 20)
      //{
      //  RBMDParallelUtil::printAtomInfo(500, linked_Cell, "ExchangeMoleculesMPI之后");
      //  RBMDParallelUtil::printAtomInfo(944, linked_Cell, "ExchangeMoleculesMPI之后");
      //}
    }
  }
}

void DomainDecompositionMPI::BalanceAndExchangeInitNonBlocking(
    bool forceRebalancing, LinkedCell *moleculeContainer) {
  // nothing to do
}

double DomainDecompositionMPI::GetBoundingBoxMin(int dimension) const {
  const auto global_length = DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(dimension);
  const auto offset = DomainServiceLocator::Instance().GetDomain()->GetOffset(dimension);
  return offset+_coords[dimension] * global_length / _gridSize[dimension];
}

double DomainDecompositionMPI::GetBoundingBoxMax(int dimension) const {
  const auto global_length = DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(dimension);
  const auto offset = DomainServiceLocator::Instance().GetDomain()->GetOffset(dimension);
  if(_coords[dimension] + 1 == _gridSize[dimension]){
    return offset+global_length;
  }
  return offset+(_coords[dimension] + 1) * global_length / _gridSize[dimension];
}

std::vector<int> DomainDecompositionMPI::GetNeighbourRanks() {
#if defined(ENABLE_MPI)
  std::vector<int> neighbours;
  if (_total_ranks == 1) {
    for (int i = 0; i < 6; i++) neighbours.push_back(_current_rank);
  } else {
    neighbours = _neighbourCommunicationScheme->get3StageNeighbourRanks();
  }
  return neighbours;
#else
  return std::vector<int>(0);
#endif
}

std::vector<std::vector<std::vector<int>>>
DomainDecompositionMPI::GetAllRanks() {
#ifdef ENABLE_MPI
  std::vector<std::vector<std::vector<int>>> ranks;
  int myRank;
  MPI_Comm_rank(_mpi_comm, &myRank);
  int numProcessors;
  MPI_Comm_size(_mpi_comm, &numProcessors);

  ranks.resize(_gridSize[0]);
  for (int i = 0; i < _gridSize[0]; i++) {
    ranks[i].resize(_gridSize[1]);
    for (int j = 0; j < _gridSize[1]; j++) {
      ranks[i][j].resize(_gridSize[2]);
    }
  }
  int coords[3];
  for (int i = 0; i < numProcessors; i++) {
    MPI_Cart_coords(_mpi_comm, i, 3, coords);
    ranks[coords[0]][coords[1]][coords[2]] = i;
  }
  return ranks;
#else
  return std::vector<std::vector<std::vector<int>>>(0);
#endif
}

void DomainDecompositionMPI::CommunicatorAllreduceSumAllowPrevious() {
  std::cout << "暂未实现" << std::endl;
}

void DomainDecompositionMPI::GetBoundingBoxMinMax(double *min, double *max) {
  for (int d = 0; d < 3; d++) {
    min[d] = GetBoundingBoxMin(d);
    max[d] = GetBoundingBoxMax(d);
  }
}

void DomainDecompositionMPI::ExchangeForces(LinkedCell *linkedCell) {
  _neighbourCommunicationScheme->exchangeMoleculesMPI(linkedCell, FORCES, false,
                                                      this);
}

bool DomainDecompositionMPI::IsPositionInThisDomain(double x, double y,
                                                    double z) const {
  return !(x < GetBoundingBoxMin(0) || x >= GetBoundingBoxMax(0) ||
           y < GetBoundingBoxMin(1) || y >= GetBoundingBoxMax(1) ||
           z < GetBoundingBoxMin(2) || z >= GetBoundingBoxMax(2));
}

int DomainDecompositionMPI::GetCurrentRank() const { return _current_rank; }

int DomainDecompositionMPI::GetTotalRank() const { return _total_ranks; }

int DomainDecompositionMPI::GetNonBlockingStageCount() {
  return _neighbourCommunicationScheme->GetCommDims();
}

size_t DomainDecompositionMPI::GetTotalDataSize(unsigned localN, float *minrnd,
                                                float *maxrnd) {
  std::vector<unsigned> moldistribution(_total_ranks);
  MPI_CHECK(MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution.data(), 1,
                          MPI_UNSIGNED, _mpi_comm));
  unsigned globalN = 0;
  for (int r = 0; r < _current_rank; r++) globalN += moldistribution[r];
  unsigned localNbottom = globalN;
  globalN += moldistribution[_current_rank];
  unsigned localNtop = globalN;
  for (int r = _current_rank + 1; r < _total_ranks; r++)
    globalN += moldistribution[r];
  *minrnd = (float)localNbottom / globalN;
  *maxrnd = (float)localNtop / globalN;
  return globalN;
}

void DomainDecompositionMPI::SetCommunicationScheme(
    const std::string &comm_scheme, const std::string &zonal_method) {
  ZonalMethod *zonalMethodP = nullptr;
  if (zonal_method == "fs") {
    zonalMethodP = new FullShell();
  } else {
    throw std::runtime_error("暂未实现的方法");
  }

  if (comm_scheme == "indirect") {
    _neighbourCommunicationScheme =
        std::make_unique<IndirectNeighbourCommunicationScheme>(zonalMethodP);
  }
}

void DomainDecompositionMPI::CommunicatorInit(MPI_Comm mpi_communicator,
                                              std::size_t capacity, int key) {
  _coll_comm_mpi->Init(mpi_communicator, capacity, key);
}

void DomainDecompositionMPI::InitCommunicationPartners(
    double cutoffRadius, LinkedCell *moleculeContainer) {
  for (int d = 0; d < DIM; ++d) {
    _neighbourCommunicationScheme->setCoverWholeDomain(d, _gridSize[d] == 1);
  }
  _neighbourCommunicationScheme->initCommunicationPartners(cutoffRadius, this,
                                                           moleculeContainer);
}
