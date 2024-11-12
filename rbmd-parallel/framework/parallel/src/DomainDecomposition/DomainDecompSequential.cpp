#include <cmath>
#include "DomainDecompSequential.h"
#include "FullShell.hpp"
#include "DomainServiceLocator.h"
#include "Forcehelper.hpp"
#include "CommunicationPartner.h"
DomainDecompSequential::~DomainDecompSequential() {}

// void DomainDecompSequential::addLeavingMolecules(std::vector<Atom> &invalidMolecules, LinkedCell *moleculeContainer)  {
//     for (auto& molecule : invalidMolecules) {
//         for (auto dim : {0, 1, 2}) {
//             auto shift = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);
//             auto r = molecule.r(dim);
//             if (r < moleculeContainer->getBoundingBoxMin(dim)) {
//                 r = r + shift;
//                 if (r >= moleculeContainer->getBoundingBoxMax(dim)) {
//                     r = std::nextafter(moleculeContainer->getBoundingBoxMax(dim), -1);
//                 }
//             } else if (r >= moleculeContainer->getBoundingBoxMax(dim)) {
//                 r = r - shift;
//                 if (r < moleculeContainer->getBoundingBoxMin(dim)) {
//                     r = moleculeContainer->getBoundingBoxMin(dim);
//                 }
//             }
//             molecule.Setr(dim, r);
//         }
//     }
//     moleculeContainer->addParticles(invalidMolecules);
//     invalidMolecules.clear();
// }

double DomainDecompSequential::GetBoundingBoxMin(int dimension) const{
    return 0;
}

double DomainDecompSequential::GetBoundingBoxMax(const int dimension) const {
    return  DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(dimension);  // 顺序版直接返回全局长度就行了
}



std::vector<CommunicationPartner> DomainDecompSequential::GetNeighboursFromHaloRegion(const HaloRegion &haloRegion,
    double cutoff) {
    std::cout << "单个进程不需要此方法，后续优化........." << std::endl;   //TODO

    throw std::runtime_error("单个进程不需要此方法");
}


void DomainDecompSequential::exchangeMolecules(LinkedCell *moleculeContainer) {
    if (moleculeContainer->IsInvalidParticleReturner()) {
        // autopas mode!
        //Log::global_log->debug() << "DDBase: Adding + shifting invalid particles." << std::endl;
        // in case the molecule container returns invalid particles using GetInvalidParticlesRef(), we have to handle them directly.
        addLeavingMolecules(moleculeContainer->GetInvalidParticlesRef(), moleculeContainer);
        // now use direct scheme to transfer the rest!
        FullShell fs;
        double rmin[3];  // lower corner
        double rmax[3];  // higher corner
        for (int d = 0; d < 3; d++) {
            rmin[d] = GetBoundingBoxMin(d);
            rmax[d] = GetBoundingBoxMax(d);
        }
        HaloRegion ownRegion = {rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0, 0.};
        bool coversWholeDomain[3];
        double cellLengthDummy[3]{};
        //Log::global_log->debug() << "DDBase: Populating halo." << std::endl;
        auto haloExportRegions =
                fs.getHaloExportForceImportRegions(ownRegion, moleculeContainer->GetCutoff(),
                                                   coversWholeDomain, cellLengthDummy);
        for (auto haloExportRegion : haloExportRegions) {
            PopulateHaloLayerWithCopiesDirect(haloExportRegion, moleculeContainer,
                                              true /*positionCheck, same as doLeavingExchange*/);
        }
    } else {
        // default ls1-mode (non-autopas, so linked-cells!)
        for (unsigned d = 0; d < 3; ++d) {
            HandleDomainLeavingParticles(d, moleculeContainer);
        }
        for (unsigned d = 0; d < 3; ++d) {
            PopulateHaloLayerWithCopies(d, moleculeContainer);
        }
    }
}

void DomainDecompSequential::ExchangeForces(LinkedCell *linkedCell) {
    for (unsigned d = 0; d < 3; ++d) {
        handleForceExchange(d,linkedCell);
    }
}



int DomainDecompSequential::GetNonBlockingStageCount() {
    return -1;
}

bool DomainDecompSequential::QueryBalanceAndExchangeNonBlocking(bool force_re_balancing, LinkedCell *linked_Cell,
                                        double etime) {
    return false;
}

void DomainDecompSequential::BalanceAndExchange(double last_traversal_time, bool force_re_balancing,
    LinkedCell *linked_Cell) {
    std::cout << "进入单线程的通信逻辑!" << std::endl;
    exchangeMolecules(linked_Cell);
}

bool DomainDecompSequential::IsPositionInThisDomain(const double x, const double y, const double z) const {
    return !(x < GetBoundingBoxMin(0) || x >= GetBoundingBoxMax(0) ||
         y < GetBoundingBoxMin(1) || y >= GetBoundingBoxMax(1) ||
         z < GetBoundingBoxMin(2) || z >= GetBoundingBoxMax(2));
}

void DomainDecompSequential::GetBoundingBoxMinMax(double *min, double *max) {
    for(int d = 0; d < 3; d++) {
        min[d] = GetBoundingBoxMin(d);
        max[d] = GetBoundingBoxMax(d);
    }
}

int DomainDecompSequential::GetCurrentRank() const {
    return _current_rank;
}

int DomainDecompSequential::GetTotalRank() const {
    return _total_ranks;
}

void DomainDecompSequential::CommunicatorInitializeMemory(const int numValues, int key) {
    _coll_comm_seq.InitializeMemory(numValues);
}

void DomainDecompSequential::CommunicatorFinalize() {
    _coll_comm_seq.Finalize();
}

void DomainDecompSequential::CommunicatorAppendInt(int intValue) {
    _coll_comm_seq.AppendInt(intValue);
}

void DomainDecompSequential::CommunicatorAppendUnsignedLong(unsigned long unsigned_long_value) {
    _coll_comm_seq.AppendUnsignedLong(unsigned_long_value);
}

void DomainDecompSequential::CommunicatorAppendDouble(double doubleValue) {
    _coll_comm_seq.AppendDouble(doubleValue);
}

void DomainDecompSequential::CommunicatorAppendLongDouble(long double longDoubleValue) {
    _coll_comm_seq.AppendLongDouble(longDoubleValue);
}

void DomainDecompSequential::CommunicatorAppendFloat(float floatValue) {
    _coll_comm_seq.AppendFloat(floatValue);
}



int DomainDecompSequential::CommunicatorGetInt() {
    return _coll_comm_seq.GetInt();
}

unsigned long DomainDecompSequential::CommunicatorGetUnsignedLong() {
    return _coll_comm_seq.GetUnsignedLong();
}

float DomainDecompSequential::CommunicatorGetFloat() {
    return _coll_comm_seq.GetFloat();
}

double DomainDecompSequential::CommunicatorGetDouble() {
    return _coll_comm_seq.GetDouble();
}

long double DomainDecompSequential::CommunicatorGetLongDouble() {
    return _coll_comm_seq.GetLongDouble();
}

void DomainDecompSequential::CommunicatorAllreduceSum() {
    _coll_comm_seq.AllReduceSum();
}

void DomainDecompSequential::CommunicatorAllreduceCustom(ReduceType type) {
    _coll_comm_seq.AllReduceCustom(type);
}

void DomainDecompSequential::CommunicatorScanSum() {
    _coll_comm_seq.ScanSum();
}

void DomainDecompSequential::CommunicatorBroadcast(int root) {
    _coll_comm_seq.Broadcast(root);
}

std::vector<int> DomainDecompSequential::GetNeighbourRanks() {
    return std::vector<int>(0);
}

std::vector<std::vector<std::vector<int>>> DomainDecompSequential::GetAllRanks() {
    return std::vector<std::vector<std::vector<int>>>(0);
}

size_t DomainDecompSequential::GetTotalDataSize(unsigned localN, float* minrnd, float* maxrnd) {
    return _coll_comm_seq.GetTotalDataSize();
}

void DomainDecompSequential::SetSendLeavingWithCopiesStrategy(bool send_together) {
    CommunicatorInitializeMemory(1,0);
    CommunicatorAppendInt(!send_together);
    CommunicatorAllreduceSum();
    _sendLeavingAndCopiesSeparately = CommunicatorGetInt();
    CommunicatorFinalize();
}

bool DomainDecompSequential::IsSendLeavingWithCopies() {
    return _sendLeavingAndCopiesSeparately == 0;
}

void DomainDecompSequential::CommunicatorInit(MPI_Comm mpi_communicator, std::size_t capacity, int key) {
    throw std::runtime_error("Seq 考虑放弃支持");
}


#ifdef ENABLE_MPI
MPI_Comm DomainDecompSequential::GetMPIComm(){
    return MPI_COMM_WORLD;
}
#endif


















