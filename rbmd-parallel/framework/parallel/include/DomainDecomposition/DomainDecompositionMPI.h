#pragma once
#include "DomainDecompBase.h"
#include <mpi.h>
#include "CommunicationPartner.h"
#include "BaseNeighbourCommunicationScheme.h"
#include "CommunicatorMPIInterface.h"
#include "IndirectNeighbourCommunicationScheme.h"
#include "CommunicatorNonBlockingScheduler.h"

class DomainDecompositionMPI : DomainDecompBase {
private:
    // mpi通信域
    MPI_Comm _mpi_comm;
    MPI_Datatype _mpiParticleType;
    MPI_Datatype _mpiParticleForceType;
    std::array<int, DIM> _gridSize{}; //!< MPI进程网格的每个维度中的进程数
    int _coords[DIM]{}; //!< MPI进程网格中进程的坐标
    std::unique_ptr<CommunicatorNonBlockingScheduler> _coll_comm_mpi; //TODO 是否正确？
    std::unique_ptr<IndirectNeighbourCommunicationScheme> _neighbourCommunicationScheme;   //TODO 用于交换邻居的！！！！1

protected:
    void InitMPIGridDims();


    /**
     * Prepares the stageNumber'th stage.
     * This includes sending particles, that have to be send in that stage.
     * @param moleculeContainer pointer to the molecule container
     * @param domain pointer to the domain
     * @param stageNumber the number of the stage, the communication is in.
     */
    virtual void prepareNonBlockingStageImpl(LinkedCell *moleculeContainer,
                                             unsigned int stageNumber, MessageType msgType,
                                             bool removeRecvDuplicates = false);

    /**
     * Finishes the stageNumber'th stage.
     * This includes receiving the particles, that have to be received in that stage.
     * @param forceRebalancing true if rebalancing should be forced
     * @param moleculeContainer pointer to the molecule container
     * @param domain pointer to the domain
     * @param stageNumber the number of the stage, the communication is in.
     */
    virtual void finishNonBlockingStageImpl(LinkedCell *moleculeContainer,
                                            unsigned int stageNumber, MessageType msgType,
                                            bool removeRecvDuplicates = false);
    MPI_Datatype getMPIParticleType() {
        return _mpiParticleType;
    }
    MPI_Datatype getMPIParticleForceType() {
        return _mpiParticleForceType;
    }

    void ExchangeMoleculesMPI(LinkedCell* moleculeContainer, MessageType msgType,
                              bool doHaloPositionCheck = true, bool removeRecvDuplicates = false);

public:
    DomainDecompositionMPI() = delete;
    DomainDecompositionMPI(int current_rank,int total_rank);
    DomainDecompositionMPI(MPI_Comm mpi_comm, const std::array<int, DIM> &gridSize,int current_rank,int total_rank);
    ~DomainDecompositionMPI() override;

#pragma region 通信方法
    void CommunicatorInit(MPI_Comm mpi_communicator, std::size_t capacity, int key = 0) override;

    //! has to call init method of a CollComm class
    void CommunicatorInitializeMemory(int numValues, int key ) override;

    //! has to call finalize method of a CollComm class
    void CommunicatorFinalize() override;

    //! has to call appendInt method of a CollComm class
    void CommunicatorAppendInt(int intValue) override;

    //! has to call appendUnsLong method of a CollComm class
    void CommunicatorAppendUnsignedLong(unsigned long unsigned_long_value) override;

    //! has to call appendFloat method of a CollComm class
    void CommunicatorAppendFloat(float floatValue) override;

    //! has to call appendDouble method of a CollComm class
    void CommunicatorAppendDouble(double doubleValue) override;

    //! has to call appendLongDouble method of a CollComm class
    void CommunicatorAppendLongDouble(long double longDoubleValue) override;

    //! has to call getInt method of a CollComm class
    int CommunicatorGetInt() override;

    //! has to call getUnsLong method of a CollComm class
    unsigned long CommunicatorGetUnsignedLong() override;

    //! has to call getFloat method of a CollComm class
    float CommunicatorGetFloat() override;

    //! has to call getDouble method of a CollComm class
    double CommunicatorGetDouble() override;

    //! has to call getLongDouble method of a CollComm class
    long double CommunicatorGetLongDouble() override;

    //! has to call allreduceSum method of a CollComm class (none in sequential version)
    void CommunicatorAllreduceSum() override;

    //! has to call allreduceSum method of a CollComm class (none in sequential version), allows for values of previous iteration.
    void CommunicatorAllreduceSumAllowPrevious() override;

    //! has to call allreduceCustom method of a CollComm class (none in sequential version)
    void CommunicatorAllreduceCustom(ReduceType type) override;

    //! has to call scanSum method of a CollComm class (none in sequential version)
    void CommunicatorScanSum() override;

    //! has to call broadcast method of a CollComm class (none in sequential version)
    void CommunicatorBroadcast(int root = 0) override;

    void SetSendLeavingWithCopiesStrategy(bool send_together)  override;

    bool IsSendLeavingWithCopies() override;

    std::vector<CommunicationPartner> GetNeighboursFromHaloRegion( const HaloRegion& haloRegion, double cutoff) override;   // TODO 要不要放别的地方去

    void GetBoundingBoxMinMax(double *min, double *max) override;

    void ExchangeForces(LinkedCell* linkedCell) override;

    bool IsPositionInThisDomain(double x, double y, double z) const override;

    int GetCurrentRank() const override;

    int GetTotalRank() const override;

    int GetNonBlockingStageCount() override;

    size_t GetTotalDataSize(unsigned localN, float* minrnd, float* maxrnd) override;

    void  SetCommunicationScheme(const std::string& comm_scheme, const std::string& zonal_method);

    MPI_Comm& GetMPIComm() override { return _mpi_comm; }
#pragma endregion

    // documentation in base class
    void PrepareNonBlockingStage(bool forceRebalancing,
                                 LinkedCell *moleculeContainer, unsigned int stageNumber);

    // documentation in base class
    void FinishNonBlockingStage(bool forceRebalancing,
                                LinkedCell *moleculeContainer, unsigned int stageNumber);

    bool QueryBalanceAndExchangeNonBlocking(bool force_re_balancing, LinkedCell *linked_Cell, double etime) override;

    void BalanceAndExchange(double last_traversal_time, bool force_re_balancing, LinkedCell *linked_Cell) override;

    void BalanceAndExchangeInitNonBlocking(bool forceRebalancing, LinkedCell *moleculeContainer);

    double GetBoundingBoxMin(int dimension) const override;

    double GetBoundingBoxMax(int dimension) const override;

    std::vector<int> GetNeighbourRanks() override;

    std::vector<std::vector<std::vector<int>>> GetAllRanks() override;

    void InitCommunicationPartners(double cutoffRadius,  LinkedCell* moleculeContainer);

};
