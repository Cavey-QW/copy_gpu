#pragma onece

#include "CommunicatorSequential.h"
#include "DomainDecompBase.h"

// 序列版 顺序版 mpirun 1 或 不启用 是dd的一种特殊情况  负责操作domain？---这个是内聚吧
class DomainDecompSequential : DomainDecompBase {
private:
    CommunicatorSequential _coll_comm_seq;

protected:
    std::vector<CommunicationPartner> GetNeighboursFromHaloRegion(const HaloRegion &haloRegion, double cutoff) override;

public:
    DomainDecompSequential() = default;

    ~DomainDecompSequential() override;

    //void addLeavingMolecules(std::vector<Atom> &invalidMolecules, LinkedCell *moleculeContainer) override;

    double GetBoundingBoxMin(int dimension) const override;

    double GetBoundingBoxMax(int dimension) const override;

    void exchangeMolecules(LinkedCell *moleculeContainer);

    void ExchangeForces(LinkedCell *linkedCell) override;

    int GetNonBlockingStageCount() override;

    bool QueryBalanceAndExchangeNonBlocking(bool force_re_balancing, LinkedCell *linked_Cell,
                                            double etime) override;

    void BalanceAndExchange(double last_traversal_time, bool force_re_balancing,
                                    LinkedCell *linked_Cell) override;


    bool IsPositionInThisDomain(double x, double y, double z) const override;

    void GetBoundingBoxMinMax(double *min, double *max) override;

    int GetCurrentRank() const override;

    int GetTotalRank() const override;

#pragma region 通信
    void CommunicatorInit(MPI_Comm mpi_communicator, std::size_t capacity, int key = 0) override;
    void CommunicatorInitializeMemory(int numValues, int key) override;
    void CommunicatorFinalize() override;
    void CommunicatorAppendInt(int intValue) override;
    void CommunicatorAppendUnsignedLong(unsigned long unsigned_long_value) override;
    void CommunicatorAppendFloat(float floatValue) override;
    void CommunicatorAppendDouble(double doubleValue) override;
    void CommunicatorAppendLongDouble(long double longDoubleValue) override;
    int CommunicatorGetInt() override;
    unsigned long CommunicatorGetUnsignedLong() override;
    float CommunicatorGetFloat() override;
    double CommunicatorGetDouble() override;
    long double CommunicatorGetLongDouble() override;
    void CommunicatorAllreduceSum() override;
    void CommunicatorAllreduceCustom(ReduceType type) override;
    void CommunicatorScanSum() override;
    void CommunicatorBroadcast(int root) override;
#pragma endregion
    std::vector<int> GetNeighbourRanks() override;
    std::vector<std::vector<std::vector<int>>> GetAllRanks() override;


#ifdef ENABLE_MPI
    /*!
     * @brief 获取MPI的通信域 Comm
     * @return MPI通信域
     */
     MPI_Comm GetMPIComm() override;
#endif

    size_t GetTotalDataSize(unsigned localN, float* minrnd, float* maxrnd) override;

    void SetSendLeavingWithCopiesStrategy(bool send_together)  override;

    bool IsSendLeavingWithCopies() override;


};
