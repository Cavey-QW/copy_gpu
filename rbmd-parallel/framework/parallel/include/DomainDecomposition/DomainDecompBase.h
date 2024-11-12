#pragma once

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <vector>
#include "Domain.h"
#include "LinkedCell.h"  //外面用的很多，这里不要作为类成员变量
#include "CommunicatorBase.h"
#include "HaloRegion.h"
#define DIM 3

class CommunicationPartner;

// Base class for domain decomposition 并行基类-------非并行的是否要在这里实现呢?   不要在这里实现，直接继承一个seq
class DomainDecompBase {
public:
    virtual ~DomainDecompBase() = default;

    virtual bool QueryBalanceAndExchangeNonBlocking(bool force_re_balancing, LinkedCell *linked_Cell,
                                                    double etime) = 0;

    virtual void BalanceAndExchange(double last_traversal_time, bool force_re_balancing,
                                    LinkedCell *linked_Cell) = 0;  //TODO 考虑Linked Cell是否要去掉 ------- 外面用的多 不要去掉

    virtual double GetBoundingBoxMin(int dimension) const = 0;

    virtual double GetBoundingBoxMax(int dimension) const = 0;

    virtual void GetBoundingBoxMinMax(double *min, double *max) = 0;

    /*!
     * @brief 获取邻居进程的编号
     * @return 26个进程，前向邻居13个后向邻居13个
     */
    virtual std::vector<int> GetNeighbourRanks() = 0;

    /*!
     * @brief 设置离开粒子与ghost粒子是否一起发送的策略
     * @param send_together 若为true则一起发送，否则单独发送
     */
    virtual void SetSendLeavingWithCopiesStrategy(bool send_together)  = 0;

    virtual void ExchangeForces(LinkedCell* linkedCell) =0;

    virtual std::vector<std::vector<std::vector<int>>> GetAllRanks() = 0;


#pragma region 通信方法
    virtual void CommunicatorInit(MPI_Comm mpi_communicator, std::size_t capacity, int key = 0) = 0;
    //! has to call init method of a CollComm class
    virtual void CommunicatorInitializeMemory(int numValues, int key ) = 0;

    //! has to call finalize method of a CollComm class
    virtual void CommunicatorFinalize() = 0;

    //! has to call appendInt method of a CollComm class
    virtual void CommunicatorAppendInt(int intValue) = 0;

    //! has to call appendUnsLong method of a CollComm class
    virtual void CommunicatorAppendUnsignedLong(unsigned long unsigned_long_value) = 0;

    //! has to call appendFloat method of a CollComm class
    virtual void CommunicatorAppendFloat(float floatValue) = 0;

    //! has to call appendDouble method of a CollComm class
    virtual void CommunicatorAppendDouble(double doubleValue) = 0;

    //! has to call appendLongDouble method of a CollComm class
    virtual void CommunicatorAppendLongDouble(long double longDoubleValue) = 0;

    //! has to call getInt method of a CollComm class
    virtual int CommunicatorGetInt() = 0;

    //! has to call getUnsLong method of a CollComm class
    virtual unsigned long CommunicatorGetUnsignedLong() = 0;

    //! has to call getFloat method of a CollComm class
    virtual float CommunicatorGetFloat() = 0;

    //! has to call getDouble method of a CollComm class
    virtual double CommunicatorGetDouble() = 0;

    //! has to call getLongDouble method of a CollComm class
    virtual long double CommunicatorGetLongDouble() = 0;

    //! has to call allreduceSum method of a CollComm class (none in sequential version)
    virtual void CommunicatorAllreduceSum() = 0;

    //! has to call allreduceSum method of a CollComm class (none in sequential version), allows for values of previous iteration.
    virtual void CommunicatorAllreduceSumAllowPrevious() = 0;

    //! has to call allreduceCustom method of a CollComm class (none in sequential version)
    virtual void CommunicatorAllreduceCustom(ReduceType type) = 0;

    //! has to call scanSum method of a CollComm class (none in sequential version)
    virtual void CommunicatorScanSum() = 0;

    //! has to call broadcast method of a CollComm class (none in sequential version)
    virtual void CommunicatorBroadcast(int root = 0) = 0;

#pragma endregion

        /*!
     * @brief 判断给定坐标是否在当前进程负责的域中
     * @param x x坐标
     * @param y y坐标
     * @param z z坐标
     * @return true表示在当前域中，否则返回false
     */
    virtual bool IsPositionInThisDomain(double x, double y, double z) const = 0;

// #ifdef ENABLE_MPI
//     /*!
//      * @brief 获取MPI的通信域 Comm
//      * @return MPI通信域
//      */
    virtual MPI_Comm& GetMPIComm() = 0;
// #endif

    /*!
     * @brief 获取当前进程编号
     * @return 当前进程编号
     */
    virtual int GetCurrentRank() const = 0;

    /*!
     * @brief 获取进程总数
     * @return 进程总数
     */
    virtual int GetTotalRank() const = 0;

    virtual MPI_Datatype getMPIParticleType() = 0;
    virtual  MPI_Datatype getMPIParticleForceType() = 0;

    /*!
     * @brief 这个函数用于DomainDecomp与力计算(may be)非阻塞进行
     * @return
     */
    virtual int GetNonBlockingStageCount() = 0;

    virtual size_t GetTotalDataSize(unsigned localN, float* minrnd, float* maxrnd) = 0;


   virtual bool IsSendLeavingWithCopies() = 0;

    virtual std::vector<CommunicationPartner> GetNeighboursFromHaloRegion(const HaloRegion& haloRegion, double cutoff) = 0;    // TODO 放在这里不合适吧


    void PopulateHaloLayerWithCopiesDirect(const HaloRegion &haloRegion,
                                         LinkedCell *moleculeContainer, bool positionCheck) const;

    void PopulateHaloLayerWithCopies(unsigned dim, LinkedCell *moleculeContainer) const;

    virtual void addLeavingMolecules(std::vector<Atom> &invalidMolecules, LinkedCell *moleculeContainer);

    void HandleDomainLeavingParticles(unsigned dim, LinkedCell *moleculeContainer) const;

    void handleForceExchange(unsigned dim, LinkedCell *moleculeContainer) const;


protected:
    //! 防止直接使用基类
    DomainDecompBase() = default;
    int _sendLeavingAndCopiesSeparately = 0;   // 不分开发送

    //!当前进程(Host or Device)的编号
    int _current_rank = 0;
    //! 总进程数
    int _total_ranks = 1;
//  //! 当前进程负责的域 是否是子域维护linkcell比较好？        --------外面不用的就直接做成员变量，外面用的就不做成员变量 ---外面用的少 直接封装成成员变量
//    Domain *domain = nullptr;     TODO 考虑是否定位器好用

    // virtual std::vector<int> getNeighbourRanksFullShell() TODO 好像弃用了



};

