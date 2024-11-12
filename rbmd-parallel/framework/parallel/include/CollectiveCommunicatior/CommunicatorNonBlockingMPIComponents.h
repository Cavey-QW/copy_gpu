#pragma once
#include <mpi.h>
#include <memory>
#include "CommunicatorBase.h"
#include "CommunicatorMPIInterface.h"
class CommunicatorNonBlockingMPIComponents: public CommunicatorMPIInterface, public CommunicatorBase {
    //TODO 先声明变量方便看明白
protected:
    //! @brief 用于存储需要通信的值的向量
    std::vector<UnionData> _values;  // 目前先用Vector,后面用模板重构

    //! @brief 用于提取已通信的值的迭代器
    std::vector<UnionData>::const_iterator _getter;
    //! Vector of the corresponding MPI types for the values stored in _values
    std::vector<MPI_Datatype> _types;
    //! MPI_Datatype which will be used in the communication command and which
    //! represents all values
    MPI_Datatype _agglomerated_type;
    //! Communicator to be used by the communication commands
    MPI_Comm _mpi_communicator;

    std::unique_ptr<MPI_Request> _request{new MPI_Request()};
    MPI_Op _agglomeratedTypeAddOperator{MPI_OP_NULL};
    bool _communicationInitiated{false};
    bool _valuesValid{false};
    /// tempValues is used for overlapped communications!
    std::vector<UnionData> _tempValues;
    bool _firstComm{true};

public:
    CommunicatorNonBlockingMPIComponents();
    CommunicatorNonBlockingMPIComponents(const CommunicatorNonBlockingMPIComponents&) = delete;
    CommunicatorNonBlockingMPIComponents& operator=(const CommunicatorNonBlockingMPIComponents&) = delete;
    ~CommunicatorNonBlockingMPIComponents() override;

    void Init(MPI_Comm mpi_communicator, std::size_t capacity, int key) override;   // 如何调用到呢  类型转换？

    void Finalize() override;

    void AppendInt(int int_value) override;

    int GetInt() override;

    void AppendUnsignedLong(unsigned long unsigned_long_value) override;

    unsigned long GetUnsignedLong() override;

    void AppendFloat(float float_value) override;

    float GetFloat() override;

    void AppendDouble(double double_value) override;

    double GetDouble() override;

    void AppendLongDouble(long double long_double_value) override;

    long double GetLongDouble() override;

    MPI_Comm GetMPIComm() override;

    void Broadcast(int root) override;

    void AllReduceSumAllowPrevious() override;

    void AllReduceCustom(ReduceType type) override;

    void AllReduceSum() override;

    void ScanSum() override;

    size_t GetTotalDataSize() override;

protected:

    void InitializeMemory(std::size_t capacity) override ;

    void SetMPIType();

    static void Add(UnionData *in_vector, UnionData *in_out_vector, int,
             MPI_Datatype *data_type);

    static void Max(UnionData *in_vector, UnionData *in_out_vector, int */*len*/,
                    MPI_Datatype *data_type);

    static void Min(UnionData *in_vector, UnionData *in_out_vector, int */*len*/,
                    MPI_Datatype *data_type);


private:
    void WaitAndUpdateData();
    void InitAllreduceSum(std::vector<UnionData> values);

};

