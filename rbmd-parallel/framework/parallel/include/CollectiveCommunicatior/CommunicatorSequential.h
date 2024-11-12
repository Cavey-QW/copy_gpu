#pragma once

#include <vector>
#include "CommunicatorBase.h"

// 通信基类，有MPI NCCL RCCL HCCL ACCL 等实现
class CommunicatorSequential : CommunicatorBase {
protected:

    //! @brief 用于存储需要通信的值的向量
    std::vector<UnionData> _values;  // 目前先用Vector,后面用模板重构

    //! @brief 用于提取已通信的值的迭代器
    std::vector<UnionData>::const_iterator _getter;
public:
    void InitializeMemory(std::size_t capacity) override ;

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


    void Broadcast(int = 0) override;

    //! Performs an all-reduce (sum)
    void AllReduceSum() override;

    // doku in base class
    void AllReduceCustom(ReduceType type) override;

    //! Performs a scan (sum)
    void ScanSum() override;

    void Finalize() override;

    size_t GetTotalDataSize() override;
};

