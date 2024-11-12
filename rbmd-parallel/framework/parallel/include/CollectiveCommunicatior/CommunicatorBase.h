#pragma once
#include <mpi.h>
#include <vector>
#include "CommunicatorBaseInterface.h"

// 通信基类，有MPI NCCL RCCL HCCL ACCL 等实现   TODO 目录结构要重新调整
class CommunicatorBase : CommunicatorBaseInterface {
public:
    ~CommunicatorBase() override = default;
protected:
    CommunicatorBase() = default;
    /*!
    * @brief 用于通信的数据包，包含不同的类型
    */
    union UnionData {
        int _int_vaule;
        unsigned long _unsigned_long_value;
        float _float_value;
        double _double_value;
        long double _long_double_value;
    };
};

