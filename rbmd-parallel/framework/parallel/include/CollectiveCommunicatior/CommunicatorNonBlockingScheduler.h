#pragma once
#include "CommunicatorMPIInterface.h"
#include "CommunicatorNonBlockingMPIComponents.h"
#include <map>
class CommunicatorNonBlockingScheduler: public CommunicatorMPIInterface {
private:
    int _currentKey;
    std::map<int, CommunicatorNonBlockingMPIComponents> _comms;
public:
    CommunicatorNonBlockingScheduler();
    ~CommunicatorNonBlockingScheduler() override = default;

    void Init(MPI_Comm mpi_communicator, std::size_t capacity, int key) override;

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

    void InitializeMemory(std::size_t capacity) override;

};

