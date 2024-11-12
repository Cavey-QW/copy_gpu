#pragma once
#include <cstddef>
//! Reduce Types of the allreduceCustom operation
enum ReduceType {
    SUM, MIN, MAX
};
class CommunicatorBaseInterface{
public:
    virtual ~CommunicatorBaseInterface()= default;

    virtual void InitializeMemory(std::size_t capacity) = 0;

    virtual void AppendInt(int int_value) = 0;

    virtual int GetInt() = 0;

    virtual void AppendUnsignedLong(unsigned long unsigned_long_value) = 0;

    virtual unsigned long GetUnsignedLong() = 0;

    virtual void AppendFloat(float float_value) = 0;

    virtual float GetFloat() = 0;

    virtual void AppendDouble(double double_value) = 0;

    virtual double GetDouble() = 0;

    virtual void AppendLongDouble(long double long_double_value) = 0;

    virtual long double GetLongDouble() = 0;


    virtual void Broadcast(int = 0) = 0;

//! Performs an all-reduce (sum)
    virtual void AllReduceSum() = 0;

// doku in base class
    virtual void AllReduceCustom(ReduceType type) = 0;

//! Performs a scan (sum)
    virtual void ScanSum() = 0;

    virtual void Finalize() = 0;

    virtual size_t GetTotalDataSize() = 0;
};

