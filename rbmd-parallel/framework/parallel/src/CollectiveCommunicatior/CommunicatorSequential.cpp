#include "CommunicatorSequential.h"

void CommunicatorSequential::InitializeMemory(std::size_t capacity) {
    _values.reserve(capacity);
    _getter = _values.begin();
}

void CommunicatorSequential::AppendInt(int int_value) {       //TODO 真的需要用吗？  为了继承故意弄的 不需要
    UnionData toPush{};
    toPush._int_vaule = int_value;
    _values.push_back(toPush);
}

int CommunicatorSequential::GetInt() {
    return (_getter++)->_int_vaule;
}

void CommunicatorSequential::AppendUnsignedLong(unsigned long unsigned_long_value) {
    UnionData toPush;
    toPush._unsigned_long_value = unsigned_long_value;
    _values.push_back(toPush);
}

unsigned long CommunicatorSequential::GetUnsignedLong() {
    return (_getter++)->_unsigned_long_value;
}

void CommunicatorSequential::AppendFloat(float float_value) {
    UnionData toPush;
    toPush._float_value = float_value;
    _values.push_back(toPush);
}

float CommunicatorSequential::GetFloat() {
    return (_getter++)->_float_value;
}

void CommunicatorSequential::AppendDouble(double double_value) {
    UnionData toPush;
    toPush._double_value = double_value;
    _values.push_back(toPush);
}

double CommunicatorSequential::GetDouble() {
    return (_getter++)->_double_value;
}

void CommunicatorSequential::AppendLongDouble(long double long_double_value) {
    UnionData toPush;
    toPush._long_double_value = long_double_value;
    _values.push_back(toPush);
}

long double CommunicatorSequential::GetLongDouble() {
    return (_getter++)->_long_double_value;
}

void CommunicatorSequential::Broadcast(int) {
}

void CommunicatorSequential::AllReduceSum() {
}

void CommunicatorSequential::AllReduceCustom(ReduceType type) {
}

void CommunicatorSequential::ScanSum() {
}

void CommunicatorSequential::Finalize() {
    _values.clear();
}

size_t CommunicatorSequential::GetTotalDataSize() {
    return _values.capacity() * sizeof(UnionData);
}










