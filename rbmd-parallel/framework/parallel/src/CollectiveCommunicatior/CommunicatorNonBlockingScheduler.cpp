#include "LogTool.h"
#include "CommunicatorNonBlockingScheduler.h"
#include <ciso646>
CommunicatorNonBlockingScheduler::CommunicatorNonBlockingScheduler() {
    _currentKey = -1;
}

void CommunicatorNonBlockingScheduler::Init(MPI_Comm mpi_communicator, std::size_t capacity, int key) {
    if (_currentKey != -1) {
        GlobalLogger::error("current key is not -1 内部错误");
    }

    _currentKey = key;

    // add the key, if it is not yet existent:
    if (_comms.count(_currentKey) == 1) {
        // this happens, if the key is already existent.
//        Log::global_log->debug() << "CollectiveCommunicationNonBlocking: key " << _currentKey
//                                 << " already existent. Reusing information." << std::endl;
    } else {
//        Log::global_log->debug() << "CollectiveCommunicationNonBlocking: key " << _currentKey
//                                 << " not existent. Cannot reuse information." << std::endl;
        // Creates the CollectiveCommunicationSingleNonBlocking object
        auto [_, inserted] = _comms.try_emplace(_currentKey);
        if (not inserted) {
//            string ster << "CollectiveCommunicationNonBlocking: key " << _currentKey
//                                     << " could not be inserted. Aborting!" << std::endl;
            GlobalLogger::error("CollectiveCommunicationNonBlocking: key could not be inserted. In CommunicatorNonBlockingScheduler::Init");
        }
    }
    _comms.at(_currentKey).Init(mpi_communicator, capacity, _currentKey);
}

void CommunicatorNonBlockingScheduler::Finalize() {
    _comms.at(_currentKey).Finalize();
    if (_currentKey == 0) {
//        Log::global_log->debug() << "CollectiveCommunicationNonBlocking: finalizing with key " << _currentKey
//                                 << ", thus the entry is removed." << std::endl;
        _comms.erase(_currentKey);
    }
    _currentKey = -1;
}

void CommunicatorNonBlockingScheduler::AppendInt(int int_value) {
    _comms.at(_currentKey).AppendInt(int_value);
}

int CommunicatorNonBlockingScheduler::GetInt() {
    return _comms.at(_currentKey).GetInt();
}

void CommunicatorNonBlockingScheduler::AppendUnsignedLong(unsigned long unsigned_long_value) {
    _comms.at(_currentKey).AppendUnsignedLong(unsigned_long_value);
}

unsigned long CommunicatorNonBlockingScheduler::GetUnsignedLong() {
    return _comms.at(_currentKey).GetUnsignedLong();
}

void CommunicatorNonBlockingScheduler::AppendFloat(float float_value) {
    _comms.at(_currentKey).AppendFloat(float_value);

}

float CommunicatorNonBlockingScheduler::GetFloat() {
    return _comms.at(_currentKey).GetFloat();
}

void CommunicatorNonBlockingScheduler::AppendDouble(double double_value) {
    _comms.at(_currentKey).AppendDouble(double_value);

}

double CommunicatorNonBlockingScheduler::GetDouble() {
    return _comms.at(_currentKey).GetDouble();
}

void CommunicatorNonBlockingScheduler::AppendLongDouble(long double long_double_value) {
    _comms.at(_currentKey).AppendLongDouble(long_double_value);

}

long double CommunicatorNonBlockingScheduler::GetLongDouble() {
    return _comms.at(_currentKey).GetLongDouble();
}

MPI_Comm CommunicatorNonBlockingScheduler::GetMPIComm() {
    return _comms.at(_currentKey).GetMPIComm();
}

void CommunicatorNonBlockingScheduler::Broadcast(int root) {
    _comms.at(_currentKey).Broadcast(root);
}

void CommunicatorNonBlockingScheduler::AllReduceSum() {
//Log::global_log->debug() << "CollectiveCommunicationNonBlocking: normal Allreduce" << std::endl;
    _comms.at(_currentKey).AllReduceSum();
}

void CommunicatorNonBlockingScheduler::AllReduceSumAllowPrevious() {
    //Log::global_log->debug() << "CollectiveCommunicationNonBlocking: nonblocking Allreduce with id "<< _currentKey << std::endl;
    if (_currentKey <= 0){
        GlobalLogger::error("currentkey <=0 in CommunicatorNonBlockingScheduler::AllReduceSumAllowPrevious");
    }
    _comms.at(_currentKey).AllReduceSumAllowPrevious();
}

void CommunicatorNonBlockingScheduler::AllReduceCustom(ReduceType type) {
    _comms.at(_currentKey).AllReduceCustom(type);
}

void CommunicatorNonBlockingScheduler::ScanSum() {
    _comms.at(_currentKey).ScanSum();

}

size_t CommunicatorNonBlockingScheduler::GetTotalDataSize() {
    size_t tmp = 0;
    for (auto& comm : _comms){
        tmp += comm.second.GetTotalDataSize();
    }
    return tmp;
}

void CommunicatorNonBlockingScheduler::InitializeMemory(std::size_t capacity) {

}






