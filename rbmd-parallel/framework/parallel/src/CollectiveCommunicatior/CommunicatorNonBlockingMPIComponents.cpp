#include "LogTool.h"
#include "CommunicatorNonBlockingMPIComponents.h"
#include "CommunicatorBaseInterface.h"
#include "MPIChecker.hpp"

CommunicatorNonBlockingMPIComponents::CommunicatorNonBlockingMPIComponents() {
    _agglomerated_type = MPI_DATATYPE_NULL;
    _mpi_communicator = MPI_COMM_NULL;
}

CommunicatorNonBlockingMPIComponents::~CommunicatorNonBlockingMPIComponents() {
    if (_request) {
        if (_communicationInitiated) {
            MPI_Wait(_request.get(), MPI_STATUS_IGNORE);
        }
    }
    if (_agglomeratedTypeAddOperator != MPI_OP_NULL) {
        MPI_CHECK(MPI_Op_free(&_agglomeratedTypeAddOperator));
    }
    if (_agglomerated_type != MPI_DATATYPE_NULL) {
        MPI_CHECK(MPI_Type_free(&_agglomerated_type));
    }
}


void CommunicatorNonBlockingMPIComponents::InitializeMemory(std::size_t capacity) {
    _values.reserve(capacity);
    _getter = _values.begin();
}

void CommunicatorNonBlockingMPIComponents::Init(MPI_Comm mpi_communicator, std::size_t capacity, int key) {
    InitializeMemory(capacity);
    _mpi_communicator = mpi_communicator;
    _types.reserve(capacity);
    if (!_firstComm) {
        if(static_cast<int>(_values.size()) != capacity){
            GlobalLogger::error("_values.size()) != capacity in CommunicatorNonBlockingMPIComponents::Init");
        }
    }
}

void CommunicatorNonBlockingMPIComponents::Finalize() {
    _values.clear();
    _types.clear();
}

void CommunicatorNonBlockingMPIComponents::AppendInt(int int_value) {
    UnionData toPush{};
    toPush._int_vaule = int_value;
    _values.push_back(toPush);
    _types.push_back(MPI_INT);
}

int CommunicatorNonBlockingMPIComponents::GetInt() {
    return (_getter++)->_int_vaule;
}

void CommunicatorNonBlockingMPIComponents::AppendUnsignedLong(unsigned long unsigned_long_value) {
    UnionData toPush{};
    toPush._unsigned_long_value = unsigned_long_value;
    _values.push_back(toPush);
    _types.push_back(MPI_UNSIGNED_LONG);
}

unsigned long CommunicatorNonBlockingMPIComponents::GetUnsignedLong() {
    return (_getter++)->_unsigned_long_value;
}

void CommunicatorNonBlockingMPIComponents::AppendFloat(float float_value) {
    UnionData toPush{};
    toPush._float_value = float_value;
    _values.push_back(toPush);
    _types.push_back(MPI_FLOAT);
}

float CommunicatorNonBlockingMPIComponents::GetFloat() {
    return (_getter++)->_float_value;
}

void CommunicatorNonBlockingMPIComponents::AppendDouble(double double_value) {
    UnionData toPush{};
    toPush._double_value = double_value;
    _values.push_back(toPush);
    _types.push_back(MPI_DOUBLE);
}

double CommunicatorNonBlockingMPIComponents::GetDouble() {
    return (_getter++)->_double_value;
}

void CommunicatorNonBlockingMPIComponents::AppendLongDouble(long double long_double_value) {
    UnionData toPush{};
    toPush._long_double_value = long_double_value;
    _values.push_back(toPush);
    _types.push_back(MPI_LONG_DOUBLE);
}

long double CommunicatorNonBlockingMPIComponents::GetLongDouble() {
    return (_getter++)->_long_double_value;
}

MPI_Comm CommunicatorNonBlockingMPIComponents::GetMPIComm() {
    return _mpi_communicator;
}

void CommunicatorNonBlockingMPIComponents::SetMPIType() {
    int numblocks = static_cast<int>(_values.size());   //TODO Narrow convert
    std::vector<int> blocklengths(numblocks);
    std::vector<MPI_Aint> disps(numblocks);
    int disp = 0;
    for (int i = 0; i < numblocks; i++) {
        blocklengths[i] = 1;
        disps[i] = disp;
        disp += sizeof(UnionData);
    }
    MPI_Datatype *startOfTypes = &(_types[0]);
    MPI_CHECK( //TODO 只考虑了mpi 3 以上的版本
            MPI_Type_create_struct(numblocks, blocklengths.data(), disps.data(),
                                   startOfTypes, &_agglomerated_type));
    MPI_CHECK(MPI_Type_commit(&_agglomerated_type));

}

void CommunicatorNonBlockingMPIComponents::Broadcast(int root) {
    SetMPIType();
    UnionData *start_of_values = _values.data();
    MPI_CHECK(
            MPI_Bcast(start_of_values, 1, _agglomerated_type, root,
                      _mpi_communicator));
    MPI_CHECK(MPI_Type_free(&_agglomerated_type));

}

void CommunicatorNonBlockingMPIComponents::Add(UnionData *in_vector, UnionData *in_out_vector, int,
                                               MPI_Datatype *data_type) {
    int num_ints;
    int num_address;
    int num_types;
    int combiner;

    MPI_CHECK(
            MPI_Type_get_envelope(*data_type, &num_ints, &num_address, &num_types,
                                  &combiner));

    std::vector<int> vector_ints(num_ints);
    std::vector<MPI_Aint> vector_address(num_address);
    std::vector<MPI_Datatype> vector_types(num_types);

    MPI_CHECK(
            MPI_Type_get_contents(*data_type, num_ints, num_address, num_types,
                                  vector_ints.data(), vector_address.data(), vector_types.data()));

    for (int i = 0; i < num_types; i++) {
        if (vector_types[i] == MPI_INT) {
            in_out_vector[i]._int_vaule += in_vector[i]._int_vaule;
        } else if (vector_types[i] == MPI_UNSIGNED_LONG) {
            in_out_vector[i]._unsigned_long_value += in_vector[i]._unsigned_long_value;
        } else if (vector_types[i] == MPI_FLOAT) {
            in_out_vector[i]._float_value += in_vector[i]._float_value;
        } else if (vector_types[i] == MPI_DOUBLE) {
            in_out_vector[i]._double_value += in_vector[i]._double_value;
        } else if (vector_types[i] == MPI_LONG_DOUBLE) {
            in_out_vector[i]._long_double_value += in_vector[i]._long_double_value;
        }
    }
}

void CommunicatorNonBlockingMPIComponents::Max(UnionData *in_vector, UnionData *in_out_vector, int *, MPI_Datatype *data_type) {
    int num_ints;
    int num_address;
    int num_types;
    int combiner;

    MPI_CHECK(
            MPI_Type_get_envelope(*data_type, &num_ints, &num_address, &num_types,
                                  &combiner));

    std::vector<int> vector_ints(num_ints);
    std::vector<MPI_Aint> vector_address(num_address);
    std::vector<MPI_Datatype> vector_types(num_types);

    MPI_CHECK(
            MPI_Type_get_contents(*data_type, num_ints, num_address, num_types,
                                  vector_ints.data(), vector_address.data(), vector_types.data()));

    for (int i = 0; i < num_types; i++) {
        if (vector_types[i] == MPI_INT) {
            in_out_vector[i]._int_vaule = std::max(in_out_vector[i]._int_vaule, in_vector[i]._int_vaule);
        } else if (vector_types[i] == MPI_UNSIGNED_LONG) {
            in_out_vector[i]._unsigned_long_value = std::max(in_out_vector[i]._unsigned_long_value,
                                                             in_vector[i]._unsigned_long_value);
        } else if (vector_types[i] == MPI_FLOAT) {
            in_out_vector[i]._float_value = std::max(in_out_vector[i]._float_value, in_vector[i]._float_value);
        } else if (vector_types[i] == MPI_DOUBLE) {
            in_out_vector[i]._double_value = std::max(in_out_vector[i]._double_value, in_vector[i]._double_value);
        } else if (vector_types[i] == MPI_LONG_DOUBLE) {
            in_out_vector[i]._long_double_value = std::max(in_out_vector[i]._long_double_value,
                                                           in_vector[i]._long_double_value);
        }
    }
}

void CommunicatorNonBlockingMPIComponents::Min(UnionData *in_vector, UnionData *in_out_vector, int *, MPI_Datatype *data_type) {
    int num_ints;
    int num_address;
    int num_types;
    int combiner;

    MPI_CHECK(
            MPI_Type_get_envelope(*data_type, &num_ints, &num_address, &num_types,
                                  &combiner));

    std::vector<int> vector_ints(num_ints);
    std::vector<MPI_Aint> vector_address(num_address);
    std::vector<MPI_Datatype> vector_types(num_types);

    MPI_CHECK(
            MPI_Type_get_contents(*data_type, num_ints, num_address, num_types,
                                  vector_ints.data(), vector_address.data(), vector_types.data()));

    for (int i = 0; i < num_types; i++) {
        if (vector_types[i] == MPI_INT) {
            in_out_vector[i]._int_vaule = std::min(in_out_vector[i]._int_vaule, in_vector[i]._int_vaule);
        } else if (vector_types[i] == MPI_UNSIGNED_LONG) {
            in_out_vector[i]._unsigned_long_value = std::min(in_out_vector[i]._unsigned_long_value,
                                                             in_vector[i]._unsigned_long_value);
        } else if (vector_types[i] == MPI_FLOAT) {
            in_out_vector[i]._float_value = std::min(in_out_vector[i]._float_value, in_vector[i]._float_value);
        } else if (vector_types[i] == MPI_DOUBLE) {
            in_out_vector[i]._double_value = std::min(in_out_vector[i]._double_value, in_vector[i]._double_value);
        } else if (vector_types[i] == MPI_LONG_DOUBLE) {
            in_out_vector[i]._long_double_value = std::min(in_out_vector[i]._long_double_value,
                                                           in_vector[i]._long_double_value);
        }
    }
}


void CommunicatorNonBlockingMPIComponents::AllReduceCustom(ReduceType type) {
    // 启用聚集归约操作。这将把所有值存储在一个数组中，并应用用户定义的归约操作，以便 MPI 归约操作只被调用一次。
    SetMPIType();
    MPI_Op agglomeratedTypeAddOperator;
    const int commutative = 1;
    UnionData *startOfValues = _values.data();

    switch (type) {
        case ReduceType::SUM:
            MPI_CHECK(
                    MPI_Op_create(
                            (MPI_User_function *) CommunicatorNonBlockingMPIComponents::Add,
                            commutative, &agglomeratedTypeAddOperator));
            break;
        case ReduceType::MAX:
            MPI_CHECK(
                    MPI_Op_create(
                            (MPI_User_function *) CommunicatorNonBlockingMPIComponents::Max,
                            commutative, &agglomeratedTypeAddOperator));
            break;
        case ReduceType::MIN:
            MPI_CHECK(
                    MPI_Op_create(
                            (MPI_User_function *) CommunicatorNonBlockingMPIComponents::Min,
                            commutative, &agglomeratedTypeAddOperator));
            break;
        default:
            exit(1);
    }

    MPI_CHECK(
            MPI_Allreduce(MPI_IN_PLACE, startOfValues, 1, _agglomerated_type, agglomeratedTypeAddOperator,
                          _mpi_communicator));
    MPI_CHECK(MPI_Op_free(&agglomeratedTypeAddOperator));
    MPI_CHECK(MPI_Type_free(&_agglomerated_type));
}

void CommunicatorNonBlockingMPIComponents::AllReduceSum() {
    AllReduceCustom(ReduceType::SUM);
}

void CommunicatorNonBlockingMPIComponents::WaitAndUpdateData() {
    MPI_CHECK(MPI_Wait(_request.get(), MPI_STATUS_IGNORE));
    // copy the temporary values to the real values!
    _values = _tempValues;

    _communicationInitiated = false;
    _valuesValid = true;
}

void CommunicatorNonBlockingMPIComponents::InitAllreduceSum(std::vector<UnionData> values) {
    // copy the values to a temporary vector:
    // all nonblocking operations have to be performed on this vector!
    // this is necessary to maintain the validity of the data from previous steps
    _tempValues = values;
    if (_agglomerated_type == MPI_DATATYPE_NULL) {
        SetMPIType();
    }
    const int commutative = 1;
    UnionData *startOfValues = &(_tempValues[0]);
    if (_agglomeratedTypeAddOperator == MPI_OP_NULL) {
        MPI_CHECK(
                MPI_Op_create((MPI_User_function *) CommunicatorNonBlockingMPIComponents::Add, commutative,
                              &_agglomeratedTypeAddOperator));
    }
    MPI_CHECK(
            MPI_Iallreduce(MPI_IN_PLACE, startOfValues, 1, _agglomerated_type, _agglomeratedTypeAddOperator,
                           _mpi_communicator, _request.get()));
    _communicationInitiated = true;
}


void CommunicatorNonBlockingMPIComponents::AllReduceSumAllowPrevious() {
    if (!_valuesValid) {
        if (_communicationInitiated) {
            // if the values are not valid and communication is already initiated, first wait for the communication to finish.
            WaitAndUpdateData();
        }
        // if the values were invalid, we have to do a proper allreduce!
        AllReduceSum();
        _tempValues = _values;  // save it somewhere safe for next iteration!
        _valuesValid = true;
    } else {

        if (_communicationInitiated) {
            // if the values are not valid and communication is already initiated, first wait for the communication to finish.
            // safe the data, that needs to be sent
            std::vector<UnionData> toSendValues = _values;
            // get the new data, this updates _values
            WaitAndUpdateData();
            // initiate the next allreduce
            InitAllreduceSum(toSendValues);
        } else {
            std::vector<UnionData> previous = _tempValues;
            // initiate the next allreduce
            InitAllreduceSum(_values);  // this updates _tempValues
            // restore saved data from previous operations.
            _values = previous;
        }

    }
}

void CommunicatorNonBlockingMPIComponents::ScanSum() {
    SetMPIType();
    MPI_Op agglomeratedTypeAddOperator;
    const int commutative = 1;
    UnionData *startOfValues = _values.data();
    MPI_CHECK(
            MPI_Op_create(
                    (MPI_User_function *) CommunicatorNonBlockingMPIComponents::Add,
                    commutative, &agglomeratedTypeAddOperator));
    MPI_CHECK(
            MPI_Scan(MPI_IN_PLACE, startOfValues, 1, _agglomerated_type, agglomeratedTypeAddOperator,
                     _mpi_communicator));
    MPI_CHECK(MPI_Op_free(&agglomeratedTypeAddOperator));
    MPI_CHECK(MPI_Type_free(&_agglomerated_type));
}

size_t CommunicatorNonBlockingMPIComponents::GetTotalDataSize() {
    return _values.capacity() * sizeof(UnionData) + _types.capacity() * sizeof(MPI_Datatype);
}


















