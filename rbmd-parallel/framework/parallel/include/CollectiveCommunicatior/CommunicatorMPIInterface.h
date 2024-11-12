#pragma once

#include <mpi.h>
#include "CommunicatorBaseInterface.h"
class CommunicatorMPIInterface: public virtual CommunicatorBaseInterface{
public:
    ~CommunicatorMPIInterface() override = default;
    virtual void Init(MPI_Comm mpi_communicator, std::size_t capacity, int key = 0)  = 0;
    virtual MPI_Comm GetMPIComm() = 0;
    virtual void  AllReduceSumAllowPrevious() = 0;
};