#pragma once
#define MPI_CHECK(x) do {                   \
    int __ret;                              \
    if (MPI_SUCCESS != (__ret = (x)))       \
    std::cout << "MPI returned with error code " << __ret << std::endl;  \
} while (0)
