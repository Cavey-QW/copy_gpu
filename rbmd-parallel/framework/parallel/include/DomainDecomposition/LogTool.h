#pragma once
//#include "mpi/mpi.h"
#include <mpi.h>
#include <iostream>
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE //必须定义这个宏,才能输出文件名和行号
inline void __rbmd_assert__(const char * expr, const char* file, int line) {
    std::cout << "Assertion \"" << expr << "\" failed at " << file << ":" << line << std::endl;
    exit(1);
}


#define rbmd_assert(EXPRESSION) ((EXPRESSION) ? (void)0 : __rbmd_assert__(#EXPRESSION, __FILE__, __LINE__))


#include <string>
//set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [thread %t] %v (%s:%#)");
namespace GlobalLogger {

    static int GetRank(){
        int current_rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&current_rank);
        return current_rank;
    }

    template<typename T>
    static void info(T& info) {
        if (GetRank() == 0) {
            //spdlog::info(info);
            std::cout << info << std::endl;
        }
    }


    template<typename T>
    static void warn(T& warn) {
        if (GetRank() == 0) {
            std::cout << warn << std::endl;
        }
    }

    template<typename T>
    static void error(T& error) {
        if (GetRank() == 0) {
            std::cout << error << std::endl;
            throw std::runtime_error("log error");
        }
    }
}


namespace Logger {

    template<typename T>
    static void info(T info) {
        std::cout << info << std::endl;

    }

    template<typename T>
    static void warn(T warn) {

        std::cout << warn << std::endl;

    }

    template<typename T>
    static void error(T error) {
        std::cout << error << std::endl;
        throw std::runtime_error("log error");
    }
}





