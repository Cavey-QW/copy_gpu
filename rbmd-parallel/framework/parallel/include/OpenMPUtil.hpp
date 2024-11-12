#pragma once
#include "LogTool.h"
#if defined(_OPENMP)
#include <omp.h>

	inline int GetThreadNum()  {return omp_get_thread_num(); }
	inline int GetNumThreads() {return omp_get_num_threads();}
	inline int GetMaxThreads() {return omp_get_max_threads();}
	typedef omp_lock_t mardyn_lock_t;
	inline void SetLock(mardyn_lock_t* l) {omp_set_lock(l);}
	inline void InitLock(mardyn_lock_t* l) {omp_init_lock(l);}
	inline void UnsetLock(mardyn_lock_t* l) {omp_unset_lock(l);}
	inline void DestoryLock(mardyn_lock_t* l) {omp_destroy_lock(l);}

#else

inline int GetThreadNum()  {return 0;}
inline int GetNumThreads() {return 1;}
inline int GetMaxThreads() {return 1;}
typedef int mardyn_lock_t;
inline void SetLock(mardyn_lock_t* l) {
    //assert(*l == 0);
    if (*l != 0){
        GlobalLogger::error("*l != 0 in SetLock");
    }
    *l = 1;
}
inline void InitLock(mardyn_lock_t* l) {
//    assert(l != nullptr);
    if (l == nullptr){
        GlobalLogger::error("l == nullptr in InitLock");
    }
}
inline void UnsetLock(mardyn_lock_t* l) {
    //assert(*l == 1);
    if (*l != 1){
        GlobalLogger::error("*l != 1 in UnsetLock");
    }
    *l = 0;
}
inline void DestoryLock(mardyn_lock_t* l) {
    //assert(l != nullptr);
    if (l == nullptr){
        GlobalLogger::error("l == nullptr in DestoryLock");

    }
}
#endif