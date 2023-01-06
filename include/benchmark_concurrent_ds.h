#pragma once
extern "C"{
    #include "curve_math.h"
	#include "types/state.h"
}
#include <tbb/concurrent_vector.h>

void benchmark_precomputation(instance_t*, shared_state_t*);
typedef tbb::concurrent_vector<CurveAndPointsSIDH> tbb_vector;
