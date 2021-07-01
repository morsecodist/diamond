#pragma once
#include "../util/enum.h"

namespace DP {

	enum class Flags { NONE = 0, TRACEBACK = 1, PARALLEL = 2, FULL_MATRIX = 4, WITH_COORDINATES = 8 };

	DEF_ENUM_FLAG_OPERATORS(Flags)

}

enum class HspValues : unsigned {
	NONE = 0,
	TRANSCRIPT = 1,
	QUERY_START = 1 << 1,
	QUERY_END = 1 << 2,
	TARGET_START = 1 << 3,
	TARGET_END = 1 << 4,
	STATS = 1 << 5,
	QUERY_COORDS = QUERY_START | QUERY_END,
	TARGET_COORDS = TARGET_START | TARGET_END,
	COORDS = QUERY_COORDS | TARGET_COORDS,
	STATS_OR_COORDS = STATS | QUERY_COORDS | TARGET_COORDS,
	STATS_OR_TRANSCRIPT = STATS | TRANSCRIPT
};

DEF_ENUM_FLAG_OPERATORS(HspValues)

static inline bool have_coords(const HspValues v) {
	return flag_any(v, HspValues::TRANSCRIPT) || flag_all(v, HspValues::QUERY_COORDS | HspValues::TARGET_COORDS);
}