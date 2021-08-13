#pragma once
#include <stdint.h>

enum class SeedEncoding { SPACED_FACTOR, HASHED, CONTIGUOUS };

struct NoFilter
{
	bool contains(uint64_t seed, uint64_t shape) const
	{
		return true;
	}
};

extern NoFilter no_filter;