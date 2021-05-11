#pragma once
#include <stdint.h>
#include <vector>

namespace Util { namespace Algo {

struct Edge {
	uint32_t v1, v2, weight;
	bool operator<(const Edge &x) const {
		return weight > x.weight;
	}
};

std::vector<int> greedy_vertex_cover(std::vector<std::vector<int>> &neighbors);

template<typename It, typename Out>
size_t merge_capped(It i0, It i1, It j0, It j1, size_t cap, Out out) {
	const ptrdiff_t m = (ptrdiff_t)cap;
	ptrdiff_t n = 0;
	size_t count = 0;
	while (n < cap) {
		if (i0 == i1) {
			const ptrdiff_t d = std::min(m - n, j1 - j0);
			std::copy(j0, j0 + d, out);
			return count + (size_t)d;
		}
		if (j0 == j1) {
			std::copy(i0, i0 + std::min(m - n, i1 - i0), out);
			return count;
		}
		if (*i0 < *j0)
			*out++ = *i0++;
		else {
			*out++ = *j0++;
			++count;
		}
		++n;
	}
	return count;
}

}}