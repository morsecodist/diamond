/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#pragma once

#include <list>
#include <vector>
#include "../../basic/config.h"

template<typename _t>
struct Deque {

	typedef std::vector<_t> Bucket;

	Deque():
		max_size_(config.deque_bucket_size / sizeof(_t))
	{
		buckets.emplace_back();
		buckets.back().reserve(max_size_);
	}

	void push_back(const _t* ptr, size_t n) {
		if (buckets.back().size() + n > max_size_) {
			buckets.emplace_back();
			buckets.back().reserve(max_size_);
		}
		buckets.back().insert(buckets.back().end(), ptr, ptr + n);
	}

	size_t size() const {
		size_t n = 0;
		for (const Bucket& b : buckets)
			n += b.size();
		return n;
	}

	void move(std::vector<_t>& dst) {
		if (buckets.size() == 1 && dst.empty())
			dst = std::move(buckets.front());
		else {
			for (const Bucket& b : buckets)
				dst.insert(dst.end(), b.begin(), b.end());
		}
		buckets.clear();
	}

private:

	std::list<Bucket> buckets;
	const size_t max_size_;

};