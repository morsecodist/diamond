#pragma once
#include <algorithm>

template<typename _t = size_t>
struct partition
{
	_t items, parts, size, remainder;
	partition() : items(0), parts(0), size(0), remainder(0)
	{ }
	partition(_t items, _t parts) : items(items), parts(std::min(parts, items))
	{
		if (this->parts > 0) {
			size = items / this->parts;
			remainder = items % this->parts;
		}
		else {
			size = 0;
			remainder = 0;
		}
	}
	_t getMin(_t i) const
	{
		_t b = std::min(i, remainder); return b * (size + 1) + (i - b) * size;
	}
	_t getMax(_t i) const
	{
		return getMin(i) + getCount(i);
	}
	_t getCount(_t i) const
	{
		return i < remainder ? (size + 1) : size;
	}
};
