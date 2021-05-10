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
#include <iterator>
#include "../parallel/mutex.h"
#include "writer.h"

template<typename T, typename Sync = ::Sync>
struct Deque {

	typedef std::vector<T> Bucket;

	Deque(size_t bucket_size = 4 * (1llu << 30)) :
		max_size_(bucket_size / sizeof(T))
	{
		new_bucket();
	}

	void reserve(size_t n) {}

	void push_back(const T& x) {
		if (buckets.back().size() >= max_size_) {
			new_bucket();
		}
		buckets.back().push_back(x);
	}

	void push_back(const T* ptr, size_t n) {
		if (buckets.back().size() + n > max_size_) {
			new_bucket();
		}
		buckets.back().insert(buckets.back().end(), ptr, ptr + n);
	}

	template<typename It>
	void push_back(It begin, It end) {
		mtx_.lock();
		if (buckets.back().size() + (end - begin) > max_size_)
			new_bucket();
		buckets.back().insert(buckets.back().end(), begin, end);
		mtx_.unlock();
	}

	size_t size() const {
		size_t n = 0;
		for (const Bucket& b : buckets)
			n += b.size();
		return n;
	}

	void move(std::vector<T>& dst) {
		if (buckets.size() == 1 && dst.empty())
			dst = std::move(buckets.front());
		else {
			for (const Bucket& b : buckets)
				dst.insert(dst.end(), b.begin(), b.end());
		}
		buckets.clear();
	}

	struct Iterator {

		using iterator_category = std::random_access_iterator_tag;
		using difference_type = ptrdiff_t;
		using value_type = T;
		using pointer = T*;
		using reference = T&;

		Iterator() {}

		Iterator(size_t i, T** data) :
			i_(i),
			data_(data)
		{
		}

		T& operator*() {
			return data_[i_ >> 29][i_ & (0x20000000 - 1)];
		}

		T* operator->() {
			return &data_[i_ >> 29][i_ & (0x20000000 - 1)];
		}

		T& operator*() const {
			return data_[i_ >> 29][i_ & (0x20000000 - 1)];
		}

		T& operator[](ptrdiff_t i) {
			ptrdiff_t j = i_ + i;
			return data_[j >> 29][j & (0x20000000 - 1)];
		}

		ptrdiff_t operator-(Iterator& it) {
			return i_ - it.i_;
		}

		ptrdiff_t operator-(const Iterator& it) const {
			return i_ - it.i_;
		}

		Iterator operator+(ptrdiff_t i) const {
			return Iterator(i_ + i, data_);
		}

		Iterator operator-(ptrdiff_t i) const {
			return Iterator(i_ - i, data_);
		}

		Iterator& operator++() {
			++i_;
			return *this;
		}

		Iterator operator++(int) {
			Iterator r(i_, data_);
			++i_;
			return r;
		}

		Iterator& operator--() {
			--i_;
			return *this;
		}

		bool operator==(const Iterator& it) const {
			return i_ == it.i_;
		}

		bool operator!=(const Iterator& it) const {
			return i_ != it.i_;
		}

		bool operator>=(const Iterator& it) const {
			return i_ >= it.i_;
		}

		bool operator<=(const Iterator& it) const {
			return i_ <= it.i_;
		}

		bool operator<(const Iterator& it) const {
			return i_ < it.i_;
		}

		bool operator>(const Iterator& it) const {
			return i_ > it.i_;
		}
		
		Iterator operator-(int i) const {
			return Iterator(i_ - i, data_);
		}

		Iterator& operator+=(ptrdiff_t i) {
			i_ += i;
			return *this;
		}

	private:

		ptrdiff_t i_;
		T** data_;

	};

	Iterator begin() {
		init();
		return Iterator(0, data_.data());
	}

	Iterator end() {
		init();
		return Iterator(total_, data_.data());
	}

private:

	void init() {
		data_.clear();
		total_ = 0;
		for (Bucket& b : buckets) {
			data_.push_back(b.data());
			total_ += b.size();
		}
	}

	void new_bucket() {
		buckets.emplace_back();
		buckets.back().reserve(max_size_);
	}

	std::list<Bucket> buckets;
	std::vector<T*> data_;
	const size_t max_size_;
	size_t total_;
	Mutex<Sync> mtx_;

};

template<typename T>
struct AsyncWriter : public Writer<T> {

	AsyncWriter(Deque<T, Async>& dst):
		dst_(&dst)
	{}

	virtual AsyncWriter& operator=(const T& v) override {
		buf_.push_back(v);
		if (buf_.size() >= BUF_SIZE) {
			dst_->push_back(buf_.begin(), buf_.end());
			buf_.clear();
		}
		return *this;
	}

	virtual ~AsyncWriter() {
		dst_->push_back(buf_.begin(), buf_.end());
	}

private:

	static const size_t BUF_SIZE = 4096;

	Deque<T, Async>* dst_;
	std::vector<T> buf_;

};