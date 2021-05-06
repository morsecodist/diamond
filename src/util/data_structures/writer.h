template<typename T>
struct Writer {
	Writer& operator*() {
		return *this;
	}
	virtual Writer& operator=(const T& v) = 0;
};