#include <iostream>
#include <unordered_map>
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/string/tsv.h"
#include "../util/string/tokenizer.h"

using std::cout;
using std::endl;
using namespace Util::Tsv;

void join() {
	TextInputFile file1(config.file1);
	TextInputFile file2(config.file2);
	std::unordered_map<string, string> values;

	while (file1.getline(), !file1.line.empty() || !file1.eof()) {
		values[Util::Tsv::column(file1.line, 0)] = Util::Tsv::column(file1.line, 1);
	}

	string key;
	while (file2.getline(), !file2.line.empty() || !file2.eof()) {		
		Util::String::Tokenizer tok(file2.line, "\t");		
		tok >> key;
		cout << key << '\t' << values.at(key) << '\t' << tok.ptr() << std::endl;
	}

	file1.close();
	file2.close();
}