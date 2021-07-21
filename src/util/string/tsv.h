#pragma once

#include <string>
#include "../io/text_input_file.h"

namespace Util { namespace Tsv {

std::string fetch_block(TextInputFile& f, std::string& buf);
std::vector<std::string> extract_column(const std::string& buf, const size_t i);

}}