/** @file inv_list_xbv_csv.h
    @brief Cross-bvector (xbv) research support for inv_list.
    @note Development-only; included by inv_list.cpp after shared options and I/O.
*/
/*
Copyright(c) 2002-2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For more information please visit:  http://bitmagic.io
*/

// Minimal CSV parsing for inv_list analysis output. The analysis CSV format
// contains no quoted fields; semicolons are used for lists inside a field.
#ifndef INV_LIST_XBV_CSV_H
#define INV_LIST_XBV_CSV_H

namespace inv_list_csv
{

inline std::vector<std::string> split(const std::string& line)
{
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;)
    {
        const std::size_t comma = line.find(',', begin);
        fields.push_back(line.substr(begin, comma - begin));
        if (comma == std::string::npos)
            break;
        begin = comma + 1;
    }
    return fields;
}

inline std::size_t column(const std::vector<std::string>& header,
                          const char* name)
{
    const auto it = std::find(header.begin(), header.end(), name);
    if (it == header.end())
        throw std::runtime_error(std::string("CSV is missing column: ") + name);
    return std::size_t(it - header.begin());
}

inline bm::id64_t parse_u64(const std::string& text, std::size_t line_no,
                            const char* field)
{
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
    {
        std::ostringstream msg;
        msg << "Invalid " << field << " at CSV line " << line_no;
        throw std::runtime_error(msg.str());
    }
    std::size_t pos = 0;
    const auto value = std::stoull(text, &pos);
    if (pos != text.size())
        throw std::runtime_error("Invalid integer in CSV");
    return bm::id64_t(value);
}

} // namespace inv_list_csv

#endif
