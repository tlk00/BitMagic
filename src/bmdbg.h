#ifndef BMDBG__H__INCLUDED__
#define BMDBG__H__INCLUDED__
/*
Copyright(c) 2002-2009 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

For more information please visit:  http://bmagic.sourceforge.net

*/

// BitMagic debugging functions (internal header)

#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <memory.h>
#include <time.h>

#include <iostream>
#include <strstream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "bmdef.h"



#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4311 4312 4127)
#endif

namespace bm
{

inline
void PrintGap(const bm::gap_word_t* gap_buf)
{
    unsigned len = (*gap_buf >> 3);
    std::cout << "[" << *gap_buf << " len=" << len << "] ";
    for (unsigned i = 0; i < len; ++i)
    {
        ++gap_buf;
        std::cout << *gap_buf << "; ";
    } 
    std::cout << std::endl;
}

inline
void PrintDGap(const bm::gap_word_t* gap_buf, unsigned gap_len=0)
{

    unsigned len = gap_len ? gap_len : (*gap_buf >> 3);
    std::cout << "[" " len=" << len << "] ";
    unsigned i = gap_len ? 0 : 1;
    for (; i < len; ++i)
    {
        std::cout << gap_buf[i] << "; ";
    } 
    std::cout << std::endl;
}

inline unsigned int iLog2(unsigned int value)
{
    unsigned int l = 0;
    while( (value >> l) > 1 ) ++l;
    return l;
}

inline
unsigned PrintGammaCode(unsigned value)
{
    unsigned bits = 0;
    // Elias gamma encode
    {
        unsigned l = iLog2(value);
        //cout << "log2=" << l << endl;
        for (unsigned i = 0; i < l; ++i)
        {
            std::cout << 0;
            ++bits;
        }
        std::cout << 1; ++bits;
        for (unsigned i = 0; i < l; ++i)
        {
            if (value & 1 << i) 
                std::cout << 1;
            else
                std::cout << 0;
            ++bits;
        }
    }
    return bits;
}

inline 
void PrintDGapGamma(const bm::gap_word_t* gap_buf, unsigned gap_len=0)
{
    unsigned total = 0;
    unsigned len = gap_len ? gap_len : (*gap_buf >> 3);
    std::cout << "[" " len=" << len << "] ";
    unsigned i = gap_len ? 0 : 1;
    for (; i < len; ++i)
    {
        unsigned v = gap_buf[i];

        unsigned bits = PrintGammaCode(v+1);
        std::cout << "; ";
        total += bits;
    } 
    std::cout << "  gamma_bits=" << total << " src_bits =" << len * 16;
    std::cout << std::endl;

}

/// Read dump file into an STL container (vector of some basic type)
///
/// @return 0 - if reading went well
///
template<class VT>
int read_dump_file(const std::string& fname, VT& data)
{
    typedef typename VT::value_type  value_type;

    std::streampos fsize;
    std::ifstream fin(fname.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good())
    {
        return -1;
    }
    fin.seekg(0, std::ios::end);
    fsize = fin.tellg();
    
    data.resize(fsize/sizeof(value_type));
    if (fsize == 0)
    {
        return 0; // empty input
    }
    fin.seekg(0, std::ios::beg);
    
    fin.read((char*) &data[0], fsize);
    if (!fin.good())
    {
        data.resize(0);
        return -2;
    }
    return 0;
    
}

template<class TBV>
void LoadBVector(const char* fname, TBV& bvector, unsigned* file_size=0)
{
    std::ifstream bv_file (fname, std::ios::in | std::ios::binary);
    if (!bv_file.good())
    {
        std::cout << "Cannot open file: " << fname << std::endl;
        exit(1);
    }
    bv_file.seekg(0, std::ios_base::end);
    unsigned length = (unsigned)bv_file.tellg();
    if (length == 0)
    {
        std::cout << "Empty file:" << fname << std::endl;
        exit(1);
    }
    if (file_size)
        *file_size = length;
        
    bv_file.seekg(0, std::ios::beg);
    
    char* buffer = new char[length];

    bv_file.read(buffer, length);
    
    bm::deserialize(bvector, (unsigned char*)buffer);
    
    delete [] buffer;
}

template<class TBV>
void SaveBVector(const char* fname, const TBV& bvector)
{
    std::ofstream bfile (fname, std::ios::out | std::ios::binary);
    if (!bfile.good())
    {
        std::cout << "Cannot open file: " << fname << std::endl;
        exit(1);
    }
    typename TBV::statistics st1;
    bvector.calc_stat(&st1);

    unsigned char* blob = new unsigned char[st1.max_serialize_mem];
    unsigned blob_size = bm::serialize(bvector, blob);


    bfile.write((char*)blob, blob_size);
    bfile.close();

    delete [] blob;
}

inline
void SaveBlob(const char* name_prefix, unsigned num, const char* ext,
              const unsigned char* blob, unsigned blob_size)
{
    std::strstream fname_str;
    fname_str << name_prefix << "-" << num << ext;
    
    char* fname = fname_str.str();
    std::ofstream bfile (fname, std::ios::out | std::ios::binary);
    if (!bfile.good())
    {
        std::cout << "Cannot open file: " << fname << std::endl;
        exit(1);
    }
    bfile.write((char*)blob, blob_size);
    bfile.close();
}


template<typename V> 
void PrintBinary(V val)
{
    for (unsigned i = 0; i < sizeof(V)*8; i++)
    {
        std::cout << (unsigned)((val >> i) & 1);
        if (i == 15 && (sizeof(V)*8 > 16)) std::cout << "-";
    }
//    cout << " :" << val;
}

inline 
void PrintBits32(unsigned val)
{
    PrintBinary(val);
}

void PrintDistanceMatrix(
   const unsigned distance[bm::set_block_plain_cnt][bm::set_block_plain_cnt])
{
    for (unsigned i = 0; i < bm::set_block_plain_cnt; ++i)
    {
        const unsigned* row = distance[i];
        std::cout << i << ": ";
        for (unsigned j = i; j < bm::set_block_plain_cnt; ++j)
        {
            std::cout << std::setw(4) << std::setfill('0') << row[j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename TM>
void PrintTMatrix(const TM& tmatrix, unsigned cols=0, bool binary = false)
{
    unsigned columns = cols ? cols : tmatrix.cols();
    for (unsigned i = 0; i < tmatrix.rows(); ++i)
    {
        const typename TM::value_type* row = tmatrix.row(i);
        std::cout << i << ": ";
        if (i < 10) std::cout << " ";
        for (unsigned j = 0; j < columns; ++j)
        {
            if (!binary)
            {
                std::cout << std::setw(4) << std::setfill('0') << row[j] << " ";
            }
            else
            {
                PrintBinary(row[j]);
            }
        }
        std::cout << std::endl;
    }
}

/// Binary code string converted to number
/// Bits are expected left to right
///
inline
unsigned BinStrLR(const char* str)
{
    unsigned value = 0;
    unsigned bit_idx = 0;
    for (; *str; ++str)
    {
        switch(*str)
        {
        case '0': 
            ++bit_idx;
            break;
        case '1':
            value |= (1 << bit_idx);
            ++bit_idx;
            break;
        default:
            assert(0);
        }
        if (bit_idx == sizeof(unsigned) * 8) 
            break;
    }    
    return value;
}

template<class BV>
void print_blocks_count(const BV& bv)
{
    const unsigned sz = 128000;
    unsigned* bc_arr = new unsigned[sz];
    for(unsigned x = 0; x < sz; ++x) bc_arr[x] = 0;


    unsigned last_block = bv.count_blocks(bc_arr);
    unsigned sum = 0;

    for (unsigned i = 0; i <= last_block; ++i)
    {
        std::cout << i << ":";

        unsigned j = 0;
        for (; i <= last_block; ++i)
        {
            std::cout << std::setw(5) << std::setfill('0') << bc_arr[i] << " ";
            sum += bc_arr[i];
            if (++j == 10) break;
        }
        std::cout << " | " << sum << std::endl;
    }
    std::cout << "Total=" << sum << std::endl;

    delete [] bc_arr;
}
inline 
void print_bc(unsigned i, unsigned count)
{
    static unsigned sum = 0;
    static unsigned row_idx = 0;
    static unsigned prev = 0;

    if (i == 0) 
    {
        sum = row_idx = 0;
    }
    else
    {
        if (prev +1 < i)
            print_bc(prev+1, 0);
        prev = i;
    }

    if (row_idx == 0)
    {
        std::cout << i << ":";
    }

    std::cout << std::setw(5) << std::setfill('0') << count << " ";
    sum += count;

    ++row_idx;
    if (row_idx == 10)
    {
        row_idx = 0;
        std::cout << " | " << sum << std::endl;
    }
}


template<class BV>
void print_stat(const BV& bv, unsigned blocks = 0)
{
    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();

    bm::id_t count = 0;
    int printed = 0;

    int total_gap_eff = 0;

    if (!blocks)
    {
        blocks = bm::set_total_blocks;
    }

    unsigned nb;
    unsigned nb_prev = 0;
    for (nb = 0; nb < blocks; ++nb)
    {
        const bm::word_t* blk = bman.get_block(nb);
        if (blk == 0)
        {
           continue;
        }

        if (IS_FULL_BLOCK(blk))
        {
           if (bman.is_block_gap(nb)) // gap block
           {
               printf("[Alert!%i]", nb);
               assert(0);
           }
           
           unsigned start = nb; 
           for(unsigned i = nb+1; i < bm::set_total_blocks; ++i, ++nb)
           {
               blk = bman.get_block(nb);
               if (IS_FULL_BLOCK(blk))
               {
                 if (bman.is_block_gap(nb)) // gap block
                 {
                     printf("[Alert!%i]", nb);
                     assert(0);
                     --nb;
                     break;
                 }

               }
               else
               {
                  --nb;
                  break;
               }
           }

           printf("{F.%i:%i}",start, nb);
           ++printed;
        }
        else
        {
            if ((nb-1) != nb_prev)
            {
                printf("..%i..", (int)nb-nb_prev);
            }

            if (bman.is_block_gap(nb)) // gap block
            {
                unsigned bc = bm::gap_bit_count(BMGAP_PTR(blk));
                /*unsigned sum = */bm::gap_control_sum(BMGAP_PTR(blk));
                unsigned level = bm::gap_level(BMGAP_PTR(blk));
                count += bc;
               unsigned len = bm::gap_length(BMGAP_PTR(blk))-1;
               unsigned raw_size=bc*2;
               unsigned cmr_len=len*2;
               int mem_eff = raw_size - cmr_len;
               total_gap_eff += mem_eff;
               
               unsigned i,j;
               bman.get_block_coord(nb, &i, &j);
                printf(" [GAP %i(%i,%i)=%i:%i-L%i(%i)] ", nb, i, j, bc, level, len, mem_eff);
                ++printed;
            }
            else // bitset
            {
                const bm::word_t* blk_end = blk + bm::set_block_size;
                unsigned bc = bm::bit_block_calc_count(blk, blk_end);

                unsigned zw = 0;
                for (unsigned i = 0; i < bm::set_block_size; ++i) 
                {
                    zw += (blk[i] == 0);
                }

                count += bc;
                printf(" (BIT %i=%i[%i]) ", nb, bc, zw);
                ++printed;                
            }
        }
        if (printed == 10)
        {
            printed = 0;
            printf("\n");
        }
        nb_prev = nb;
    } // for nb
    printf("\n");
    printf("gap_efficiency=%i\n", total_gap_eff); 

}



template<class SV>
void print_svector_stat(const SV& svect)
{
    typename SV::statistics st;
    svect.calc_stat(&st);
    
    std::cout << "size = " << svect.size() << std::endl;
    std::cout << "Bit blocks:       " << st.bit_blocks << std::endl;
    std::cout << "Gap blocks:       " << st.gap_blocks << std::endl;
    std::cout << "Max serialize mem:" << st.max_serialize_mem << " "
              << (st.max_serialize_mem / (1024 * 1024)) << "MB" << std::endl;
    std::cout << "Memory used:      " << st.memory_used << " "
              << (st.memory_used / (1024 * 1024))       << "MB" << std::endl;
    
    std::cout << "Projected mem usage for vector<int>:" << sizeof(int) * svect.size() << std::endl;
    std::cout << "Projected mem usage for vector<long long>:" << sizeof(long long) * svect.size() << std::endl;
    
    std::cout << "Plains:" << std::endl;
    
    for (unsigned i = 0; i < svect.plains(); ++i)
    {
        const typename SV::bvector_type* bv_plain = svect.plain(i);
        std::cout << i << ":";
            if (bv_plain == 0)
            {
                std::cout << "NULL";
                bool any_else = false;
                for (unsigned j = i+1; j < svect.plains(); ++j) // look ahead
                {
                    if (svect.plain(i))
                    {
                        any_else = true;
                        break;
                    }
                }
                if (!any_else)
                {
                    break;
                }
            }
            else
            {
                std::cout << bv_plain->count();
            }
            std::cout << std::endl;
        
    } // for i
}

// save sparse_vector dump to disk
//
template<class SV>
int file_save_svector(const SV& sv, const std::string& fname, size_t* sv_blob_size=0)
{
    BM_ASSERT(!fname.empty());
    
    bm::sparse_vector_serial_layout<SV> sv_lay;
    
    BM_DECLARE_TEMP_BLOCK(tb)
    auto res = bm::sparse_vector_serialize(sv, sv_lay, tb);
    if (res != 0)
    {
        return res;
    }

    std::ofstream fout(fname.c_str(), std::ios::binary);
    if (!fout.good())
    {
        return -1;
    }
    const char* buf = (char*)sv_lay.buf();
    fout.write(buf, sv_lay.size());
    if (!fout.good())
    {
        return -1;
    }
    
    fout.close();
    
    if (sv_blob_size)
    {
        *sv_blob_size = sv_lay.size();
    }
    return 0;
}


} // namespace

#include "bmundef.h"


#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif
