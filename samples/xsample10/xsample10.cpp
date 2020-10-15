/*
Copyright(c) 2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example xsample10.cpp

    Example to show case scrolling decompression with range deserailization.

    @sa bm::sparse_vector
    @sa bm::str_sparse_vector
    @sa bm::sparse_vector_deserializer
    @sa bm::sparse_vector_deserializer::deserialize_range
    @sa bm::sparse_vector_serializer
    @sa bm::sparse_vector_serializer::set_bookmarks
    @sa bm::sparse_vector_serializer::enable_xor_compression
    @sa bm::str_sparse_vector::clear_range
    @sa bm::str_sparse_vector::keep_range
    @sa bm::sparse_vector::clear_range
    @sa bm::sparse_vector::keep_range
*/

/*! \file xsample10.cpp
    \brief Example: scrolling decompression with range deserailization
*/


#include <iostream>
#include <memory>
#include <vector>
#include <stdexcept>

using namespace std;

#include "bm.h"
#include "bmtimer.h"
#include "bmsparsevec.h"
#include "bmstrsparsevec.h"


#include "bmdbg.h"

#include "bmundef.h" /* clear the pre-proc defines from BM */

// ----------------------------------------------------
// Global parameters and types
// ----------------------------------------------------

const unsigned test_size = 25000000;  // number of records to generate
const unsigned page_size = 1000000;


typedef bm::bvector<>                                       bvector_type;
typedef bvector_type::size_type                             bv_size_type;
typedef bm::sparse_vector<unsigned, bvector_type>           sparse_vector_u32;
typedef bm::str_sparse_vector<char, bvector_type, 64>       str_sv64_type;

/// Test data frame (does not represent any particular use case)
///
struct data_frame
{
    str_sv64_type        id;
    sparse_vector_u32    pos;
    str_sv64_type        seq;
};

struct compressed_data_frame
{
    std::vector<char> id_buf;
    std::vector<char> pos_buf;
    std::vector<char> seq_buf;
};



typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32>  rsc_sparse_vector_u32;
typedef std::vector<std::pair<bv_size_type, bv_size_type> > bv_ranges_vector;



// timing storage for benchmarking
bm::chrono_taker::duration_map_type  timing_map;

/// Generate random test data, resembles bio-data but does not mean anything
/// in particular
///
static
void generate_test_dataframe(data_frame& df, unsigned size)
{
    {
        auto bit_id  = df.id.get_back_inserter();
        auto bit_pos = df.pos.get_back_inserter();
        auto bit_seq = df.seq.get_back_inserter();

        string id_str;
        for (unsigned i = 0; i < size; ++i)
        {
            unsigned id = (unsigned)rand();
            string prefix;
            switch (id%16)
            {
                case 1: prefix = "csv"; break;
                case 2: prefix = "snv"; break;
                default: prefix = "rs"; break;
            } // switch
            id_str = prefix;
            id_str.append(std::to_string(id));

            bit_id = id_str;
            bit_pos = i;
            const char* s = "TATATA";
            switch (i%16)
            {
                case 0: s = "ATGCC"; break;
                case 1: s = "GCTG"; break;
                case 3: s = "TAGC"; break;
                case 4: s = "GCTAA"; break;
                case 5: s = "CGATA"; break;
                case 6: s = "GCGGTA"; break;
                default: break;
            } // switch
            bit_seq = s;
        } // for i

        bit_id.flush();
        bit_pos.flush();
        bit_seq.flush();
    }

    BM_DECLARE_TEMP_BLOCK(tb)
    df.id.remap();
    df.id.optimize(tb);

    df.pos.optimize(tb);

    df.seq.remap();
    df.seq.optimize(tb);
}

static
void check_range(const data_frame& df_r, const data_frame& df,
                 unsigned from, unsigned to)
{
    auto it_r_id = df_r.id.get_const_iterator(from);
    auto it_id = df.id.get_const_iterator(from);

    for (auto i = from; i <= to; ++i)
    {
        if (!it_r_id.valid() || !it_id.valid())
        {
            break;
        }
        const  char* id_r_str = *it_r_id;
        const  char* id_str = *it_id;
        int cmp = strcmp(id_r_str, id_str);
        if (cmp != 0)
        {
            cerr << "Range comparison check failed! at:" << i << endl;
            exit(1);
        }
    } // for

}

static
void scroll_benchmark1(const compressed_data_frame& df_cz, unsigned size,
                       const data_frame& df)
{
    data_frame df_r; // range data frame (partial)

    bm::sparse_vector_deserializer<str_sv64_type> sv_deserial_str;
    bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial_u32;


    unsigned cnt = 0;
    for (unsigned i = 0; i < size; i+=page_size/3, ++cnt)
    {
        unsigned from(i), to(i+page_size-1);
        sv_deserial_str.deserialize_range(df_r.id, (const unsigned char*)df_cz.id_buf.data(), from, to);
        sv_deserial_str.deserialize_range(df_r.seq, (const unsigned char*)df_cz.seq_buf.data(), from, to);
        sv_deserial_u32.deserialize_range(df_r.pos, (const unsigned char*)df_cz.pos_buf.data(), from, to);

//cout << from << ".." << to << endl;
        check_range(df_r, df, from, to);
    } // for
//    cout << "Total=" << cnt << endl;
}

static
void scroll_benchmark2(const compressed_data_frame& df_cz, unsigned size,
                        const data_frame& df)
{
    data_frame df_r;

    bm::sparse_vector_deserializer<str_sv64_type> sv_deserial_str;
    bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial_u32;


    unsigned from(0), to(0);
    unsigned i = 0;
    for (; i < size; i+=page_size/3)
    {
        data_frame df_tmp;
        // clear previous range
        if (i)
        {
            df_r.id.clear_range(from, i-1);
            df_r.seq.clear_range(from, i-1);
            df_r.pos.clear_range(from, i-1);
        }

        auto prev_to = i ? to+1 : 0;
        from = i; to = i+page_size-1;
        assert(from < prev_to || !i);

        // deserialize into a temp (empty) dataframe
        sv_deserial_str.deserialize_range(df_tmp.id, (const unsigned char*)df_cz.id_buf.data(), prev_to, to);
        sv_deserial_str.deserialize_range(df_tmp.seq, (const unsigned char*)df_cz.seq_buf.data(), prev_to, to);
        sv_deserial_u32.deserialize_range(df_tmp.pos, (const unsigned char*)df_cz.pos_buf.data(), prev_to, to);

        // merge temp dataframe
        df_r.id.merge(df_tmp.id);
        df_r.seq.merge(df_tmp.seq);
        df_r.pos.merge(df_tmp.pos);

        check_range(df_r, df, from, to);

    } // for
}


static
void scroll_benchmark3(const compressed_data_frame& df_cz, unsigned size,
                        const data_frame& df)
{
    data_frame df_r;

    bm::sparse_vector_deserializer<str_sv64_type> sv_deserial_str;
    bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial_u32;


    unsigned from(0), to(0);
    unsigned i = 0;
    for (; i < size; i+=page_size/3)
    {
        data_frame df_tmp;
        auto prev_to = i ? to+1 : 0;
        from = i; to = i+page_size-1;
        assert(from < prev_to || !i);

        // deserialize into a temp (empty) dataframe
        sv_deserial_str.deserialize_range(df_tmp.id, (const unsigned char*)df_cz.id_buf.data(), prev_to, to);
        df_r.id.merge(df_tmp.id);
        df_r.id.keep_range(from, to);

        sv_deserial_str.deserialize_range(df_tmp.seq, (const unsigned char*)df_cz.seq_buf.data(), prev_to, to);
        df_r.seq.merge(df_tmp.seq);
        df_r.seq.keep_range(from, to);

        sv_deserial_u32.deserialize_range(df_tmp.pos, (const unsigned char*)df_cz.pos_buf.data(), prev_to, to);
        df_r.pos.merge(df_tmp.pos);
        df_r.pos.keep_range(from, to);

        check_range(df_r, df, from, to);

    } // for
}



int main(void)
{
    try
    {
        data_frame df;
        compressed_data_frame df_cz;

        generate_test_dataframe(df, test_size);

        // serialize into compressed data frame
        {
            bm::sparse_vector_serializer<str_sv64_type> sv_serializer;
            bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer_u32;

            // configure with bookmarks for fast range deserialization
            //
            sv_serializer.set_bookmarks(true, 32);
            sv_serializer.enable_xor_compression();
            sv_serializer_u32.set_bookmarks(true, 32);
            sv_serializer_u32.enable_xor_compression();

            {
                bm::sparse_vector_serial_layout<str_sv64_type> sv_lay_str;
                {
                    sv_serializer.serialize(df.id, sv_lay_str);
                    const unsigned char* buf = sv_lay_str.buf();
                    auto sz = sv_lay_str.size();
                    cout << sz << endl;

                    df_cz.id_buf.resize(sz);
                    ::memcpy(df_cz.id_buf.data(), buf, sz);
                }

                {
                    sv_serializer.serialize(df.seq, sv_lay_str);
                    const unsigned char* buf = sv_lay_str.buf();
                    auto sz = sv_lay_str.size();
                    cout << sz << endl;

                    df_cz.seq_buf.resize(sz);
                    ::memcpy(df_cz.seq_buf.data(), buf, sz);
                }
            }

            {
                bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay_u32;
                sv_serializer_u32.serialize(df.pos, sv_lay_u32);
                const unsigned char* buf = sv_lay_u32.buf();
                auto sz = sv_lay_u32.size();
                    cout << sz << endl;

                df_cz.pos_buf.resize(sz);
                ::memcpy(df_cz.pos_buf.data(), buf, sz);
            }
        }

        {
            bm::chrono_taker tt1("01. Scrolling test 1", 1, &timing_map);
            scroll_benchmark1(df_cz, test_size, df);
        }

        {
            bm::chrono_taker tt1("02. Scrolling test 2", 1, &timing_map);
            scroll_benchmark2(df_cz, test_size, df);
        }

        {
            bm::chrono_taker tt1("02. Scrolling test 3", 1, &timing_map);
            scroll_benchmark3(df_cz, test_size, df);
        }

        cout << endl;
        bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_ops_per_sec);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;

}
