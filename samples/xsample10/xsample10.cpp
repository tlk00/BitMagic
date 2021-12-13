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

unsigned bookmark_blocks = 16;
bool is_check = false;

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
bm::chrono_taker<std::ostream>::duration_map_type  timing_map;

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

/// Test-compare two data frames for match in [from..to] ranges
///
static
void check_range(const data_frame& df_r, const data_frame& df,
                 unsigned from, unsigned to)
{
    {
        auto it_r_id = df_r.id.get_const_iterator(from);
        auto it_id = df.id.get_const_iterator(from);
        for (auto i = from; i <= to; ++i)
        {
            if (!it_r_id.valid() || !it_id.valid())
                break;
            const  char* id_r_str = *it_r_id;
            const  char* id_str = *it_id;
            int cmp = strcmp(id_r_str, id_str);
            if (cmp != 0)
            {
                cerr << "Range Id comparison check failed! at:" << i << endl;
                exit(1);
            }
            ++it_r_id;
            ++it_id;
        } // for
    }
    {
        auto it_r_pos = df_r.pos.get_const_iterator(from);
        auto it_pos = df.pos.get_const_iterator(from);
        for (auto i = from; i <= to; ++i)
        {
            if (!it_r_pos.valid() || !it_pos.valid())
                break;
            auto pos = *it_pos;
            auto r_pos = *it_r_pos;
            if (pos != r_pos)
            {
                cerr << "Range Pos comparison check failed! at:" << i << endl;
                exit(1);
            }
            ++it_r_pos;
            ++it_pos;
        } // for
    }
    {
        auto it_r_seq = df_r.seq.get_const_iterator(from);
        auto it_seq = df.seq.get_const_iterator(from);
        for (auto i = from; i <= to; ++i)
        {
            if (!it_r_seq.valid() || !it_seq.valid())
                break;
            const  char* seq_r_str = *it_r_seq;
            const  char* seq_str = *it_seq;
            int cmp = strcmp(seq_r_str, seq_str);
            if (cmp != 0)
            {
                cerr << "Range Seq comparison check failed! at:" << i << endl;
                exit(1);
            }
            ++it_r_seq;
            ++it_seq;
        } // for
    }

}

/// Forward scrolling 1
///
/// The simplest forward scroll, uses no overlap check, so last
/// just decoded range just ignored.
/// Slow but simple.
static
void scroll_benchmark_range(
            const compressed_data_frame& df_cz, unsigned size,
            const data_frame& df,
            bm::sparse_vector_deserializer<str_sv64_type>& sv_deserial_str,
            bm::sparse_vector_deserializer<sparse_vector_u32>& sv_deserial_u32)
{
    data_frame df_r; // range data frame (partial)

    unsigned cnt = 0;
    for (unsigned i = 0; i < size; i+=page_size/3, ++cnt)
    {
        unsigned from(i), to(i+page_size-1); // define range variables

        // deserialize range for each vector in the data frame
        //
        sv_deserial_str.deserialize_range(df_r.id,
                    (const unsigned char*)df_cz.id_buf.data(), from, to);
        sv_deserial_str.deserialize_range(df_r.seq,
                    (const unsigned char*)df_cz.seq_buf.data(), from, to);
        sv_deserial_u32.deserialize_range(df_r.pos,
                    (const unsigned char*)df_cz.pos_buf.data(), from, to);

        // as a payload to simulate real work
        // use test comparison with the orginal data frame
        //
        if (is_check)
            check_range(df_r, df, from, to);
    } // for
}

/// Forward scrolling 2
///
/// Faster but more complex:
/// - explicitly clears the not needed range
/// - desralization always goes to a temp vector and then merged into scrolled
/// range data frame. Merge is a fast operation analog of STL move()
///
static
void scroll_benchmark_clear_range_merge(const compressed_data_frame& df_cz,
                                        unsigned size,
                                        const data_frame& df,
            bm::sparse_vector_deserializer<str_sv64_type>& sv_deserial_str,
            bm::sparse_vector_deserializer<sparse_vector_u32>& sv_deserial_u32)
{
    data_frame df_r;

    unsigned from(0), to(0);
    unsigned i = 0;
    for (; i < size; i+=page_size/3)
    {
        {
        data_frame df_tmp;
        // clear the range evicted out of the scrolling window
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
        sv_deserial_str.deserialize_range(df_tmp.id,
                (const unsigned char*)df_cz.id_buf.data(), prev_to, to);
        sv_deserial_str.deserialize_range(df_tmp.seq,
                (const unsigned char*)df_cz.seq_buf.data(), prev_to, to);
        sv_deserial_u32.deserialize_range(df_tmp.pos,
                (const unsigned char*)df_cz.pos_buf.data(), prev_to, to);

        // merge temp dataframe into the target dataframe
        df_r.id.merge(df_tmp.id);
        df_r.seq.merge(df_tmp.seq);
        df_r.pos.merge(df_tmp.pos);
        }

        if (is_check)
            check_range(df_r, df, from, to);

    } // for
}


/// Forward scrolling 3
///
/// About as fast as 3 but a bit simpler:
/// - desralization always goes to a temp vector and then merged into scrolled
/// range data frame.
/// - uses keep_range() to free the memory evicted out of the scrolling window
///
static
void scroll_benchmark_merge_keep_range(const compressed_data_frame& df_cz,
                       unsigned size,
                       const data_frame& df,
            bm::sparse_vector_deserializer<str_sv64_type>& sv_deserial_str,
            bm::sparse_vector_deserializer<sparse_vector_u32>& sv_deserial_u32)
{
    data_frame df_r;

    unsigned from(0), to(0);
    unsigned i = 0;
    for (; i < size; i+=page_size/3)
    {
        data_frame df_tmp;
        auto prev_to = i ? to+1 : 0;
        from = i; to = i+page_size-1;
        assert(from < prev_to || !i);

        // deserialize into a temp (empty) dataframe
        sv_deserial_str.deserialize_range(df_tmp.id,
                (const unsigned char*)df_cz.id_buf.data(), prev_to, to);
        df_r.id.merge(df_tmp.id);
        df_r.id.keep_range(from, to);

        sv_deserial_str.deserialize_range(df_tmp.seq,
                (const unsigned char*)df_cz.seq_buf.data(), prev_to, to);
        df_r.seq.merge(df_tmp.seq);
        df_r.seq.keep_range(from, to);

        sv_deserial_u32.deserialize_range(df_tmp.pos,
                (const unsigned char*)df_cz.pos_buf.data(), prev_to, to);
        df_r.pos.merge(df_tmp.pos);
        df_r.pos.keep_range(from, to);

        if (is_check)
            check_range(df_r, df, from, to); // check the range

    } // for
}

/// Display help
static void show_help()
{
    std::cerr
        << "BitMagic range deserialization example (c) 2020" << std::endl
        << "-h                    -- print help" << std::endl
        << "-check                -- run test within benchmark " << std::endl
        << "-bookm number         -- bookmark parameter (256, 128, 64, 32, 16...)" << std::endl
      ;
}

/// cmd-line arguments parser
///
static int parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help"))
        {
            show_help(); exit(0);
        }
        if (arg == "-check" || arg == "--check")
        {
            is_check = true;
            continue;
        }
        if (arg == "-bookm" || arg == "--bookm")
        {
            if (i + 1 < argc)
                bookmark_blocks = (unsigned)::atoi(argv[++i]);
            else
            {
                std::cerr << "Error: -bookm requires a number" << std::endl;
                return 1;
            }
            continue;
        }
    } // for i
    return 0;
}

int main(int argc, char *argv[])
{
    try
    {
        data_frame            df;    // sample data frame (columnar)
        compressed_data_frame df_cz; // sample data frame compressed into BLOBs


        auto ret = parse_args(argc, argv);
        if (ret != 0)
        {
            cerr << "cmd-line parse error. " << endl;
            return ret;
        }


        generate_test_dataframe(df, test_size);

        // serialize into compressed data frame
        {
            size_t size_total  = 0;
            bm::sparse_vector_serializer<str_sv64_type> sv_serializer;
            bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer_u32;

            // configure with bookmarks for fast range deserialization
            //
            if (bookmark_blocks)
            {
                sv_serializer.set_bookmarks(true, bookmark_blocks);
                sv_serializer_u32.set_bookmarks(true, bookmark_blocks);
                cout << "  Serialization bookmark at:" << bookmark_blocks << endl;
            }
            if (is_check)
                cout << "  Testing is ON" << endl;

            sv_serializer.enable_xor_compression();
            sv_serializer_u32.enable_xor_compression();

            {
                bm::sparse_vector_serial_layout<str_sv64_type> sv_lay_str;
                {
                    sv_serializer.serialize(df.id, sv_lay_str);
                    const unsigned char* buf = sv_lay_str.buf();
                    auto sz = sv_lay_str.size();
                    size_total += sz;
                    cout << "  ID BLOB size = " << sz << endl;

                    df_cz.id_buf.resize(sz);
                    ::memcpy(df_cz.id_buf.data(), buf, sz);
                }
                {
                    sv_serializer.serialize(df.seq, sv_lay_str);
                    const unsigned char* buf = sv_lay_str.buf();
                    auto sz = sv_lay_str.size();
                    size_total += sz;
                    cout << "  SEQ BLOB size = " << sz << endl;

                    df_cz.seq_buf.resize(sz);
                    ::memcpy(df_cz.seq_buf.data(), buf, sz);
                }
            }
            {
                bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay_u32;
                sv_serializer_u32.serialize(df.pos, sv_lay_u32);
                const unsigned char* buf = sv_lay_u32.buf();
                auto sz = sv_lay_u32.size();
                size_total += sz;
                cout << "  POS BLOB size = " << sz << endl;

                df_cz.pos_buf.resize(sz);
                ::memcpy(df_cz.pos_buf.data(), buf, sz);
            }
            cout << "Total = " <<  size_total << endl;

        }

        // create deserializers, one per sparse vector type
        //
        // deserializers are heavy objects so it is best to re-use objects
        //
        bm::sparse_vector_deserializer<str_sv64_type> sv_deserial_str;
        bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial_u32;


        {
            bm::chrono_taker tt1(cout, "01. Scrolling test range", 1, &timing_map);
            scroll_benchmark_range(df_cz, test_size, df,
                                   sv_deserial_str, sv_deserial_u32);
        }

        {
            bm::chrono_taker tt1(cout, "02. Scrolling test clear/range/merge", 1, &timing_map);
            scroll_benchmark_clear_range_merge(df_cz, test_size, df,
                                           sv_deserial_str, sv_deserial_u32);
        }

        {
            bm::chrono_taker tt1(cout, "03. Scrolling merge/keep_range", 1, &timing_map);
            scroll_benchmark_merge_keep_range(df_cz, test_size, df,
                                           sv_deserial_str, sv_deserial_u32);
        }

        cout << endl;
        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_ops_per_sec);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;

}

