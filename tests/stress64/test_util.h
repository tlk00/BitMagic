/*
Copyright(c) 2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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


/// Load bit-vector using ref syntax
///
template<typename BV, typename VT, typename TOut>
void load_BV_set_ref(TOut& tout, BV& bv, const VT& vect, bool print_stat = true)
{
    for (auto it = vect.begin(); it != vect.end(); ++it)
    {
        auto v = *it;
        bv[v] = true;
    }
    assert(bv.count() == vect.size());
    if (print_stat)
        print_bvector_stat(tout, bv);
}

/// Load bit-vector using ref syntax
///
template<typename BV, typename VT, typename TOut>
void clear_BV_set_ref(TOut& tout, BV& bv, const VT& vect, bool print_stat = true)
{
    for (auto it = vect.begin(); it != vect.end(); ++it)
    {
        auto v = *it;
        bv[v] = false;
    }
    if (print_stat)
        print_bvector_stat(tout, bv);
}



/// CMP bit-vector using ref syntax
///
template<typename BV, typename VT>
void compare_BV_set_ref(const BV& bv, const VT& vect, bool compare_count=true)
{
    for (auto it = vect.begin(); it != vect.end(); ++it)
    {
        auto v = *it;
        bool b = bv[v];
        if (!b)
        {
            cerr << "Error! Vector(ref) comparison failed. v=" << v
                 << endl;
            assert(0);
            exit(1);
        }
    }
    if (compare_count)
    {
        auto count = bv.count();
        if (count != vect.size())
        {
            cerr << "Error! Vector(ref) size cmp failed. vect.size()=" << vect.size()
                << " bv.count()=" << count << endl;
            assert(0);
            exit(1);
        }
    }
}

/// CMP bit-vector using enumerator
///
template<typename BV, typename VT>
void compare_BV(const BV& bv, const VT& vect, bool compare_count = true)
{
    typename BV::enumerator en = bv.first();
    auto prev_id = 0ULL;
    for (auto it = vect.begin(); it != vect.end(); ++it, ++en)
    {
        auto v0 = *it;
        if (!en.valid())
        {
            cerr << "Error! Vector(en) comparison failed. enumerator invalid at value=" << v0
                 << endl;
            assert(0); exit(1);
        }
        auto v1 = *en;
        if (v1 != v0)
        {
            cerr << "Error! Vector(en) comparison failed. v0=" << v0
                 << " v1=" << v1
                 << endl;
            assert(0); exit(1);
        }
        if ((v1 != prev_id) && compare_count)
        {
            auto r = bv.count_range(prev_id, v1);
            if (r != 2ULL)
            {
                cerr << "Error! Vector(en) comparison failed. count_range() = " << r
                     << " [" << prev_id << ", " << v1 << "]"
                     << endl;
                assert(0); exit(1);
            }
            bool b = bv.any_range(prev_id, v1);
            if (r)
            {
                assert(b);
            }
            (void)b;
        }
        prev_id = v1;
    }
    if (compare_count)
    {
        auto count = bv.count();
        if (count != vect.size())
        {
            cerr << "Error! Vector(en) size cmp failed. vect.size()=" << vect.size()
                << " bv.count()=" << count << endl;
            assert(0);
            exit(1);
        }
    }
}


template<class BV>
void IntervalsCheck(const BV& bv)
{
    BV bv_inv(bv);
    bv_inv.invert();

    typename BV::size_type intervals = bm::count_intervals(bv);
    typename BV::size_type intervals_c = 1;

    typename BV::enumerator en1 = bv.get_enumerator(0);
    typename BV::enumerator en2 = bv_inv.get_enumerator(0);

    while (en1.valid())
    {
        typename BV::size_type from = *en1;
        typename BV::size_type to = *en2;
        assert(from != to);

        bool all_one, any_one, is_interval;
        if (to == bm::id_max)
        {
            all_one = bv.is_all_one_range(from, to-1);
            assert(all_one);
            any_one = bv.any_range(from, to - 1);
            assert(any_one);
            break;
        }
        else
        {
            all_one = bv.is_all_one_range(from, to);
            assert(!all_one);
            auto cnt = bv.count_range(from, to);
            any_one = bv.any_range(from, to);
            assert(any_one == (cnt>0));
        }

        if (to == bm::id_max)
        {}
        else
        {
            ++intervals_c;
        }
        if (to < from)
        {
            --from;
            assert(!bv.test(to));
            all_one = bv.is_all_one_range(from, to);
            assert(!all_one);
            typename BV::size_type cnt = bv.count_range(to, from);
            assert(!cnt);
            any_one = bv.any_range(to, from);
            assert(!any_one);
            is_interval = bm::is_interval(bv, to, from);
            assert(!is_interval);

            en2.go_to(from+1);
            if (!en2.valid())
                break;
        }
        else
        {
            --to;
            assert(bv.test(to));
            all_one = bv.is_all_one_range(from, to);
            assert(all_one);
            typename BV::size_type cnt = bv.count_range(from, to);
            assert(cnt == (to - from + 1));
            en1.go_to(to+1);
            any_one = bv.any_range(from, to);
            assert(any_one);
            is_interval = bm::is_interval(bv, from, to);
            assert(is_interval);
        }

    } // while
    if (intervals != intervals_c)
    {
        typename BV::size_type diff;
        diff = std::max(intervals, intervals_c) - std::min(intervals, intervals_c);
        if (diff > 1)
        {
            cerr << "Intervals difference:" << diff << endl;
            assert(0);
            exit(1);
        }
    }
}


template<class BV>
void interval_copy_range(BV& bv, const BV& bv_src,
                         typename BV::size_type from, typename BV::size_type to)
{
    bv.clear();

    if (from > to)
        bm::xor_swap(from, to);

    bm::interval_enumerator<BV> ien(bv_src, from, false);
    while (ien.valid())
    {
        auto st = ien.start();
        assert(st >= from);
        if (st > to)
            break;
        auto end = ien.end();
        if (end > to)
            end = to;

        bv.set_range(st, end);
        if (!ien.advance())
            break;
    } // while
}

template<typename BV>
void IntervalsEnumeratorCheck(const BV& bv)
{
    typename BV::allocator_pool_type pool;

    typename BV::size_type f, l, m;
    auto b = bv.find_range(f, l);
    if (!b)
    {
        assert(bv.count() == 0);
        return;
    }
    m = l - f;
    if (!m)
        m = l;

    bool eq;
    // Full vector
    {
        BV bv2; BV bv2_c;
        typename BV::mem_pool_guard g1(pool, bv2);
        typename BV::mem_pool_guard g2(pool, bv2_c);

        bv2_c.copy_range(bv, 0, bm::id_max-1);

        interval_copy_range(bv2, bv, 0, bm::id_max - 1);
        eq = bv2.equal(bv2_c);
        assert(eq);
    }
    // 0 -> frist
    {
        BV bv2, bv2_c;
        typename BV::mem_pool_guard g1(pool, bv2);
        typename BV::mem_pool_guard g2(pool, bv2_c);

        bv2_c.copy_range(bv, 0, l);

        interval_copy_range(bv2, bv, 0, l);
        eq = bv2.equal(bv2_c);
        assert(eq);
    }
    // [first..last]
    {
        BV bv2, bv2_c;
        typename BV::mem_pool_guard g1(pool, bv2);
        typename BV::mem_pool_guard g2(pool, bv2_c);

        bv2_c.copy_range(bv, f, l);

        interval_copy_range(bv2, bv, f, l);
        eq = bv2.equal(bv2_c);
        assert(eq);
    }
    // [last..]
    {
        BV bv2, bv2_c;
        typename BV::mem_pool_guard g1(pool, bv2);
        typename BV::mem_pool_guard g2(pool, bv2_c);

        bv2_c.copy_range(bv, l, bm::id_max-1);

        interval_copy_range(bv2, bv, l, bm::id_max - 1);
        eq = bv2.equal(bv2_c);
        assert(eq);
    }
    // [mid..last]
    {
        BV bv2, bv2_c;
        typename BV::mem_pool_guard g1(pool, bv2);
        typename BV::mem_pool_guard g2(pool, bv2_c);

        bv2_c.copy_range(bv, m, l);

        interval_copy_range(bv2, bv, m, l);
        eq = bv2.equal(bv2_c);
        assert(eq);
    }
    // [first..mid]
    {
        BV bv2, bv2_c;
        typename BV::mem_pool_guard g1(pool, bv2);
        typename BV::mem_pool_guard g2(pool, bv2_c);

        bv2_c.copy_range(bv, f, m);

        interval_copy_range(bv2, bv, f, m);
        eq = bv2.equal(bv2_c);
        assert(eq);
    }

}



template<class SV, class Vect>
bool CompareSparseVector(const SV& sv, const Vect& vect, bool interval_filled = false)
{
    if (vect.size() != sv.size())
    {
        cerr << "Sparse vector size test failed!" << vect.size() << "!=" << sv.size() << endl;
        return false;
    }
    
    if (sv.is_nullable())
    {
        const typename SV::bvector_type* bv_null = sv.get_null_bvector();
        assert(bv_null);
        auto non_null_cnt = bv_null->count();
        if (vect.size() != non_null_cnt)
        {
            if (!interval_filled)
            {
                cerr << "NULL vector count failed." << non_null_cnt << " size=" << vect.size() << endl;
                assert(0); exit(1);
            }
        }
    }
    
    {
        typename SV::const_iterator it = sv.begin();
        typename SV::const_iterator it_end = sv.end();

        for (unsigned i = 0; i < vect.size(); ++i)
        {
            typename Vect::value_type v1 = vect[i];
            typename SV::value_type v2 = sv[i];
            typename SV::value_type v3 = *it;
            
            int cmp = sv.compare(i, v1);
            assert(cmp == 0);
            if (v1 > 0)
            {
                cmp = sv.compare(i, v1-1);
                assert(cmp > 0);
            }

            if (v1 != v2)
            {
                cerr << "SV discrepancy:" << "sv[" << i << "]=" << v2
                     <<  " vect[" << i << "]=" << v1
                     << endl;
                return false;
            }
            if (v1 != v3)
            {
                cerr << "SV discrepancy:" << "sv[" << i << "]=" << v2
                     <<  " *it" << v3
                     << endl;
                return false;
            }
            assert(it < it_end);
            ++it;
        } // for
        if (it != it_end)
        {
            cerr << "sv const_iterator discrepancy!" << endl;
            assert(0);return false;
        }
    }
    
    // extraction comparison
    {
        std::vector<typename SV::value_type> v1(sv.size());
        std::vector<typename SV::value_type> v1r(sv.size());
        sv.extract(&v1[0], sv.size(), 0);
        sv.extract_range(&v1r[0], sv.size(), 0);
        for (unsigned i = 0; i < sv.size(); ++i)
        {
            if (v1r[i] != v1[i] || v1[i] != vect[i])
            {
                cerr << "TestEqualSparseVectors Extract 1 failed at:" << i
                     << " v1[i]=" << v1[i] << " v1r[i]=" << v1r[i]
                     << endl;
                assert(0);exit(1);
            }
        } // for
    }
    
    // serialization comparison
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::sparse_vector_serial_layout<SV> sv_lay;
    bm::sparse_vector_serialize<SV>(sv, sv_lay, tb);
    SV sv2;
    const unsigned char* buf = sv_lay.buf();
    int res = bm::sparse_vector_deserialize(sv2, buf, tb);
    if (res != 0)
    {
        cerr << "De-Serialization error" << endl;
        assert(0);exit(1);
    }
    if (sv.is_nullable() != sv2.is_nullable())
    {
        cerr << "Serialization comparison of two svectors failed (NULL vector)" << endl;
        assert(0);exit(1);
    }
    const typename SV::bvector_type* bv_null = sv.get_null_bvector();
    const typename SV::bvector_type* bv_null2 = sv.get_null_bvector();
    
    if (bv_null != bv_null2 && (bv_null == 0 || bv_null2 == 0))
    {
        cerr << "Serialization comparison (NUUL vector missing)!" << endl;
        assert(0);exit(1);
    }
    if (bv_null)
    {
        if (bv_null->compare(*bv_null2) != 0)
        {
            cerr << "Serialization comparison of two svectors (NULL vectors unmatch)!" << endl;
            assert(0);exit(1);
        }
    }
    if (!sv.equal(sv2) )
    {
        cerr << "Serialization comparison of two svectors failed" << endl;
        assert(0); exit(1);
    }
    return true;
}


template<typename SV, typename VT>
void load_SV_set_ref(SV* sv, const VT& vect)
{
    for (auto it = vect.begin(); it != vect.end(); ++it)
    {
        auto v = *it;
        sv->set(v, v);
    }
}


template<typename SV, typename VT>
void compare_SV_set_ref(const SV& sv, const VT& vect)
{
    for (size_t i = 0; i != vect.size(); ++i)
    {
        auto v = vect[i];
        auto vv = sv[v];
        if (v != vv)
        {
            std::cerr << "SV compare failed at:" << i
                      << " v=" << v << " sv[]=" << vv << std::endl;
            vv = sv[v];
            assert(v == vv);
            exit(1);
        }
    }
}


template<typename SV, typename VT>
void bulk_load_SV_set_ref(SV* sv, const VT& vect)
{
    assert(vect.size());
    typename SV::back_insert_iterator bi(sv->get_back_inserter());
    auto v_prev = vect[0];
    if (v_prev)
        bi.add_null(v_prev);
    *bi = v_prev;
    for (auto it = vect.begin(); it != vect.end(); ++it)
    {
        auto v = *it;
        if (v == v_prev)
            continue;
        assert(v > v_prev);
        typename SV::size_type diff = v - v_prev;
        if (diff > 1)
        {
            bi.add_null(diff-1);
        }
        *bi = v;
        v_prev = v;
    }
    bi.flush();
}

template<typename CSV, typename SV>
void DetailedCompareSparseVectors(const CSV& csv, const SV& sv)
{
    SV   sv_s(bm::use_null);  // sparse vector decompressed

    // de-compression test
    csv.load_to(sv_s);
    /*
    if (!sv.equal(sv_s))
    {
        cerr << "compressed vector load_to (decompression) failed!" << endl;
        exit(1);
    }
    */


    size_t csv_size = csv.size();
    size_t sv_size = sv.size();
    size_t sv_s_size = sv_s.size();

    const typename SV::bvector_type* bv_null_sv = sv.get_null_bvector();
    const typename SV::bvector_type* bv_null_sv_s = sv_s.get_null_bvector();
    const typename SV::bvector_type* bv_null_csv = csv.get_null_bvector();

    if (csv_size != sv_size || sv_s_size != sv_size)
    {
        assert(bv_null_sv != bv_null_csv);

        auto cnt_sv = bv_null_sv->count();
        auto cnt_sv_s = bv_null_sv_s->count();
        auto cnt_csv = bv_null_csv->count();

        if (cnt_sv != cnt_csv)
        {
            cerr << "Sparse compressed vector comparison failed (size check):"
                << "csv.size()=" << csv_size
                << "sv.size()=" << sv_size
                << "cnt sv = " << cnt_sv
                << "cnt csv = " << cnt_csv
                << endl;
            assert(0); exit(1);
        }
        if (cnt_sv_s != cnt_csv)
        {
            cerr << "Restored Sparse vector comparison failed (size check):"
                << "csv.size()=" << csv_size
                << "sv_s.size()=" << sv_s_size
                << "cnt sv = " << cnt_sv
                << "cnt csv = " << cnt_csv
                << endl;
            assert(0); exit(1);
        }
    }
    if (!sv_size)
        return;
 
    {
        int cmp;
        cmp = bv_null_sv->compare(*bv_null_sv_s);
        assert(cmp == 0);
        cmp = bv_null_sv->compare(*bv_null_csv);
        assert(cmp == 0);
    }

    typename SV::size_type sv_first, sv_last;
    typename SV::size_type sv_s_first, sv_s_last;
    typename SV::size_type csv_first, csv_last;

    bool found_sv = bv_null_sv->find_range(sv_first, sv_last);
    bool found_sv_s = bv_null_sv_s->find_range(sv_s_first, sv_s_last);
    bool found_csv = bv_null_csv->find_range(csv_first, csv_last);

    if (found_sv)
    {
        assert(found_sv_s && found_csv);
    }
    else
    {
        assert(!found_sv_s);
        assert(!found_csv);
        sv_first = sv_last = sv_s_first = sv_s_last = csv_first = csv_last = 0;
    }

    if (csv_first != sv_first)
    {
        cerr << csv_first << "!=" << sv_first << endl;
        assert(0); exit(1);
    }
    if (csv_last != sv_last)
    {
        cerr << sv_s_last << "!=" << sv_last << endl;
        assert(0); exit(1);
    }

    if (sv_s_first != sv_first)
    {
        cerr << sv_s_first << "!=" << sv_first << endl;
        assert(0); exit(1);
    }
    if (sv_s_last != sv_last)
    {
        cerr << sv_s_last << "!=" << sv_last << endl;
        assert(0); exit(1);
    }

    assert(sv_first == sv_s_first && sv_first == csv_first);
    assert(sv_last == sv_s_last && sv_last == csv_last);


    cout << "detailed compare from=" << sv_first << " to=" << sv_size << "..." << flush;
    for (typename SV::size_type i = sv_first; i < sv_size; ++i)
    {
        bool is_null_sv = sv.is_null(i);
        bool is_null_sv_s = sv_s.is_null(i);
        bool is_null_csv = csv.is_null(i);
        if (is_null_sv != is_null_csv || is_null_sv != is_null_sv_s)
        {
            cerr << "Detailed csv check failed (null mismatch) at i=" << i
                << " sv=" << is_null_sv
                << " sv_s=" << is_null_sv_s
                << " csv=" << is_null_csv
                << endl;
            int cmp = bv_null_sv->compare(*bv_null_csv);
            if (cmp != 0)
            {
                cerr << "1. cmp=" << cmp << endl;
                exit(1);
            }
            cmp = bv_null_sv->compare(*bv_null_sv_s);
            if (cmp != 0)
            {
                cerr << "2. cmp=" << cmp << endl;
                exit(1);
            }

            assert(0); exit(1);
        }
        if (!is_null_sv)
        {
            auto v1 = sv[i];
            auto v1_s = sv_s[i];
            auto v2 = csv[i];

            if (v1 != v2 || v1_s != v1)
            {
                cerr << "Detailed csv check failed (value mismatch) at i=" << i
                    << " v1=" << v1
                    << " v1_s=" << v1_s
                    << " v2=" << v2
                    << endl;
                assert(0);  exit(1);
            }
        }
    } // for
    cout << "OK" << endl;


    {
        BM_DECLARE_TEMP_BLOCK(tb)
        bm::sparse_vector_serial_layout<CSV> sv_lay;
        bm::sparse_vector_serialize<CSV>(csv, sv_lay, tb);

        CSV csv1;
        const unsigned char* buf = sv_lay.buf();
        bm::sparse_vector_deserialize(csv1, buf, tb);

        if (!csv.equal(csv1))
        {
            cerr << "Conpressed sparse vector serialization comparison failed!" << endl;
            assert(0); exit(1);
        }
    }

}



template<typename CSV>
void CheckCompressedDecode(const CSV& csv,
    typename CSV::size_type from, typename CSV::size_type size)
{
    std::vector<typename CSV::value_type> vect, vect2, vect_tmp;
    vect.resize(size);
    vect2.resize(size);
    vect_tmp.resize(size);

    typename CSV::size_type sz = csv.decode(&vect[0], from, size);
    typename CSV::size_type sz2 = csv.decode_buf(&vect2[0], &vect_tmp[0], from, size);
    assert(sz == sz2);

    typename CSV::const_iterator it = csv.get_const_iterator(from);
    typename CSV::size_type ex_idx = 0;
    for (typename CSV::size_type i = from; i < from + sz; ++i)
    {
        auto v = csv.get(i);
        auto vx = vect[ex_idx];
        auto vx2 = vect[ex_idx];
        auto vx_it = *it;

        if (v != vx || v != vx2 || v != vx_it)
        {
            cerr << "compressed vector decode mismatch from="
                << from << " idx=" << i
                << " v=" << v << " vx=" << vx
                << endl;
            assert(0);  exit(1);
        }
        ++ex_idx;
        ++it;
    }
}

template<typename CSV>
void DetailedCheckCompressedDecode(const CSV& csv)
{
    auto size = csv.size();
    cout << endl;

    {
        typename CSV::size_type size1 = 100;
        for (typename CSV::size_type i = 0; i < size1; )
        {
            CheckCompressedDecode(csv, i, size);
            if (i % 128 == 0)
                cout << "\r" << i << "/" << size1 << flush;
            i++;
        }
    }
    cout << endl;

    {
        typename CSV::size_type size1 = 100000;
        for (typename CSV::size_type i = 0; i < size1; )
        {
            CheckCompressedDecode(csv, i, size1);
            cout << "\r" << i << "/" << size1 << flush;
            i += (unsigned)rand() % 3;
            size1 -= (unsigned)rand() % 5;
        }
    }
    cout << endl;

    {
        typename CSV::size_type size1 = size;
        for (typename CSV::size_type i = size - size / 2; i < size1; )
        {
            CheckCompressedDecode(csv, i, size1);
            cout << "\r" << i << "/" << size1 << flush;
            i += (1 + i);
        }
    }
    cout << endl;

    for (typename CSV::size_type i = size - size / 2; i < size; )
    {
        CheckCompressedDecode(csv, i, size);
        cout << "\r" << i << "/" << size << flush;
        i += (unsigned)rand() % 25000;
    }
    cout << endl;

    for (typename CSV::size_type i = size - size / 2; i < size; )
    {
        if (size <= i)
            break;
        CheckCompressedDecode(csv, i, size);
        cout << "\r" << i << "/" << size << flush;
        i += (unsigned)rand() % 25000;
        size -= (unsigned)rand() % 25000;;
    }
    cout << endl;

}

template<class SV>
bool TestEqualSparseVectors(const SV& sv1, const SV& sv2, bool detailed = true)
{
    if (sv1.size() != sv2.size())
    {
        cerr << "TestEqualSparseVectors failed incorrect size" << endl;
        exit(1);
    }
    
    if (sv1.is_nullable() == sv2.is_nullable())
    {
        bool b = sv1.equal(sv2);
        if (!b)
        {
            cerr << "sv1.equal(sv2) failed" << endl;
            return b;
        }
        const typename SV::bvector_type* bv_null1 = sv1.get_null_bvector();
        const typename SV::bvector_type* bv_null2 = sv2.get_null_bvector();
        
        if (bv_null1 != bv_null2)
        {
            int r = bv_null1->compare(*bv_null2);
            if (r != 0)
            {
                cerr << "sparse NULL-vectors comparison failed" << endl;
                exit(1);
            }
        }
    }
    else  // NULLable does not match
    {
        detailed = true; // simple check not possible, use slow, detailed
    }
    

    // test non-offset extraction
    //
    {
        std::vector<unsigned> v1(sv1.size());
        std::vector<unsigned> v1r(sv1.size());
        std::vector<unsigned> v1p(sv1.size());
        
        sv1.extract(&v1[0], sv1.size(), 0);
        sv1.extract_range(&v1r[0], sv1.size(), 0);
        sv1.extract_planes(&v1p[0], sv1.size(), 0);
        
        for (typename SV::size_type i = 0; i < sv1.size(); ++i)
        {
            if (v1r[i] != v1[i] || v1p[i] != v1[i])
            {
                cerr << "TestEqualSparseVectors Extract 1 failed at:" << i
                     << " v1[i]=" << v1[i] << " v1r[i]=" << v1r[i] << " v1p[i]=" << v1p[i]
                     << endl;
                exit(1);
            }
        } // for
    }

    // test offset extraction
    //
    {
        std::vector<unsigned> v1(sv1.size());
        std::vector<unsigned> v1r(sv1.size());
        std::vector<unsigned> v1p(sv1.size());
        
        typename SV::size_type pos = sv1.size() / 2;
        
        sv1.extract(&v1[0], sv1.size(), pos);
        sv1.extract_range(&v1r[0], sv1.size(), pos);
        sv1.extract_planes(&v1p[0], sv1.size(), pos);
        
        for (typename SV::size_type i = 0; i < sv1.size(); ++i)
        {
            if (v1r[i] != v1[i] || v1p[i] != v1[i])
            {
                cerr << "TestEqualSparseVectors Extract 1 failed at:" << i
                     << " v1[i]=" << v1[i] << " v1r[i]=" << v1r[i] << " v1p[i]=" << v1p[i]
                     << endl;
                exit(1);
            }
        } // for
    }

    {
        SV svv1(sv1);
        SV svv2(sv2);
        
        bm::null_support is_null = (sv1.is_nullable() == sv2.is_nullable()) ? bm::use_null : bm::no_null;
        
        bool b = svv1.equal(svv2, is_null);
        if (!b)
        {
            cerr << "Equal, copyctor comparison failed" << endl;
            return b;
        }

        svv1.swap(svv2);
        b = svv1.equal(svv2, is_null);
        if (!b)
        {
            cerr << "Equal, copyctor-swap comparison failed" << endl;
            return b;
        }
    }

    // comparison using elements assignment via reference
    if (detailed)
    {
        SV sv3;
        sv3.resize(sv1.size());
        for (typename SV::size_type i = 0; i < sv1.size(); ++i)
        {
            sv3[i] = sv1[i];
            unsigned v1 = sv1[i];
            unsigned v2 = sv3[i];
            if (v1 != v2)
            {
                cerr << "1. sparse_vector reference assignment validation failed" << endl;
                return false;
            }
        }
        bm::null_support is_null = (sv1.is_nullable() == sv3.is_nullable()) ? bm::use_null : bm::no_null;
        bool b = sv1.equal(sv3, is_null);
        if (!b)
        {
            cerr << "2. sparse_vector reference assignment validation failed" << endl;
            return b;
        }
    }
    
    // comparison via const_iterators
    //
    {{
        typename SV::const_iterator it1 = sv1.begin();
        typename SV::const_iterator it2 = sv2.begin();
        typename SV::const_iterator it1_end = sv1.end();
        
        for (; it1 < it1_end; ++it1, ++it2)
        {
            if (*it1 != *it2)
            {
                cerr << "1. sparse_vector::const_iterator validation failed" << endl;
                return false;
            }
        }
    }}

    // comparison through serialization
    //
    {{
        int res;
        bm::sparse_vector_serial_layout<SV> sv_lay;
        bm::sparse_vector_serialize(sv1, sv_lay);
        
        // copy buffer to check if serialization size is actually correct
        const unsigned char* buf = sv_lay.buf();
        size_t buf_size = sv_lay.size();
        
        vector<unsigned char> tmp_buf(buf_size);
        ::memcpy(&tmp_buf[0], buf, buf_size);
        
        SV sv3;
        res = bm::sparse_vector_deserialize(sv3, &tmp_buf[0]);
        if (res != 0)
        {
            cerr << "De-Serialization error in TestEqualSparseVectors()" << endl;
            exit(1);
        }
        
        const typename SV::bvector_type* bv_null1 = sv1.get_null_bvector();
        const typename SV::bvector_type* bv_null2 = sv2.get_null_bvector();
        const typename SV::bvector_type* bv_null3 = sv3.get_null_bvector();
        
        if (bv_null1 && bv_null3)
        {
            int r = bv_null1->compare(*bv_null3);
            if (r != 0)
            {
                cerr << "2. NULL bvectors comparison failed" << endl;
                exit(1);
            }
        }
        if (bv_null1 && bv_null2)
        {
            int r = bv_null1->compare(*bv_null2);
            if (r != 0)
            {
                cerr << "3. NULL bvectors comparison failed" << endl;
                exit(1);
            }
        }

        bm::null_support is_null = (sv1.is_nullable() == sv3.is_nullable()) ? bm::use_null : bm::no_null;
        if (!sv1.equal(sv3, is_null) )
        {
            cerr << "Serialization comparison of two svectors failed (1)" << endl;
            exit(1);
        }
        is_null = (sv2.is_nullable() == sv3.is_nullable()) ? bm::use_null : bm::no_null;
        if (!sv2.equal(sv3, is_null))
        {
            cerr << "Serialization comparison of two svectors failed (2)" << endl;
            exit(1);
        }
        
    
    }}
    return true;
}



template<typename SSV>
void CompareStrSparseVector(const SSV& str_sv,
                            const std::vector<string>& str_coll)
{
    assert(str_sv.size() == str_coll.size());
    
    string str_h = "z";
    string str_l = "A";
    
    typedef typename SSV::bvector_type bvect;

    bm::sparse_vector_scanner<bm::str_sparse_vector<char, bvect, 32> > scanner;

    typename SSV::const_iterator it = str_sv.begin();
    string str;
    for (typename SSV::size_type i = 0; i < str_sv.size(); ++i, ++it)
    {
        assert (it.valid());
        assert (it != str_sv.end());
        
        str_sv.get(i, str);
        const string& str_control = str_coll[i];
        if (str != str_control)
        {
            std::cerr << "String mis-match at:" << i << std::endl;
            exit(1);
        }
        {
            const char* s = *it;
            int cmp = ::strcmp(s, str_control.c_str());
            if (cmp != 0)
            {
                cerr << "Iterator comparison failed! " << s << " != " << str_control
                     << endl;
                exit(1);
            }
            typename SSV::const_iterator it2 = str_sv.get_const_iterator(i);
            assert(it == it2);
            s = *it2;
            cmp = ::strcmp(s, str_control.c_str());
            if (cmp != 0)
            {
                cerr << "2. Iterator comparison failed! " << s << " != " << str_control
                     << endl;
                exit(1);
            }
        }
        int cmp = str_sv.compare(i, str_control.c_str());
        if (cmp != 0)
        {
            std::cerr << "String comparison failure at:" << i << std::endl;
            exit(1);
        }
        if (!str_sv.is_remap()) // re-mapped vectors can give incorrect compare
        {
            cmp = str_sv.compare(i, str_h.c_str());
            if (cmp < 0)
            {
                assert(str < str_h);
            }
            if (cmp > 0)
            {
                assert(str > str_h);
            }

            cmp = str_sv.compare(i, str_l.c_str());
            if (cmp < 0)
            {
                assert(str < str_l);
            }
            if (cmp > 0)
            {
                assert(str > str_l);
            }
        }
        
       typename SSV::size_type pos;
       bool found = scanner.find_eq_str(str_sv, str_control.c_str(), pos);
       if (!found)
       {
            cerr << "Scanner search failed! " << str_control << endl;
            exit(1);
       }
       assert(pos == i);
        if (i % 100000 == 0)
        {
            cout << "\r" << i << " / " << str_sv.size() << flush;
        }
    } // for
    cout << endl;
}


template<typename BV>
void TestFindDiff(const BV& bv1, BV& bv2)
{
    bool f;
    typename BV::size_type pos, pos_c, pos_l;
    f = bv1.find_first_mismatch(bv2, pos);
    BV bv_x;
    bv_x.bit_xor(bv1, bv2, BV::opt_compress);
    if (!f)
    {
        auto a = bv_x.any();
        assert(!a);
        return;
    }
    else // found
    {
        bool f2 = bv1.find_first_mismatch(bv2, pos_l, pos);
        assert(f2 == f);
        assert(pos_l == pos);
        if (pos)
        {
            f2 = bv1.find_first_mismatch(bv2, pos_l, pos-1);
            assert(!f2);
        }
    }

    bool cf = bv_x.find(pos_c);
    assert(f == cf);
    assert(pos == pos_c);

    f = bv2.find_first_mismatch(bv1, pos);
    assert(f == cf);
    assert(pos == pos_c);
}
