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
template<typename BV, typename VT>
void load_BV_set_ref(BV& bv, const VT& vect, bool print_stat = true)
{
    for (auto it = vect.begin(); it != vect.end(); ++it)
    {
        auto v = *it;
        bv[v] = true;
    }
    assert(bv.count() == vect.size());
    if (print_stat)
        print_bvector_stat(bv);
}

/// Load bit-vector using ref syntax
///
template<typename BV, typename VT>
void clear_BV_set_ref(BV& bv, const VT& vect, bool print_stat = true)
{
    for (auto it = vect.begin(); it != vect.end(); ++it)
    {
        auto v = *it;
        bv[v] = false;
    }
    if (print_stat)
        print_bvector_stat(bv);
}



/// CMP bit-vector using ref syntax
///
template<typename BV, typename VT>
void compare_BV_set_ref(const BV& bv, const VT& vect)
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
    auto count = bv.count();
    if (count != vect.size())
    {
        cerr << "Error! Vector(ref) size cmp failed. vect.size()=" << vect.size()
             << " bv.count()=" << count << endl;
        assert(0);
        exit(1);
    }
}

/// CMP bit-vector using enumerator
///
template<typename BV, typename VT>
void compare_BV(const BV& bv, const VT& vect)
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
        if ((v1 != prev_id))
        {
            auto r = bv.count_range(prev_id, v1);
            if (r != 2ULL)
            {
                cerr << "Error! Vector(en) comparison failed. count_range() = " << r
                     << " [" << prev_id << ", " << v1 << "]"
                     << endl;
                assert(0); exit(1);
            }
        }
        prev_id = v1;
    }
    auto count = bv.count();
    if (count != vect.size())
    {
        cerr << "Error! Vector(en) size cmp failed. vect.size()=" << vect.size()
             << " bv.count()=" << count << endl;
        assert(0);
        exit(1);
    }
}
