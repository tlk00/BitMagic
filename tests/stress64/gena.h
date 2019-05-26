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



/// generate certain simple benchmark values
///
template<typename VT>
void generate_vect_simpl0(VT& vect)
{
    VT v_tmp {0, 10, 31, 32, 62, 63,
             (5 * bm::bits_in_array), (5 * bm::bits_in_array)+1,
             bm::id_max32-1, bm::id_max32, bm::id64_t(bm::id_max32)+1,
             bm::id_max48-1
            };
    std::swap(vect, v_tmp);
}


// generate pseudo-random bit-vector, mix of blocks
//
template<typename BV>
void generate_bvector(BV& bv, typename BV::size_type vector_max, bool optimize)
{
    typename BV::size_type i, j;
    for (i = 0; i < vector_max;)
    {
        // generate bit-blocks
        for (j = 0; j < 65535*8; i += 10, j++)
        {
            bv.set(i);
        }
        if (i > vector_max)
            break;
        // generate GAP (compressed) blocks
        for (j = 0; j < 65535; i += 120, j++)
        {
            unsigned len = rand() % 64;
            bv.set_range(i, i + len);
            i += len;
            if (i > vector_max)
                break;
        }
    }
    if (optimize)
        bv.optimize();
}
