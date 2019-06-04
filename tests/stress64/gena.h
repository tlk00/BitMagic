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

template<typename VT, typename SIZE_TYPE>
void generate_test_vectors(VT &v1,
                           VT &v2,
                           VT &v3,
                           SIZE_TYPE from,
                           SIZE_TYPE to)
{
    SIZE_TYPE j;
    for (j = from; j < to; j += 2)
        v1.push_back(j);
    for (j = from; j < to; j += 5)
        v2.push_back(j);
    for (j = from; j < to; j += 120)
        v3.push_back(j);
}


template<typename BV>
void SimpleGapFillSets(BV&   bv0,
                       BV&   bv1,
                       typename BV::size_type min,
                       typename BV::size_type max,
                       unsigned fill_factor)
{
    typename BV::bulk_insert_iterator bii1(bv1);
    for (typename BV::size_type i = min; i < max; i += fill_factor)
    {
        bv0.set(i);
        bii1 = i;
    } // for i
}

//
// Interval filling.
// 111........111111........111111..........11111111.......1111111...
//

template<typename BVMINI, typename BV, typename SZT>
void FillSetsIntervals(BVMINI* bvect_min,
    BV& bvect_full,
    SZT min,
    SZT max,
    SZT fill_factor,
    bool set_flag = true)
{

    while (fill_factor == 0)
    {
        fill_factor = rand() % 10;
    }
    bvect_full.init();

    cout << "Intervals filling. Factor="
        << fill_factor << endl << endl;

    SZT i, j;
    SZT factor = 70 * fill_factor;
    for (i = min; i < max; ++i)
    {
        unsigned len;
        SZT end;

        do
        {
            len = unsigned(rand()) % factor;
            end = i + len;

        } while (end >= max);
        if (i < end)
        {
            bvect_full.set_range(i, end - 1, set_flag);
        }

        for (j = i; j < end; ++j)
        {
            if (set_flag)
            {
                if (bvect_min)
                    bvect_min->set_bit(j);
            }
            else
            {
                if (bvect_min)
                    bvect_min->clear_bit(j);
            }
        } // j
        i = end;
        len = unsigned(rand()) % (factor * 10 * bm::gap_max_bits);
        if (len % 2)
        {
            len *= unsigned(rand()) % (factor * 10);
        }

        i += len;
        if ((len % 6) == 0)
        {
            for (unsigned k = 0; k < 1000 && i < max; k += 3, i += 3)
            {
                if (set_flag)
                {
                    if (bvect_min)
                        bvect_min->set_bit(i);
                    bvect_full.set_bit_no_check(i);
                }
                else
                {
                    if (bvect_min)
                        bvect_min->clear_bit(j);
                    bvect_full.clear_bit(j);
                }
            }
        }
    } // for i
}

template<typename BV, typename SZT>
void FillSetsIntervals(
    BV& bvect_full,
    SZT min,
    SZT max,
    SZT fill_factor,
    bool set_flag = true)
{
    while (fill_factor == 0)
    {
        fill_factor = rand() % 10;
    }
    bvect_full.init();

    cout << "Intervals filling. Factor="
        << fill_factor << endl << endl;

    SZT i;
    SZT factor = 70 * fill_factor;
    for (i = min; i < max; ++i)
    {
        unsigned len;
        SZT end;

        do
        {
            len = unsigned(rand()) % factor;
            end = i + len;

        } while (end >= max);
        if (i < end)
        {
            bvect_full.set_range(i, end - 1, set_flag);
        }

        i = end;

        len = unsigned(rand()) % (factor * 10 * bm::gap_max_bits);
        if (len % 2)
        {
            len *= unsigned(rand()) % (factor * 10);
        }

        i += len;

        if ((len % 6) == 0)
        {
            for (unsigned k = 0; k < 1000 && i < max; k += 3, i += 3)
            {
                if (set_flag)
                {
                    bvect_full.set_bit_no_check(i);
                }
                else
                {
                    bvect_full.clear_bit(i);
                }
            }
        }
    } // for i
}


template<typename SZT>
SZT random_minmax(SZT min, SZT max)
{
    SZT r = (unsigned(rand()) << 16u) | unsigned(rand());
    if (sizeof(SZT) == 8)
    {
        SZT r2 = (unsigned(rand()) << 16u) | unsigned(rand());
        r |= (r2 << 32);
    }
    return r % (max - min) + min;
}

template <typename BVMINI, typename BV, typename SZT>
void FillSets(BVMINI* bvect_min,
              BV*     bvect_full,
              SZT     min,
              SZT     max,
              SZT     fill_factor)
{
    SZT i;
    SZT id;

    //Random filling
    if (fill_factor == 0)
    {
        SZT n_id = (max - min) / 100;
        cout << "random filling : " << n_id << endl;
        for (i = 0; i < n_id; i++)
        {
            id = random_minmax(min, max);
            bvect_min->set_bit(id);
            bvect_full->set_bit(id);
        }
        cout << endl;
    }
    else
    {
        cout << "fill_factor random filling : factor = "
             << fill_factor << std::endl;

        for (i = 0; i < fill_factor; i++)
        {
            unsigned k = unsigned(rand()) % 10;
            if (k == 0)
                k += 2;

            //Calculate start
            SZT start = min + (max - min) / (fill_factor * k);

            //Randomize start
            start += random_minmax(1ULL, (max - min) / (fill_factor * 10));

            if (start > max)
            {
                start = min;
            }

            //Calculate end 
            SZT end = start + (max - start) / (fill_factor * 2);

            //Randomize end
            end -= random_minmax(1ULL, (max - start) / (fill_factor * 10));

            if (end > max)
            {
                end = max;
            }

            typename BV::bulk_insert_iterator iit = bvect_full->inserter();
            if (fill_factor > 1)
            {
                for (; start < end;)
                {
                    unsigned r = unsigned(rand()) % 8;

                    if (r > 7)
                    {
                        unsigned inc = unsigned(rand()) % 3;
                        ++inc;
                        SZT end2 = start + rand() % 1000;
                        if (end2 > end)
                            end2 = end;
                        while (start < end2)
                        {
                            bvect_min->set_bit(start);
                            iit = start;
                            start += inc;
                        }
                        continue;
                    }

                    if (r)
                    {
                        bvect_min->set_bit(start);
                        iit = start;
                        ++start;
                    }
                    else
                    {
                        start += r;
                        bvect_min->set_bit(start);
                        iit = start;
                    }
                }
            }
            else
            {
                unsigned c = unsigned(rand()) % 15;
                if (c == 0)
                    ++c;
                for (; start < end; ++start)
                {
                    bvect_min->set_bit(start);
                    iit = start;
                    if (start % c)
                    {
                        start += c;
                    }
                }
            }
            cout << endl;
        }
    }
}

template <typename BVMINI, typename BV, typename SZT>
void FillSetClearIntervals(BVMINI* bvect_min,
                           BV* bvect_full,
                           SZT min,
                           SZT max,
                           SZT fill_factor)
{
    FillSetsIntervals(bvect_min, *bvect_full, min, max, fill_factor, true);
    FillSetsIntervals(bvect_min, *bvect_full, min, max, fill_factor, false);
}

template <typename BVMINI, typename BV, typename SZT>
void FillSetsRandomOne(BVMINI* bvect_min,
                       BV* bvect_full,
                       SZT min,
                       SZT max)
{
    SZT range = max - min;
    SZT bit_idx = SZT(rand()) % range;
    bvect_min->set_bit(bit_idx);
    bvect_full->set_bit(bit_idx);
    cout << "Bit_idx=" << bit_idx << endl;
}

template <typename BVMINI, typename BV, typename SZT>
void FillSetsRandom(BVMINI* bvect_min,
                    BV* bvect_full,
                    SZT min,
                    SZT max,
                    SZT fill_factor)
{
    bvect_full->init();
    SZT diap = max - min;
    SZT count;

    switch (fill_factor)
    {
    case 0:
        count = diap / 1000;
        break;
    case 1:
        count = diap / 100;
        break;
    default:
        count = diap / 10;
        break;

    }

    for (unsigned i = 0; i < count; ++i)
    {
        SZT bn = SZT(rand()) % count;
        bn += min;
        if (bn > max)
        {
            bn = max;
        }
        bvect_min->set_bit(bn);
        bvect_full->set_bit_no_check(bn);
    }
    cout << "Ok" << endl;

}

template <typename BVMINI, typename BV, typename SZT>
void FillSetsRegular(BVMINI* bvect_min,
                     BV* bvect_full,
              SZT /*min*/,
              SZT max,
              SZT /*fill_factor*/)
{
    typename BV::bulk_insert_iterator iit = bvect_full->inserter();
    SZT step = rand() % 4;
    if (step < 2) ++step;
    for (SZT i = 0; i < max; i+=step)
    {
        bvect_min->set_bit(i);
        iit = i;
    }
    cout << "Ok" << endl;
}



//
//  Quasi random filling with choosing randomizing method.
//
//
template <typename BVMINI, typename BV, typename SZT>
void FillSetsRandomMethod(BVMINI* bvect_min,
                          BV* bvect_full,
                          SZT min,
                          SZT max,
                          int optimize = 0,
                          int method = -1)
{
    if (method == -1)
    {
        method = rand() % 7;
    }
    SZT factor;
///method = 3;
    switch (method)
    {

    case 0:
        cout << "Random filling: method - FillSets - factor(0)" << endl;
        FillSets(bvect_min, bvect_full, min, max, 0ull);
        break;

    case 1:
        cout << "Random filling: method - FillSets - factor(random)" << endl;
        factor = rand()%3;
        FillSets(bvect_min, bvect_full, min, max, factor?factor:1);
        break;

    case 2:
        cout << "Random filling: method - Set-Clear Intervals - factor(random)" << endl;
        factor = rand()%10;
        FillSetClearIntervals(bvect_min, bvect_full, min, max, factor);
        break;
    case 3:
        cout << "Random filling: method - FillRandom - factor(random)" << endl;
        factor = rand()%3;
        FillSetsRandom(bvect_min, bvect_full, min, max, factor?factor:1);
        break;
    case 4:
        cout << "Random set one bit" << endl;
        FillSetsRandomOne(bvect_min, bvect_full, min, max);
        break;
    case 5:
        cout << "Regular pattern filling" << endl;
        FillSetsRegular(bvect_min, bvect_full, min, max, 2ull);
        break;
    default:
        cout << "Random filling: method - Set Intervals - factor(random)" << endl;
        factor = rand()%10;
        FillSetsIntervals(bvect_min, *bvect_full, min, max, factor);
        break;

    } // switch

    if (optimize && (method <= 1))
    {
        cout << "Vector optimization..." << flush;
        BM_DECLARE_TEMP_BLOCK(tb)
        bvect_full->optimize(tb);
        cout << "OK" << endl;
    }
}



template<typename BV>
void generate_sparse_bvector(BV& bv,
                             typename BV::size_type min,
                             typename BV::size_type max = 40000000,
                             unsigned fill_factor = 65536)
{
    typename BV::bulk_insert_iterator iit(bv);
    unsigned ff = fill_factor / 10;
    for (typename BV::size_type i = min; i < max; i+= ff)
    {
        iit = i;
        ff += ff / 2;
        if (ff > fill_factor)
            ff = fill_factor / 10;
    }
    iit.flush();
}


template<typename VECT>
void GenerateShiftTestCollection(VECT* target,
                            unsigned count,
                            unsigned long long vector_max,
                            bool optimize)
{
    assert(target);
    typename VECT::value_type bv_common; // sub-vector common for all collection
    generate_sparse_bvector(bv_common, vector_max/10, vector_max, 250000);
    
    unsigned cnt1 = (count / 2);
    unsigned i = 0;
    
    for (i = 0; i < cnt1; ++i)
    {
        std::unique_ptr<typename VECT::value_type> bv (new typename VECT::value_type);
        generate_bvector(*bv, vector_max, optimize);
        *bv |= bv_common;
        if (optimize)
            bv->optimize();
        target->push_back(std::move(*bv));
    } // for
    
    unsigned long long fill_factor = 10;
    for (; i < count; ++i)
    {
        std::unique_ptr<typename VECT::value_type> bv (new typename VECT::value_type);
        FillSetsIntervals(*bv, vector_max/ 10, vector_max, fill_factor);
        *bv |= bv_common;

        target->push_back(std::move(*bv));
    } // for
}
