/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/*
    Benchmark application, simulates snowflake query accelaration from TPC-H.
*/

#include <iostream>
#include <chrono>
#include <time.h>
#include <stdio.h>


#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4996)
#endif

#include <vector>
#include <chrono>
#include <map>

#include "bm.h"
#include "bmalgo.h"
#include "bmserial.h"
#include "bmrandom.h"
#include "bmsparsevec.h"


//#include "bmdbg.h"
#include "bmtimer.h"


void show_help()
{
    std::cerr
      << "BitMagic benchmark (analytical) (c) 2017."  << std::endl
      << std::endl
      ;
}

// Benchmark parameters
//

bool is_timing = false;

unsigned nations_cnt = 200;
unsigned nations_top_cnt = 10;
unsigned suppliers_cnt   = 100000;
unsigned customers_cnt   = 1500000;
unsigned orders_cnt      = customers_cnt * 5;
//unsigned lineitem_cnt    = orders_cnt * 6;



int parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help"))
        {
            show_help();
            return 0;
        }
    } // for i
    return 0;
}


// Globals
//

typedef  bm::bvector<>   TBVector;
typedef  std::map<unsigned, TBVector> TIDMap;
typedef  std::map<uint64_t, TBVector> TID64Map;
typedef  std::map<unsigned, std::vector<char> > TIDSMap;
typedef  std::map<uint64_t, std::vector<char> > TID64SMap;


bm::chrono_taker::duration_map_type  timing_map;

// collection of indexes for Supplier table
//
struct Suppliers
{
    Suppliers()
    {
        suppliers_total_bv = new TBVector(bm::BM_GAP);
    }
    
    ~Suppliers()
    {
        delete suppliers_total_bv;
    }

    TBVector*  suppliers_total_bv;
    TIDMap     suppliers_nations_bvmap; ///< NATION_KEY to SUPPLIER
};

// collection of indexes for Customer table
//
struct Customer
{
    Customer()
    {
        customers_total_bv = new TBVector(bm::BM_GAP);
    }
    
    ~Customer()
    {
        delete customers_total_bv;
    }

    TBVector*  customers_total_bv;
    TIDMap     customers_nations_bvmap; ///< NATION_KEY to customer
};

// collection of indexes for Order table
//
struct Order
{
    Order()
    {
        orders_total_bv = new TBVector(bm::BM_GAP);
    }
    
    ~Order()
    {
        delete orders_total_bv;
    }

    TBVector*  orders_total_bv;
    TIDSMap    orders_customer_smap; // ORDER to customer
};

// collection of indexes for LineItem table
//
struct LineItem
{
    LineItem()
    {
        lineitem_total_bv = new TBVector(bm::BM_GAP);
    }
    
    ~LineItem()
    {
        delete lineitem_total_bv;
    }

    TBVector*  lineitem_total_bv;
    TID64Map   lineitem_shipdate_bvmap; // shipdate index
    TID64SMap  lineitem_shipdate_smap; // shipdate index (compressed)

    TIDMap     lineitem_supplier_bvmap; // supplier index
    TIDSMap    lineitem_supplier_smap;  // supplier index (compressed)
    
    TIDMap     lineitem_order_bvmap;    // order index
    TIDSMap    lineitem_order_smap;     // order index (compressed)
};


template<typename TM>
void OptimizeIDMap(TM& id_map, bool opt_gap = false)
{
    for (typename TM::iterator it = id_map.begin();
         it != id_map.end();
         ++it)
    {
        it->second.optimize();
        if (opt_gap)
        {
            it->second.optimize_gap_size();
        }
    } // for
}


static
size_t ComputeIndexSize(const TIDSMap& sm)
{
    size_t s = 0;
    for (TIDSMap::const_iterator it = sm.begin();
         it != sm.end();
         ++it)
    {
        const std::vector<char>& v = it->second;
        if (v.size()==0)
        {
            std::cerr << "Empty vector found!" << std::endl;
            exit(1);
        }
        s += v.size();
    } // for
    return s;
}

/// bit-vector serializer works over temp buffer but the final result gets
/// saved in a smaller target buffer
/// std::vector<char> is used as a simple dynamic buffer class
///
static
void SerializeBVector(bm::serializer<TBVector>& bvs,
                                TBVector& bv,
                                std::vector<char>& temp_buf_vect,
                                std::vector<char>& buf_vect
                                )
{
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::bvector<>::statistics st;
    bv.optimize(tb, bm::bvector<>::opt_compress, &st);
    
    temp_buf_vect.resize(st.max_serialize_mem);

    unsigned len = bvs.serialize(bv,
                                 (unsigned char*)&temp_buf_vect[0],
                                 st.max_serialize_mem);

    buf_vect.resize(len);
    
    ::memcpy(&buf_vect[0], &temp_buf_vect[0], len);
    
    // integrity check
    #ifdef DEBUG
    TBVector bv1;
    bm::deserialize(bv1, (unsigned char*)&buf_vect[0]);
    if (bv.compare(bv1)!=0)
    {
        std::cerr << "Deserialization check failed!" << std::endl;
        exit(1);
    }
    #endif
}

/// Function assumes that src_vect contains a serialized BVector
/// it applies merge/join/OR between serialized and argument bv
/// then performs serialization
///
static
void SerializeORBVector(bm::serializer<TBVector>& bvs,
                                TBVector& bv,
                                std::vector<char>& temp_buf_vect,
                                std::vector<char>& buf_vect,
                                const std::vector<char>& src_vect
                                )
{
    // bv = bv OR src_vect
    if (src_vect.size() != 0)
        bm::deserialize(bv, (unsigned char*)&src_vect[0]);
    
    SerializeBVector(bvs, bv, temp_buf_vect, buf_vect);
}

template<typename TM1, typename TM2>
void SerializeMergeIDMap(bm::serializer<TBVector>& bvs,
                         std::vector<char>&        temp_buf_vect,
                         TM1&                      id_map,
                         TM2&                      id_smap)
{
    for (typename TM1::iterator it = id_map.begin();
         it != id_map.end();
         ++it)
    {
        TBVector &bv = it->second;
        typename TM1::key_type id = it->first;
        
        std::vector<char>& buf_vect = id_smap[id];
        SerializeORBVector(bvs, bv, temp_buf_vect, buf_vect, buf_vect);
    } // for
    
    id_map.clear();
}




// Generate bit vector indexes for SUPPLIER table
//
void GenerateSuppliersIdx(Suppliers& sup)
{
    bm::random_subset<TBVector> rsub; // random subset sampler
    unsigned i;
    
    
    // fill in total vector of suppliers
    // as a monotonically increasing id
    //
    TBVector& bv_supp = *sup.suppliers_total_bv;
    bv_supp.set_range(0, suppliers_cnt-1, true);

    
    // randomly sample 50% of suppliers
    // assume they belong to 10 top nations
    //
    TBVector supp50p_bv;
    rsub.sample(supp50p_bv, bv_supp, suppliers_cnt / 2);


    unsigned big10_cnt = supp50p_bv.count() / 10; // target count for top 10 members
    
    TBVector assigned_supp_bv;
    
    // generate nations to supplier inverted vector
    //
    for (i = 0; i < nations_top_cnt; ++i)
    {
        unsigned nation_id = i;
        {
            TBVector supp50_10_bv;
            rsub.sample(supp50_10_bv, supp50p_bv, big10_cnt);
            
            //std::cout << "Nation:" << nation_id << " " << supp50_10_bv.count() << std::endl;
            assigned_supp_bv |= supp50_10_bv;
            supp50p_bv -= supp50_10_bv;
            //std::cout << "Remained:" << supp50p_bv.count() << std::endl;
            
            sup.suppliers_nations_bvmap[nation_id] = supp50_10_bv;
        }
        if (!supp50p_bv.any())
        {
            break;
        }
    } // for i
    
    std::cout << "Assigned suppliers to Big10:" << assigned_supp_bv.count()
              << std::endl;
    
    // un-assigned = total MINUS assigned
    TBVector un_assigned_supp_bv = bv_supp - assigned_supp_bv;
    
    unsigned minor_nations_cnt = nations_cnt - nations_top_cnt;
    unsigned minor_nation_sample = un_assigned_supp_bv.count() / minor_nations_cnt;
    
    std::cout << "Suppliers per minor nation: " << minor_nation_sample << std::endl;
    
    for (i = nations_top_cnt; i < nations_cnt; ++i)
    {
        unsigned nation_id = i;
        {
            TBVector minor_supp_bv;
            rsub.sample(minor_supp_bv, un_assigned_supp_bv, minor_nation_sample);
            
            un_assigned_supp_bv -= minor_supp_bv;
            sup.suppliers_nations_bvmap[nation_id] = minor_supp_bv;
        }
        if (!un_assigned_supp_bv.any())
        {
            break;
        }
    } // for i
    
    if (un_assigned_supp_bv.any()) // in case of a rounding error
    {
        sup.suppliers_nations_bvmap[rand()%nations_cnt] |= un_assigned_supp_bv;
    }

    std::cout << "Nations suppliers index size = " << sup.suppliers_nations_bvmap.size()
              << std::endl;
    OptimizeIDMap(sup.suppliers_nations_bvmap);
}

// Generate bit vector indexes for CUSTOMER table
//
void GenerateCustomersIdx(Customer& cust)
{
    bm::random_subset<TBVector> rsub; // random subset sampler
    unsigned i;
    
    // fill in total vector of customers as a monotonically increasing id
    //
    TBVector& bv_cust = *cust.customers_total_bv;
    bv_cust.set_range(0, customers_cnt-1, true);
    
    // randomly sample 50% of suppliers
    // assume they belong to 10 top nations
    //
    TBVector cust50p_bv;
    rsub.sample(cust50p_bv, bv_cust, customers_cnt / 2);


    unsigned big10_cnt = cust50p_bv.count() / 10; // target count for top 10 members
    
    TBVector assigned_cust_bv;
    
    // generate nations to supplier inverted vector
    //
    for (i = 0; i < nations_top_cnt; ++i)
    {
        unsigned nation_id = i;
        {
            TBVector cust50_10_bv;
            rsub.sample(cust50_10_bv, cust50p_bv, big10_cnt);
            
            //std::cout << "Nation:" << nation_id << " " << cust50_10_bv.count() << std::endl;
            assigned_cust_bv |= cust50_10_bv;
            cust50p_bv -= cust50_10_bv;
            //std::cout << "Remained:" << supp50p_bv.count() << std::endl;
            
            cust.customers_nations_bvmap[nation_id] = cust50_10_bv;
        }
        if (!cust50p_bv.any())
        {
            break;
        }
    } // for i
    
    std::cout << "Assigned customers to Big10:" << assigned_cust_bv.count()
              << std::endl;

    // un-assigned = total MINUS assigned
    TBVector un_assigned_cust_bv = bv_cust - assigned_cust_bv;
    
    unsigned minor_nations_cnt = nations_cnt - nations_top_cnt;
    unsigned minor_nation_sample = un_assigned_cust_bv.count() / minor_nations_cnt;
    
    std::cout << "Customers per minor nation: " << minor_nation_sample << std::endl;
    
    for (i = nations_top_cnt; i < nations_cnt; ++i)
    {
        unsigned nation_id = i;
        {
            TBVector minor_cust_bv;
            rsub.sample(minor_cust_bv, un_assigned_cust_bv, minor_nation_sample);
            
            un_assigned_cust_bv -= minor_cust_bv;
            cust.customers_nations_bvmap[nation_id] = minor_cust_bv;
        }
        if (!un_assigned_cust_bv.any())
        {
            break;
        }
    } // for i
    
    if (un_assigned_cust_bv.any()) // in case of a rounding error
    {
        cust.customers_nations_bvmap[rand()%nations_cnt] |= un_assigned_cust_bv;
    }

    std::cout << "Nations customers index size = " << cust.customers_nations_bvmap.size()
              << std::endl;
    
    OptimizeIDMap(cust.customers_nations_bvmap);
}

// Generate bit vector indexes for ORDER table
//
void GenerateOrdersIdx(Order& ord, const Customer& cust)
{
    bm::serializer<bm::bvector<> > bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);
    bvs.set_compression_level(4);
    
    std::vector<char> temp_buf_vect;
    std::vector<char> buf_vect;


    bm::random_subset<TBVector> rsub; // random subset sampler
    unsigned i;
    
    // fill in total vector of customers as a monotonically increasing id
    //
    TBVector& bv_ord = *ord.orders_total_bv;
    bv_ord.set_range(0, orders_cnt-1, true);
    
    // lets assume uniform order distribution
    unsigned orders_per_cust = orders_cnt / customers_cnt;
    TBVector ord_bv(bm::BM_GAP, bm::gap_len_table_min<true>::_len, orders_cnt + 65535);

    TBVector undistr_ord_bv(bv_ord);
    for (i = 0; i < customers_cnt; ++i)
    {
        unsigned cust_id = i;
        
        {
            ord_bv.clear();
            rsub.sample(ord_bv, undistr_ord_bv, orders_per_cust);
            undistr_ord_bv -= ord_bv;
            
            std::vector<char> buf_vect;
            SerializeBVector(bvs, ord_bv, temp_buf_vect, buf_vect);
            
            ord.orders_customer_smap[cust_id] = buf_vect;
            
#ifdef DEBUG
            // integrity check
            std::vector<char>& bufv = ord.orders_customer_smap[cust_id];
            if (bufv.size() != buf_vect.size())
            {
                std::cout << "Problem!" << std::endl;
                exit(1);
            }
            TBVector bv1;
            bm::deserialize(bv1, (unsigned char*)&bufv[0]);
            if (bv1.compare(ord_bv)!=0)
            {
                std::cerr << "Deserialization check failed!" << std::endl;
                exit(1);
            }
#endif

        }
        if ((i % 10000)==0)
        {
            std::cout << "\r" << i << "/" << customers_cnt << " " << std::flush;
        }
    } // for
    
    if (undistr_ord_bv.any()) // in case of a rounding error
    {
        //ord.orders_customer_smap[rand()%customers_cnt] |= undistr_ord_bv;
    }
    std::cout << std::endl;
    
    std::cout << "Orders count = " << ord.orders_total_bv->count() << std::endl;
    std::cout << "Orders customers index size = " << ord.orders_customer_smap.size()
              << std::endl;
    
    size_t buf_sum = ComputeIndexSize(ord.orders_customer_smap);
    std::cout << "Orders customers index mem.size = " << buf_sum
              << std::endl;

}

// Generate bit vector indexes for ORDER table
//
void GenerateLineItemIdx(LineItem& litem, const Order& ord, const Customer& cust)
{
    bm::random_subset<TBVector> rsub; // random subset sampler
    
    bm::serializer<bm::bvector<> > bvs; // bit vector serialization utility
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);
    bvs.set_compression_level(4);

    std::vector<char> temp_buf_vect;

    unsigned i, j;
    
    // fill in total vector of customers as a monotonically increasing id
    //
    TBVector& bv_litem = *litem.lineitem_total_bv;
    //bv_litem.set_range(0, lineitem_cnt-1, true);
    
    unsigned li_id = 0;
    unsigned year_from = 1994;
    
    // for each order, generate line items
    //
    for (i = 0; i < orders_cnt; ++i)
    {
        unsigned order_id = i;
        // lets assume that each order has 20 to 1 line items
        //
        unsigned li_cnt = rand() % 7;
        if (li_cnt == 0)
            li_cnt = 1;
        
        // pick a random year from 10
        uint64_t order_year = year_from + (rand()%10);
        
        TBVector& bv_order = litem.lineitem_order_bvmap[order_id];

        if (li_cnt < 4) // small order
        {
            // day of a year, assume all items shipped the same date
            unsigned li_day = rand() % 365;
            // pick a random supplier (one for all items)
            unsigned supp_id = rand() % suppliers_cnt;
            uint64_t li_date = (order_year << 32) | li_day;
            
            TBVector& bv_date = litem.lineitem_shipdate_bvmap[li_date];
            TBVector& bv_supp = litem.lineitem_supplier_bvmap[supp_id];

            for (j = 0; j < li_cnt; ++j)
            {
                bv_date[li_id] = true;
                bv_supp[li_id] = true;
                bv_litem[li_id] = true;
                bv_order[li_id] = true;
                ++li_id;
            }
        }
        else
        {
            // day of a year, assume all items shipped the same date
            unsigned li_day = rand() % 365;
            
            // pick a random supplier (one for all items)
            unsigned supp_id = rand() % suppliers_cnt;
            
            for (j = 0; j < li_cnt; ++j)
            {
                uint64_t li_date = (order_year << 32) | li_day;
                TBVector& bv_date = litem.lineitem_shipdate_bvmap[li_date];
                TBVector& bv_supp = litem.lineitem_supplier_bvmap[supp_id];
                
                bv_date[li_id] = true;
                bv_supp[li_id] = true;
                bv_litem[li_id] = true;
                bv_order[li_id] = true;

                if (rand()%3 == 0)
                {
                    ++li_day;
                    supp_id = rand() % suppliers_cnt;
                }
                if (li_day > 365) { li_day = 1; ++order_year; }
                ++li_id;
            } // for j
        }
        
        // periodic index compression
        if ((i % 200000) == 0)
        {
            std::cout << "\r" << i << "[OPT]" << std::flush;
            SerializeMergeIDMap(bvs,
                                temp_buf_vect,
                                litem.lineitem_shipdate_bvmap,
                                litem.lineitem_shipdate_smap);
            SerializeMergeIDMap(bvs,
                                temp_buf_vect,
                                litem.lineitem_supplier_bvmap,
                                litem.lineitem_supplier_smap);
            
            SerializeMergeIDMap(bvs,
                                temp_buf_vect,
                                litem.lineitem_order_bvmap,
                                litem.lineitem_order_smap);
            std::cout << "\r" << i << "      " << std::flush;
        }
        
    } // for i
    std::cout << std::endl;
    
    SerializeMergeIDMap(bvs,
                        temp_buf_vect,
                        litem.lineitem_shipdate_bvmap,
                        litem.lineitem_shipdate_smap);

    SerializeMergeIDMap(bvs,
                        temp_buf_vect,
                        litem.lineitem_supplier_bvmap,
                        litem.lineitem_supplier_smap);
    
    SerializeMergeIDMap(bvs,
                        temp_buf_vect,
                        litem.lineitem_order_bvmap,
                        litem.lineitem_order_smap);

    std::cout
      << "Lineitems count = " << litem.lineitem_total_bv->count() << std::endl
      << "Lineitems shipdate index size = " << litem.lineitem_shipdate_smap.size() << std::endl
      << "Lineitems supplier index size = " << litem.lineitem_supplier_smap.size() << std::endl
      << "Lineitems order index size = " << litem.lineitem_order_smap.size()
      << std::endl;

}

int main(int argc, char *argv[])
{
    Suppliers supp;
    Customer  cust;
    Order     ord;
    LineItem  lineitem;
/*
    if (argc < 3)
    {
        show_help();
        return 1;
    }
*/
    try
    {
    /*
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;
    */
        
        GenerateSuppliersIdx(supp);
        GenerateCustomersIdx(cust);
        GenerateOrdersIdx(ord, cust);
        GenerateLineItemIdx(lineitem, ord, cust);
        
        getchar();

        
        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Timings (ms):" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif



