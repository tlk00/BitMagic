/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

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
unsigned orders_cnt      = customers_cnt * 10;
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
    TID64Map   lineitem_shipdate_bvmap; // lineitem shipdate index
    TIDMap     lineitem_order_bvmap;    // lineitem shipdate index
};


static
void OptimizeIDMap(TIDMap& id_map)
{
    for (TIDMap::iterator it = id_map.begin();
         it != id_map.end();
         ++it)
    {
        it->second.optimize();
    } // for
}

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
            
            std::cout << "Nation:" << nation_id << " " << supp50_10_bv.count() << std::endl;
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
            
            std::cout << "Nation:" << nation_id << " " << cust50_10_bv.count() << std::endl;
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
    
    std::cout << "Orders customers index size = " << ord.orders_customer_smap.size()
              << std::endl;
}

// Generate bit vector indexes for ORDER table
//
void GenerateLineItemIdx(LineItem& litem, const Order& ord, const Customer& cust)
{
    bm::random_subset<TBVector> rsub; // random subset sampler
    unsigned i, j;
    
    // fill in total vector of customers as a monotonically increasing id
    //
    TBVector& bv_litem = *litem.lineitem_total_bv;
    //bv_litem.set_range(0, lineitem_cnt-1, true);
    
    unsigned li_id = 0;
    unsigned year_from = 1994;
    unsigned date = 1;
    
    // for each order, generate line items
    for (i = 0; i < orders_cnt; ++i)
    {
        // lets assume that each order has 20 to 1 line items
        //
        unsigned li_cnt = rand() % 10;
        if (li_cnt == 0)
            li_cnt = 1;
        
        // pick a random year from 10
        unsigned order_year = year_from + (rand()%10);
        
        if (li_cnt < 3) // small order
        {
            // day of a year, assume all items shipped the same date
            unsigned li_day = rand() % 365;
            for (j = 0; j < li_cnt; ++j)
            {
                // pick a supplier
                unsigned supp_id = rand() % suppliers_cnt;
                
            }
        }
        
    } // for i

    
    for (i = 0; i < lineitem_cnt; ++i)
    {
    
    }
    
}

int main(int argc, char *argv[])
{
    Suppliers supp;
    Customer  cust;
    Order     ord;
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



