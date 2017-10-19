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

#include <bitset>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <sstream>

//#define BMSSE2OPT
#define BMSSE42OPT


#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4996)
#endif

#include <vector>

#include "bm.h"
#include "bmalgo.h"
#include "bmserial.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"

//#include "bmdbg.h"

#include <math.h>

using namespace std;

const unsigned BSIZE = 150000000;
const unsigned int REPEATS = 300;

typedef  bitset<BSIZE>  test_bitset;

unsigned platform_test = 1;


class TimeTaker
{
public:

    TimeTaker(const char* test_name, unsigned repeats) 
        : test_name_(test_name), repeats_(repeats) 
    {
        start_ = clock();
    }

    ~TimeTaker()
    {
        finish_ = clock();
        clock_t elapsed_clocks = finish_ - start_;
        double duration = (double)(finish_ - start_) / CLOCKS_PER_SEC;

        cout << test_name_ << " ; ";
        if (platform_test) 
        {
            cout << duration << endl;
            return;
        }

        cout <<  elapsed_clocks << ";" << duration << ";";
        if (repeats_)
        {
            double ops_per_sec = (double)repeats_ / duration;
            cout << ops_per_sec;
        }
        cout << endl;
    }

private:
    const char*  test_name_;
    clock_t      start_;
    clock_t      finish_;
    unsigned     repeats_;
};

typedef bm::bvector<> bvect;



void SimpleFillSets(test_bitset& bset, 
                       bvect& bv,
                       unsigned min, 
                       unsigned max,
                       unsigned fill_factor,
                       bool set_flag=true)
{
    for (unsigned i = min; i < max; i+=fill_factor)
    {
        bset[i] = set_flag;
        bv[i] = set_flag;
    } // for i
}


//
// Interval filling.
// 111........111111........111111..........11111111.......1111111...
//

void FillSetsIntervals(test_bitset& bset, 
                       bvect& bv,
                       unsigned min, 
                       unsigned max,
                       unsigned fill_factor,
                       bool set_flag=true)
{
    while(fill_factor==0)
    {
        fill_factor=rand()%10;
    }

    unsigned i, j;
    unsigned factor = 10 * fill_factor;
    for (i = min; i < max; ++i)
    {
        unsigned len, end; 

        do
        {
            len = rand() % factor;
            end = i+len;
            
        } while (end >= max);
        for (j = i; j < end; ++j)
        {
            if (set_flag)
            {
                bset[j] = true;
                bv[j]= true;
            }
            else
            {
                bset[j] = false;
                bv[j] = false;
            }
                           
        } // j

        i = end;


        len = rand() % 10;

        i+=len;

        {
            for(unsigned k=0; k < 1000 && i < max; k+=3,i+=3)
            {
                if (set_flag)
                {
                    bset[i] = true;
                    bv[i] = true;
                }
                else
                {
                    bset[j] = false;
                    bv[j] = false;
                }

            }
        }

    } // for i

}


void MemCpyTest()
{
    unsigned* m1 = new unsigned[BSIZE/32];
    unsigned* m2 = new unsigned[BSIZE/32];
    
    unsigned int i,j;

    if (!platform_test)
    {
    TimeTaker tt("Memory ADD transfer test", REPEATS * 4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        for (j = 0; j < BSIZE/32; j+=4)
        {
            m1[j+0] += m2[j+0];
            m1[j+1] += m2[j+1];
            m1[j+2] += m2[j+2];
            m1[j+3] += m2[j+3];
        }
    }
    }
    
    if (!platform_test)
    {
    TimeTaker tt("memcpy transfer test", REPEATS * 4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        memcpy(m1, m2, BSIZE/32 * sizeof(unsigned));
    }
    }
    
    delete [] m1;
    delete [] m2;
}


void BitCountTest()
{
    {
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned value = 0;

    FillSetsIntervals(*bset, *bv, 0, BSIZE, 10);

    if (!platform_test)
    {
    TimeTaker tt("BitCount. Random bitvector", REPEATS*2);
    for (unsigned i = 0; i < REPEATS*2; ++i)
    {    
        value+=bv->count();
    }
    }

    volatile unsigned* p = &value;
    unsigned c1;
    c1 = value = 0;

    if (!platform_test)
    {
    TimeTaker tt("BitCount. Random bitvector (STL)", REPEATS*2);
    for (unsigned i = 0; i < REPEATS*2; ++i)
    {    
        value += (unsigned)bset->count();
    }
    }

    c1 = *p;
    c1 = value = 0;
	stringstream s;
	s << value << c1; // to fool the optimization

    delete bset;
    delete bv;

    }
}


void BitForEachTest()
{
    // setup the test data
    //
    unsigned* test_arr = new unsigned[65536];
	for (unsigned j = 0; j < 65536; ++j)
	{
        test_arr[j] = j;		
	}

    if (!platform_test)
    {
	unsigned bit_list[32];
    TimeTaker tt("BitList algorithm. Conventional (AND based check)", REPEATS*10);
    

	for (unsigned i = 0; i < REPEATS*10; ++i)
    {    
		for (unsigned j = 0; j < 65536; ++j)
		{
			bm::bit_list(i*test_arr[j], bit_list);
		}
	}
	}
	
    {
	unsigned bit_list[32];
    TimeTaker tt("BitList4 algorithm(sub-octet+switch)", REPEATS*20);

	for (unsigned i = 0; i < REPEATS*100; ++i)
    {    
		for (unsigned j = 0; j < 65536; ++j)
		{
			bm::bit_list_4(i*test_arr[j], bit_list);
		}
	}
	}

	{
		unsigned bit_list[32];
		TimeTaker tt("BitScan on bitcount algorithm", REPEATS * 20);

		for (unsigned i = 0; i < REPEATS * 100; ++i)
		{
			for (unsigned j = 0; j < 65536; ++j)
			{
				bm::bitscan_popcnt(i*test_arr[j], bit_list);
			}
		}
	}

	delete [] test_arr;
}


void BitCountSparseTest()
{
    {
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned value = 0, c1;
    volatile unsigned* p = &value;

    SimpleFillSets(*bset, *bv, 0, BSIZE, 2500);

    {
    TimeTaker tt("BitCount: Sparse bitset ", REPEATS*10);
    for (unsigned i = 0; i < REPEATS*10; ++i)
    {    
        value += bv->count();
    }
    }

    if (!platform_test)
    {
    TimeTaker tt("BitCount: Sparse bitset (STL)", REPEATS*10);
    for (unsigned int i = 0; i < REPEATS*10; ++i)
    {    
        value += bset->count();
    }
    }

    c1 = *p;
    value = c1 = 0;
    
    BM_DECLARE_TEMP_BLOCK(tb)
    bv->optimize(tb);

    {
    TimeTaker tt("BitCount: GAP Sparse bitset", REPEATS*1000);
    for (unsigned i = 0; i < REPEATS*1000; ++i)
    {    
        value += bv->count();
    }
    delete bv;
    delete bset;
    }
    c1 = *p;
    value = c1 = 0;

    }

}



void BitCompareTest()
{
    {
    bvect*  bv1 = new bvect();
    bvect*  bv2 = new bvect();
    test_bitset*  bset = new test_bitset();
    int value = 0;

    SimpleFillSets(*bset, *bv1, 0, BSIZE, 10);
    SimpleFillSets(*bset, *bv2, 0, BSIZE, 10);

    {
    TimeTaker tt("BitCompare: Random bitvector", REPEATS*10);
    for (unsigned int i = 0; i < REPEATS*10; ++i)
    {    
        value+=bv1->compare(*bv2);
    }
    }

    delete bset;
    delete bv1;
    delete bv2;

    }

    if (platform_test) return;

    unsigned cnt = REPEATS * 100000;
    unsigned* arr1 = new unsigned[cnt];
    unsigned* arr2 = new unsigned[cnt];

    unsigned i;
    for (i = 0; i < cnt; ++i)
    {
        if ((rand() % 10) == 0)
        {
            arr1[i] = 0;
        }
        else 
        {
            arr1[i] = rand();
            arr2[i] = rand();   
        }
    }

    {
    TimeTaker tt("wordcmp complex: Random words comparison", cnt);

    for (i = 0; i < cnt; ++i)
    {    
        int res2 = bm::wordcmp(arr1[i], arr2[i]);
        int res = bm::wordcmp0(arr1[i], arr2[i]);

        if (res != res2)
        {
            cerr << "Incorrect result ! " << arr1[i] 
                 << "<=>" << arr2[i] << " res=" << res <<
                 endl;
            exit(1);
        }
    }
    }

    int c = 0;
    volatile void* p = &c;

    {
    TimeTaker tt("wordcmp0. Random words comparison", cnt);
    for (i = 0; i < cnt; ++i)
    {    
        c += bm::wordcmp0(arr1[i], arr2[i]);
    }
    }


    {
    TimeTaker tt("wordcmp. Random words comparison", cnt);
    for (i = 0; i < cnt; ++i)
    {    
        c += bm::wordcmp(arr1[i], arr2[i]);
    }
    }

    c = 0;;

    delete [] arr1;
    delete [] arr2;

    char buf[256];
    sprintf(buf, "%p", p);
}

void EnumeratorTest()
{
    bvect  bv;
    test_bitset*  bset = new test_bitset();
    unsigned value = 0;

    FillSetsIntervals(*bset, bv, 0, BSIZE, 10);

    unsigned cnt1 = bv.count();
    //unsigned cnt2 = bset->count();

    
    unsigned i;

    {
    TimeTaker tt("bvector<>::enumerator", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {    
        bvect::enumerator en = bv.first();

        for (;en.valid();++en)
        {
            value = *en;
        }
    }
    }


    // -----------------------------------------------

    unsigned cnt = 0;
    {
        TimeTaker tt("bvector<>::get_next()", REPEATS);
        for (i = 0; i < REPEATS; ++i)
        {
            if (bv.any())
            {
                unsigned v = bv.get_first();
                do
                {
                    v = bv.get_next(value);
                    cnt += v;
                } while(v);
            }
        }
    }

    delete bset;

    char buf[256];
    sprintf(buf, "%i %i ", cnt, cnt1);//, cnt2);

}


void EnumeratorTestGAP()
{
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned i;

    SimpleFillSets(*bset, *bv, 0, BSIZE, 2500);
    bv->count();

    for (int k = 0; k < 2; ++k)
    {

    {
    unsigned v = 0;
    TimeTaker tt("Sparse bvector (enumerator)", REPEATS*10*(k+1));
    for (i = 0; i < REPEATS*10*(k+1); ++i)
    {    
        bvect::enumerator en = bv->first();
        bvect::enumerator bend = bv->end();

        while (en < bend)
        {
            v += *en;
            ++en;
        }
    }

	stringstream s;
	s << v << endl; // attempt to fool optimization

    }

    // -----------------------------------------------

    unsigned cnt = 0;
    {
    TimeTaker tt("Sparse bvector (get_next())", REPEATS*10*(k+1));

    for (i = 0; i < REPEATS*10*(k+1); ++i)
    {
        if (bv->any())
        {
            unsigned v = bv->get_first();
            do
            {
                v = bv->get_next(v);
                cnt += v;
            } while (v);
        }
    }
    }
    char buf[256];
    sprintf(buf, "%i", cnt); // to fool some smart compilers like ICC

    {

    BM_DECLARE_TEMP_BLOCK(tb)
    bv->optimize(tb);
    }

    if (!platform_test) 
    {
        cout << "Testing optimized vectors." << endl;
    }
    }

    delete bv;
    delete bset;
    // -----------------------------------------------

}

void SerializationTest()
{
	bvect bv_sparse;
    // stack declaration of temp block instead of re-allocation makes things faster
    BM_DECLARE_TEMP_BLOCK(tb)


	// prepare a test bitset with a small number of bits set somewhere
	// far from the beginning
	for (unsigned i = 0; i < 5; ++i)
	{
		bv_sparse[BSIZE/2 + i * 3] = true;		
	}
	bv_sparse[100] = true;
	bv_sparse[70000] = true;
	bv_sparse[200000] = true;

	bv_sparse.optimize(tb);

	unsigned cnt = bv_sparse.count();
    bvect::statistics st;
    bv_sparse.calc_stat(&st);
    unsigned char*  buf = new unsigned char[st.max_serialize_mem];

    unsigned len, id_size;
    len = id_size = 0;
    {
	TimeTaker tt("Small bvector serialization", REPEATS*70000);
	for (unsigned i = 0; i < REPEATS*70000; ++i)
	{
		len += bm::serialize(bv_sparse, buf, tb, bm::BM_NO_BYTE_ORDER|bm::BM_NO_GAP_LENGTH);
		id_size += cnt * (unsigned)sizeof(unsigned);
	}
	}
	
	delete [] buf; buf = 0;
		
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned value = 0;

    SimpleFillSets(*bset, *bv, 0, BSIZE, 4);
    
	cnt = bv->count();
    bv->calc_stat(&st);
    buf = new unsigned char[st.max_serialize_mem];
    
    {
	TimeTaker tt("Large bvector serialization", REPEATS*4);
	for (unsigned i = 0; i < REPEATS*4; ++i)
	{
		len += bm::serialize(*bv, buf, tb, bm::BM_NO_BYTE_ORDER|bm::BM_NO_GAP_LENGTH);
		id_size += cnt * (unsigned)sizeof(unsigned);
	}
	}
    
	char cbuf[256];
	sprintf(cbuf, "%i %i %i", id_size, len, value);
    
    
    delete bv;
    delete bset;	
    delete [] buf;
}

void InvertTest()
{
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned i;
    //unsigned value = 0;

    SimpleFillSets(*bset, *bv, 0, BSIZE, 2500);
    {
    TimeTaker tt("Invert bvector", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv->flip();    
    }
    }

    if (!platform_test)
    {
    TimeTaker tt("Invert bvector (STL)", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bset->flip();    
    }
    }

    delete bv;
    delete bset;
}


void AndTest()
{
    bvect*  bv1 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    bvect*  bv2 = new bvect();
    unsigned i;
    //unsigned value = 0;

    SimpleFillSets(*bset1, *bv1, 0, BSIZE, 100);
    SimpleFillSets(*bset1, *bv2, 0, BSIZE, 100);
    {
    TimeTaker tt("AND bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        *bv1 &= *bv2;
    }
    }

    if (!platform_test)
    {
    TimeTaker tt("AND bvector test(STL)", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        *bset1 &= *bset2;
    }
    }

    delete bv1;
    delete bv2;

    delete bset1;
    delete bset2;
}


void SubTest()
{
    bvect*  bv1 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    bvect*  bv2 = new bvect();
    unsigned i;

    SimpleFillSets(*bset1, *bv1, 0, BSIZE, 100);
    SimpleFillSets(*bset1, *bv2, 0, BSIZE, 100);
    delete bset1;

    {
    TimeTaker tt("SUB bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        *bv1 -= *bv2;
    }
    }
        
    delete bv1;
    delete bv2;
}


void XorCountTest()
{
    bvect*  bv1 = new bvect();
    bvect*  bv2 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    unsigned i;

    SimpleFillSets(*bset1, *bv1, 0, BSIZE, 400);
    SimpleFillSets(*bset2, *bv2, 0, BSIZE, 500);

    unsigned count1 = 0;
    unsigned count2 = 0;
    unsigned test_count = 0;

    if (!platform_test)
    {
    bvect bv_tmp;
    TimeTaker tt("XOR COUNT bvector test with TEMP vector", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_tmp.clear(false);
        bv_tmp |= *bv1;
        bv_tmp ^= *bv2;
        count1 += bv_tmp.count();
    }
    }

    if (!platform_test)
    {
    test_bitset*  bset_tmp = new test_bitset();
    TimeTaker tt("XOR COUNT bvector test with TEMP vector (STL)", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bset_tmp->reset();
        *bset_tmp |= *bset1;
        *bset_tmp ^= *bset2;
        test_count += (unsigned)bset_tmp->count();
    }
    }


    {
    TimeTaker tt("XOR COUNT bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        count2 += bm::count_xor(*bv1, *bv2);
    }
    }

    
    if (!platform_test)    
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            cout << count1 << " " << count2 << " " << test_count << endl;
            exit(1);
        }
    count1 = count2 = 0;
    
    // -----------------------------------------
    if (!platform_test)
    {
        cout << "One optimized vector" << endl;
    }
    BM_DECLARE_TEMP_BLOCK(tb)
    bv2->optimize(tb);

    if (!platform_test)
    {
    bvect bv_tmp;
    TimeTaker tt("XOR COUNT bvector test with TEMP vector", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_tmp.clear(false);
        bv_tmp |= *bv1;
        bv_tmp ^= *bv2;
        count1 += (unsigned)bv_tmp.count();
    }
    }

    {
    TimeTaker tt("XOR COUNT bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        count2 += (unsigned)bm::count_xor(*bv1, *bv2);
    }
    }

    if (!platform_test)
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            exit(1);
        }
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "Both vectors optimized" << endl;
    }
    bv1->optimize(tb);
    //bv1->stat();
    if (!platform_test)
    {
    bvect bv_tmp;
    TimeTaker tt("XOR COUNT bvector test with TEMP vector", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_tmp.clear(false);
        bv_tmp |= *bv1;
        bv_tmp ^= *bv2;
        count1 += (unsigned)bv_tmp.count();
    }
    }

    {
    TimeTaker tt("XOR COUNT bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        count2 += (unsigned)bm::count_xor(*bv1, *bv2);
    }
    }
    if (!platform_test)
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            exit(1);
        }
    count1 = count2 = 0;


    delete bv1;
    delete bv2;
    
    delete bset1;
    delete bset2;    
}


void TI_MetricTest()
{
    bvect*  bv1 = new bvect();
    bvect*  bv2 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    unsigned i;

    SimpleFillSets(*bset1, *bv1, 0, BSIZE, 500);
    SimpleFillSets(*bset2, *bv2, 0, BSIZE, 250);

    unsigned count1 = 0;
    unsigned count2 = 0;
    unsigned countA=0, countB=0, test_countA=0, test_countB=0;
    unsigned test_count = 0;
    double ti1=0, ti2=0;

    {
    TimeTaker tt("Tversky Index bvector test vector", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        count1 = bm::count_and(*bv1, *bv2);
        
        countA = bm::count_sub(*bv1, *bv2);
        countB = bm::count_sub(*bv2, *bv1);
        
        ti1 = double(count1) / double(0.4*countA + 0.5*countB + count1);
    }
    }


    if (!platform_test)
    {
    test_bitset*  bset_tmp = new test_bitset();
    double test_dice = 0;
    TimeTaker tt("Dice bvector test with TEMP vector(STL)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bset_tmp->reset();
        *bset_tmp |= *bset1;
        *bset_tmp &= *bset2;
        test_count += (unsigned)bset_tmp->count();
        
        test_countA += (unsigned)bset1->count();
        test_countB += (unsigned)bset2->count();
        
        test_countA += (unsigned)bset1->count();
        test_countB += (unsigned)bset2->count();
        
        test_dice += double(2*test_count) / double(test_countA + test_countB);
    }
    }


    {
    bm::distance_metric_descriptor dmd[3];
    dmd[0].metric = bm::COUNT_AND;
    dmd[1].metric = bm::COUNT_SUB_AB;
    dmd[2].metric = bm::COUNT_SUB_BA;    
    
    TimeTaker tt("Tversky Index bvector test (pipeline)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bm::distance_operation(*bv1, *bv2, &dmd[0], (&dmd[0])+3);
                
        ti2 = double(dmd[0].result) / double(0.4*dmd[1].result + 0.5*dmd[2].result + dmd[0].result);
        
        dmd[0].result = dmd[1].result = dmd[2].result = 0;
    }
    }

    
    if (fabs(ti2 - ti1) > 0.1)
    {
        cout << "Check failed ! error=" << fabs(ti2 - ti1) << endl;
        cout << ti1 << " " << ti2 << endl;
        exit(1);
    }
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "One optimized vector" << endl;
    }
    BM_DECLARE_TEMP_BLOCK(tb)
    bv2->optimize(tb);
    bv1->count(); // trying to fool the CPU cache

    
    {
    TimeTaker tt("Dice metric bvector test", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        count1 = bm::count_and(*bv1, *bv2);
        
        countA = bm::count_sub(*bv1, *bv2);
        countB = bm::count_sub(*bv2, *bv1);
        
        ti1 = double(count1) / double(0.4*countA + 0.5*countB + count1);
    }
    }



    {
    bm::distance_metric_descriptor dmd[3];
    dmd[0].metric = bm::COUNT_AND;
    dmd[1].metric = bm::COUNT_SUB_AB;
    dmd[2].metric = bm::COUNT_SUB_BA;    
    
    TimeTaker tt("Tversky Index bvector test(pipeline)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bm::distance_operation(*bv1, *bv2, &dmd[0], (&dmd[0])+3);
                
        ti2 = double(dmd[0].result) / double(0.4*dmd[1].result + 0.5*dmd[2].result + dmd[0].result);
        
        dmd[0].result = dmd[1].result = dmd[2].result = 0;
    }
    }


    if (fabs(ti2 - ti1) > 0.1)
    {
        cout << "Check failed !" << endl;
        cout << ti1 << " " << ti2 << endl;
        exit(1);
    }
    count1 = count2 = 0;
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "Both vectors optimized" << endl;
    }
    bv1->optimize(tb);

    {
    TimeTaker tt("Tversky index bvector test", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        count1 = bm::count_and(*bv1, *bv2);
        
        countA = bm::count_sub(*bv1, *bv2);
        countB = bm::count_sub(*bv2, *bv1);
        
        ti1 = double(count1) / double(0.4*countA + 0.5*countB + count1);
    }
    }

    {
    bm::distance_metric_descriptor dmd[3];
    dmd[0].metric = bm::COUNT_AND;
    dmd[1].metric = bm::COUNT_SUB_AB;
    dmd[2].metric = bm::COUNT_SUB_BA;    
    
    TimeTaker tt("Tversky Index bvector test (pipeline)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bm::distance_operation(*bv1, *bv2, &dmd[0], (&dmd[0])+3);
                
        ti2 = double(dmd[0].result) / double(0.4*dmd[1].result + 0.5*dmd[2].result + dmd[0].result);
        
        dmd[0].result = dmd[1].result = dmd[2].result = 0;
    }
    }

    if (fabs(ti2 - ti1) > 0.1)
    {
        cout << "Check failed !" << endl;
        cout << ti1 << " " << ti2 << endl;
        exit(1);
    }


    delete bv1;
    delete bv2;
    
    delete bset1;
    delete bset2;    
}

void BitBlockTransposeTest()
{
	bm::word_t BM_ALIGN16 block1[bm::set_block_size] BM_ALIGN16ATTR = { 0, };
    //bm::word_t BM_ALIGN16 block2[bm::set_block_size] = {0xFF,};
    unsigned   BM_ALIGN16 tmatrix1[32][bm::set_block_plain_size] BM_ALIGN16ATTR;

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        block1[i] = 1 | (1 << 5) | (7 << 15) | (3 << 22);
    }

    const unsigned blocks_count = 70000;
    bm::word_t* blocks[blocks_count];
    for (unsigned k = 0; k < blocks_count; ++k)
    {
        blocks[k] = bm::block_allocator::allocate(bm::set_block_size, 0);
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            blocks[k][i] = 1 | (1 << 5) | (7 << 15) | (3 << 22);
        }
    }

    unsigned cnt=0;
/*
    {
    TimeTaker tt("Bit-block transpose.", REPEATS*1000);
    for (unsigned i = 0; i < REPEATS*1000; ++i)
    {
        bm::bit_block_transpose(block1, tmatrix1);
    }
    }

    {
    TimeTaker tt("Bit-block trestore.", REPEATS*1000);
    for (unsigned i = 0; i < REPEATS*1000; ++i)
    {
        bm::bit_block_trestore(tmatrix1, block2);
        cnt += block2[10];

    }
    }

    {
    TimeTaker tt("Bit-block transpose distance.", REPEATS*1000);
    unsigned distance[bm::set_block_plain_cnt][bm::set_block_plain_cnt];
    for (unsigned i = 0; i < REPEATS*1000; ++i)
    {
        bm::bit_block_tmatrix_distance(tmatrix1, distance);
        cnt += distance[1][1];
    }
    }
    printf("", cnt);
*/

    unsigned d2[bm::set_block_plain_cnt][bm::set_block_plain_cnt];
    {
    TimeTaker tt("Bit-block transpose+distance", 100000);
    unsigned distance[bm::set_block_plain_cnt][bm::set_block_plain_cnt];
    unsigned idx = 0;
    for (unsigned i = 0; i < 100000; ++i)
    {
        bm::vect_bit_transpose<unsigned, 
                               bm::set_block_plain_cnt, 
                               bm::set_block_plain_size>
                               (blocks[idx], bm::set_block_size, tmatrix1);
        bm::tmatrix_distance<unsigned, 
                             bm::set_block_plain_cnt, 
                             bm::set_block_plain_size>
                             (tmatrix1, distance);
    
        cnt += distance[1][1];
        ++idx;
        if (idx >= blocks_count) idx = 0;
        memcpy(d2, distance, sizeof(distance));
    }

    }
    
    char cbuf[256];
    sprintf(cbuf, "%i %i", cnt, d2[10][10]);

    for (unsigned i = 0; i < blocks_count; ++i)
    {
        bm::block_allocator::deallocate(blocks[i], 0);
    }

}

void ptest()
{
    bvect*  bv_small = new bvect(bm::BM_GAP);
    bvect*  bv_large = new bvect(bm::BM_GAP);

    test_bitset*  bset = new test_bitset();

    FillSetsIntervals(*bset, *bv_large, 0, 2000000000, 10);

    for (int i = 0; i < 2000; ++i)
    {
        bv_small->set(i*10000);
    }

    {
    TimeTaker tt("Operation &= test", REPEATS * 10);
    unsigned count = 0;
    for (unsigned i = 0; i < REPEATS*10; ++i)
    {
        bvect t1(bm::BM_GAP);
        t1 = *bv_small;
        t1 &= *bv_large;
        count += t1.count();
    }
    }



    {
    TimeTaker tt("Operation &= with enumerator test", REPEATS * 10);
    unsigned count = 0;
    for (unsigned i = 0; i < REPEATS*10; ++i)
    {
        bvect t1(bm::BM_GAP);
        bvect t2(bm::BM_GAP);
        t1 = *bv_small;
        
        for (bvect::enumerator it = t1.first(); it != t1.end(); ++it) {
            if ((*bv_large)[*it]) {
                t2.set_bit(*it);
            }
        }
        count += t2.count();
    }
    }


}


typedef bm::sparse_vector<unsigned, bvect> svect;

// create a benchmark vector with a few dufferent distribution patterns
//

void FillSparseIntervals(svect& sv)
{
    sv.resize(250000000);
    
    unsigned i;
    for (i = 256000; i < 256000 * 2; ++i)
    {
        sv.set(i, 0xFFE);
    }
    for (i = 256000 * 3; i < 256000 * 5; ++i)
    {
        sv.set(i, i);
    }
    
    for (i = 180000000; i < 190000000; ++i)
    {
        sv.set(i, rand() % 128000);
    }
    
    for (i = 200000000; i < 210000000; ++i)
    {
        sv.set(i, rand() % 128000);
    }
}


void SparseVectorAccessTest()
{
    std::vector<unsigned> target;
    svect   sv1;

    FillSparseIntervals(sv1);
    BM_DECLARE_TEMP_BLOCK(tb)
    sv1.optimize(tb);
    target.resize(250000000);

    unsigned long long cnt = 0;
    {
    TimeTaker tt("sparse_vector random element access test", REPEATS/10 );
    for (unsigned i = 0; i < REPEATS/10; ++i)
    {
        for (unsigned j = 256000; j < 190000000/2; ++j)
        {
            unsigned v = sv1[j];
            cnt += v;
        }
    }
    }
    {
    TimeTaker tt("sparse_vector extraction access test", REPEATS/10 );
    for (unsigned i = 0; i < REPEATS/10; ++i)
    {
        unsigned target_size = 190000000/2 - 256000;
        sv1.extract(&target[0], 256000, target_size);
    }
    }
    
    char buf[256];
    sprintf(buf, "%i", (int)cnt); // to fool some smart compilers like ICC
    
}


int main(void)
{
//    ptest();

    TimeTaker tt("TOTAL", 1);

    MemCpyTest();

    BitCountTest();

    BitForEachTest();

    BitCountSparseTest();

    BitCompareTest();

    BitBlockTransposeTest();

    EnumeratorTest();

    EnumeratorTestGAP();

    AndTest();
    SubTest();  

    InvertTest();  

    XorCountTest();

    TI_MetricTest();

    SerializationTest();
    
    SparseVectorAccessTest();
    
    return 0;
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif



