#ifndef BMSPARSEVEC_COMPR_H__INCLUDED__
#define BMSPARSEVEC_COMPR_H__INCLUDED__
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

#include <memory.h>

#ifndef BM_NO_STL
#include <stdexcept>
#endif

#include "bmsparsevec.h"

namespace bm
{


/*!
   \brief compressed sparse vector for NULL-able sparse vectors
 
   \ingroup svector
*/
template<class Val, class SV>
class compressed_sparse_vector
{
public:
    typedef Val                                      value_type;
    typedef bm::id_t                                 size_type;
    typedef SV                                       sparse_vector_type;
    typedef typename SV::bvector_type                bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
    typedef const value_type&                        const_reference;
    typedef typename bvector_type::allocator_type    allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
public:
    /*! Statistical information about  memory allocation details. */
    struct statistics : public bv_statistics
    {};
public:
    compressed_sparse_vector(allocation_policy_type ap = allocation_policy_type(),
                             size_type bv_max_size = bm::id_max,
                             const allocator_type&   alloc  = allocator_type());
    /*! copy-ctor */
    compressed_sparse_vector(const compressed_sparse_vector<Val, SV>& csv);
    
    /*! copy assignmment operator */
    compressed_sparse_vector<Val,SV>& operator = (const compressed_sparse_vector<Val, SV>& csv)
    {
        if (this != &csv)
        {
            sv_ = csv.sv_;
        }
        return *this;
    }
    
    /*!
        \brief set specified element with bounds checking and automatic resize
     
        Method cannot insert elements, so every new idx has to be greater or equal
        than what it used before.
     
        \param idx - element index
        \param v   - element value
    */
    void push_back(size_type idx, value_type v);


protected:
private:
    sparse_vector_type                    sv_;
    bm::id_t                              max_id_;
};

//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<class Val, class SV>
compressed_sparse_vector<Val, SV>::compressed_sparse_vector(allocation_policy_type ap,
                                                            size_type bv_max_size,
                                                            const allocator_type&   alloc)
    : sv_(bm::use_null, ap, bv_max_size, alloc),
      max_id_(0)
{
}

//---------------------------------------------------------------------

template<class Val, class SV>
compressed_sparse_vector<Val, SV>::compressed_sparse_vector(
                          const compressed_sparse_vector<Val, SV>& csv)
: sv_(csv.sv_),
  max_id_(csv.max_id_)
{
}

//---------------------------------------------------------------------

template<class Val, class SV>
void compressed_sparse_vector<Val, SV>::push_back(size_type idx, value_type v)
{
    if (sv_.empty())
    {
    }
    else
    if (idx <= max_id_)
    {
        sv_.throw_range_error("compressed sparse vector push_back() range error");
    }
    
    bvector_type* bv_null = sv_.get_null_bvect();
    BM_ASSERT(bv_null);
    bv_null->set(idx);
    sv_.push_back(v);
    
    max_id_ = idx;
}

//---------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"


#endif
