#ifndef BMSERIAL__H__INCLUDED__
#define BMSERIAL__H__INCLUDED__
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

/*! \file bmserial.h
    \brief Serialization / compression of bvector<>. Set operations on compressed BLOBs.
*/

/*! 
    \defgroup bvserial bvector<> serialization
    Serialization for bvector<> container
    \ingroup bvector
*/

#ifndef BM__H__INCLUDED__
#define BM__H__INCLUDED__

#include "bm.h"

#endif

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4311 4312 4127)
#endif



#include "encoding.h"
#include "bmdef.h"
#include "bmfunc.h"
#include "bmtrans.h"
#include "bmalgo_impl.h"
#include "bmutil.h"
#include "bmbuffer.h"

//#include "bmgamma.h"


namespace bm
{


// Serialization stream markup constants


const unsigned char set_block_end               = 0;   //!< End of serialization
const unsigned char set_block_1zero             = 1;   //!< One all-zero block
const unsigned char set_block_1one              = 2;   //!< One block all-set (1111...)
const unsigned char set_block_8zero             = 3;   //!< Up to 256 zero blocks
const unsigned char set_block_8one              = 4;   //!< Up to 256 all-set blocks
const unsigned char set_block_16zero            = 5;   //!< Up to 65536 zero blocks
const unsigned char set_block_16one             = 6;   //!< UP to 65536 all-set blocks
const unsigned char set_block_32zero            = 7;   //!< Up to 4G zero blocks
const unsigned char set_block_32one             = 8;   //!< UP to 4G all-set blocks
const unsigned char set_block_azero             = 9;   //!< All other blocks zero
const unsigned char set_block_aone              = 10;  //!< All other blocks one
const unsigned char set_block_bit               = 11;  //!< Plain bit block
const unsigned char set_block_sgapbit           = 12;  //!< SGAP compressed bitblock
const unsigned char set_block_sgapgap           = 13;  //!< SGAP compressed GAP block
const unsigned char set_block_gap               = 14;  //!< Plain GAP block
const unsigned char set_block_gapbit            = 15;  //!< GAP compressed bitblock 
const unsigned char set_block_arrbit            = 16;  //!< List of bits ON
const unsigned char set_block_bit_interval      = 17; //!< Interval block
const unsigned char set_block_arrgap            = 18;  //!< List of bits ON (GAP block)
const unsigned char set_block_bit_1bit          = 19; //!< Bit block with 1 bit ON
const unsigned char set_block_gap_egamma        = 20; //!< Gamma compressed GAP block
const unsigned char set_block_arrgap_egamma     = 21; //!< Gamma compressed delta GAP array
const unsigned char set_block_bit_0runs         = 22; //!< Bit block with encoded zero intervals
const unsigned char set_block_arrgap_egamma_inv = 23; //!< Gamma compressed inverted delta GAP array
const unsigned char set_block_arrgap_inv        = 24;  //!< List of bits OFF (GAP block)


/// \internal
/// \ingroup bvserial 
enum serialization_header_mask {
    BM_HM_DEFAULT = 1,
    BM_HM_RESIZE  = (1 << 1), ///< resized vector
    BM_HM_ID_LIST = (1 << 2), ///< id list stored
    BM_HM_NO_BO   = (1 << 3), ///< no byte-order
    BM_HM_NO_GAPL = (1 << 4)  ///< no GAP levels
};



#define SER_NEXT_GRP(enc, nb, B_1ZERO, B_8ZERO, B_16ZERO, B_32ZERO) \
   if (nb == 1) \
      enc.put_8(B_1ZERO); \
   else if (nb < 256) \
   { \
      enc.put_8(B_8ZERO); \
      enc.put_8((unsigned char)nb); \
   } \
   else if (nb < 65536) \
   { \
      enc.put_8(B_16ZERO); \
      enc.put_16((unsigned short)nb); \
   } \
   else \
   {\
      enc.put_8(B_32ZERO); \
      enc.put_32(nb); \
   }


#define BM_SET_ONE_BLOCKS(x) \
    {\
         unsigned end_block = i + x; \
         for (;i < end_block; ++i) \
            bman.set_block_all_set(i); \
    } \
    --i




/**
    Bit-vector serialization class.
    
    Class designed to convert sparse bit-vectors into a single block of memory
    ready for file or database storage or network transfer.
    
    Reuse of this class for multiple serializations may offer some performance
    advantage.
    
    @ingroup bvserial 
*/
template<class BV>
class serializer
{
public:
    typedef BV                                                bvector_type;
    typedef typename bvector_type::allocator_type             allocator_type;
    typedef typename bvector_type::blocks_manager_type        blocks_manager_type;
    typedef typename bvector_type::statistics                 statistics_type;

    typedef byte_buffer<allocator_type> buffer;
public:
    /**
        Construct serializer
        
        \param alloc - memory allocator
        \param temp_block - temporary block for various operations
               (if NULL it will be allocated and managed by serializer class)
        Temp block is used as a scratch memory during serialization,
        use of external temp block allows to avoid unnecessary re-allocations.
     
        Temp block attached is not owned by the class and NOT deallocated on
        destruction.
    */
    serializer(const allocator_type&   alloc  = allocator_type(),
              bm::word_t*  temp_block = 0);
    
    serializer(bm::word_t*  temp_block);

    ~serializer();

    /**
        Set compression level. Higher compression takes more time to process.
        @param clevel - compression level (0-4)
    */
    void set_compression_level(unsigned clevel);

    /**
        Get compression level
    */
    unsigned get_compression_level() const;
    
    /**
        Bitvector serilization into memory block
        
        @param bv - input bitvector
        @param buf - out buffer (pre-allocated)
           No range checking is done in this method. 
           It is responsibility of caller to allocate sufficient 
           amount of memory using information from calc_stat() function.        
        
        @param buf_size - size of the output buffer
        
       @return Size of serialization block.
       @sa calc_stat     
    */
    unsigned serialize(const BV& bv, 
                       unsigned char* buf, size_t buf_size);
    
    /**
        Bitvector serilization into buffer object (it gets resized automatically)
     
        @param bv       - input bitvector
        @param buf      - output buffer object
        @param bv_stat  - input (optional) bit-vector statistics object
                          if NULL, serizlize will compute statistics
    */
    void serialize(const BV& bv, typename serializer<BV>::buffer& buf, const statistics_type* bv_stat);

    
    /**
        Set GAP length serialization (serializes GAP levels of the original vector)
                
        @param value - when TRUE serialized vector includes GAP levels parameters
    */
    void gap_length_serialization(bool value);
    /**
        Set byte-order serialization (for cross platform compatibility)
        
        @param value - TRUE serialization format includes byte-order marker
    */
    void byte_order_serialization(bool value);

protected:
    /**
        Encode serialization header information
    */
    void encode_header(const BV& bv, bm::encoder& enc);
    
    /**
        Encode GAP block
    */
    void encode_gap_block(bm::gap_word_t* gap_block, bm::encoder& enc);

    /**
        Encode GAP block with Elias Gamma coder
    */
    void gamma_gap_block(bm::gap_word_t* gap_block, bm::encoder& enc);

    /**
        Encode GAP block as delta-array with Elias Gamma coder
    */
    void gamma_gap_array(const bm::gap_word_t* gap_block, 
                         unsigned              arr_len, 
                         bm::encoder&          enc,
                         bool                  inverted = false);

    /**
        Encode BIT block with repeatable runs of zeroes
    */
    void encode_bit_interval(const bm::word_t* blk, 
                             bm::encoder&      enc,
                             unsigned          size_control);

private:
    serializer(const serializer&);
    serializer& operator=(const serializer&);

private:

    typedef bm::bit_out<bm::encoder>  bit_out_type;
    typedef bm::gamma_encoder<bm::gap_word_t, bit_out_type> gamma_encoder_func;

private:
    allocator_type alloc_;
    bool           gap_serial_;
    bool           byte_order_serial_;
    bm::word_t*    temp_block_;
    unsigned       compression_level_;
    bool           own_temp_block_;
};

/**
    Base deserialization class
    \ingroup bvserial
*/
template<class DEC>
class deseriaizer_base
{
public:
    typedef DEC decoder_type;
protected:
    deseriaizer_base(){}

    /// Read GAP block from the stream
    void read_gap_block(decoder_type&   decoder, 
                        unsigned        block_type, 
                        bm::gap_word_t* dst_block,
                        bm::gap_word_t& gap_head);

	/// Read list of bit ids
	///
	/// @return number of ids
	unsigned read_id_list(decoder_type&   decoder, 
                          unsigned        block_type, 
                          bm::gap_word_t* dst_arr);

protected:
    bm::gap_word_t   id_array_[bm::gap_equiv_len * 2];
};

/**
    Deserializer for bit-vector
    \ingroup bvserial 
*/
template<class BV, class DEC>
class deserializer : protected deseriaizer_base<DEC>
{
public:
    typedef BV bvector_type;
    typedef typename deseriaizer_base<DEC>::decoder_type decoder_type;
public:
    deserializer() : temp_block_(0) {}
    
    unsigned deserialize(bvector_type&        bv, 
                         const unsigned char* buf, 
                         bm::word_t*          temp_block);
protected:
   typedef typename BV::blocks_manager_type blocks_manager_type;
   typedef typename BV::allocator_type allocator_type;

protected:
   void deserialize_gap(unsigned char btype, decoder_type& dec, 
                        bvector_type&  bv, blocks_manager_type& bman,
                        unsigned i,
                        bm::word_t* blk);
protected:
    bm::gap_word_t   gap_temp_block_[bm::gap_equiv_len * 4];
    bm::word_t*      temp_block_;
};


/**
    Iterator to walk forward the serialized stream.

    \internal
    \ingroup bvserial 
*/
template<class BV, class SerialIterator>
class iterator_deserializer
{
public:
    typedef BV              bvector_type;
    typedef SerialIterator  serial_iterator_type;
public:
    static
    unsigned deserialize(bvector_type&         bv, 
                         serial_iterator_type& sit, 
                         bm::word_t*           temp_block,
                         set_operation         op = bm::set_OR,
                         bool                  exit_on_one = false);

    /// experimental 3 way deserialization 
    /// target = mask %OPERATION% BLOB
    ///
    static
    void deserialize(bvector_type&         bv_target,
                     const bvector_type&   bv_mask,
                     serial_iterator_type& sit, 
                     bm::word_t*           temp_block,
                     set_operation         op);

private:
    typedef typename BV::blocks_manager_type blocks_manager_type;

    /// load data from the iterator of type "id list"
    static
    void load_id_list(bvector_type&         bv, 
                      serial_iterator_type& sit,
                      unsigned              id_count,
                      bool                  set_clear);

    /// Finalize the deserialization (zero target vector tail or bit-count tail)
    static
    unsigned finalize_target_vector(blocks_manager_type& bman,
                                    set_operation        op,
                                    unsigned             bv_block_idx);

    /// Process (obsolete) id-list serialization format
    static
    unsigned process_id_list(bvector_type&         bv, 
                             serial_iterator_type& sit,
                             set_operation         op);


};

/*!
    @brief Serialization stream iterator

    Iterates blocks and control tokens of serialized bit-stream
    \ingroup bvserial 
*/
template<class DEC>
class serial_stream_iterator : protected deseriaizer_base<DEC>
{
public:
    typedef typename deseriaizer_base<DEC>::decoder_type decoder_type;
public:
    serial_stream_iterator(const unsigned char* buf);

    /// serialized bitvector size
    unsigned bv_size() const { return bv_size_; }

    /// Returns true if end of bit-stream reached 
    bool is_eof() const { return end_of_stream_; }

    /// get next block
    void next();

	/// skip all zero or all-one blocks
	void skip_mono_blocks();

    /// read bit block, using logical operation
    unsigned get_bit_block(bm::word_t*       dst_block, 
                           bm::word_t*       tmp_block,
                           set_operation     op);


    /// Read gap block data (with head)
    void get_gap_block(bm::gap_word_t* dst_block);

    /// Return current decoder size
    unsigned dec_size() const { return decoder_.size(); }

    /// Get low level access to the decoder (use carefully)
    decoder_type& decoder() { return decoder_; }

    /// iterator is a state machine, this enum encodes 
    /// its key value
    ///
    enum iterator_state 
    {
        e_unknown = 0,
        e_list_ids,     ///< plain int array
        e_blocks,       ///< stream of blocks
        e_zero_blocks,  ///< one or more zero bit blocks
        e_one_blocks,   ///< one or more all-1 bit blocks
        e_bit_block,    ///< one bit block
        e_gap_block     ///< one gap block

    };

    /// Returns iterator internal state
    iterator_state state() const { return this->state_; }

    iterator_state get_state() const { return this->state_; }
    /// Number of ids in the inverted list (valid for e_list_ids)
    unsigned get_id_count() const { return this->id_cnt_; }

    /// Get last id from the id list
    bm::id_t get_id() const { return this->last_id_; }

    /// Get current block index 
    unsigned block_idx() const { return this->block_idx_; }

public:
    /// member function pointer for bitset-bitset get operations
    /// 
    typedef 
        unsigned (serial_stream_iterator<DEC>::*get_bit_func_type)
                                                (bm::word_t*,bm::word_t*);

    unsigned 
    get_bit_block_ASSIGN(bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_OR    (bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_AND   (bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_SUB   (bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_XOR   (bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT (bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT_AND(bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT_OR(bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT_XOR(bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT_SUB_AB(bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT_SUB_BA(bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT_A(bm::word_t* dst_block, bm::word_t* tmp_block);
    unsigned 
    get_bit_block_COUNT_B(bm::word_t* dst_block, bm::word_t* tmp_block)
    {
        return get_bit_block_COUNT(dst_block, tmp_block);
    }

    /// Get array of bits out of the decoder into bit block
    /// (Converts inverted list into bits)
    /// Returns number of words (bits) being read
    unsigned get_arr_bit(bm::word_t* dst_block, 
                         bool clear_target=true);

	/// Get current block type
	unsigned get_block_type() const { return block_type_; }

	unsigned get_bit();

protected:
    get_bit_func_type  bit_func_table_[bm::set_END];

    decoder_type       decoder_;
    bool               end_of_stream_;
    unsigned           bv_size_;
    iterator_state     state_;
    unsigned           id_cnt_;  ///< Id counter for id list
    bm::id_t           last_id_; ///< Last id from the id list
    gap_word_t         glevels_[bm::gap_levels]; ///< GAP levels

    unsigned           block_type_;     ///< current block type
    unsigned           block_idx_;      ///< current block index
    unsigned           mono_block_cnt_; ///< number of 0 or 1 blocks

    gap_word_t         gap_head_;
};

/**
    Deserializer, performs logical operations between bit-vector and
    serialized bit-vector. This utility class potentially provides faster
    and/or more memory efficient operation than more conventional deserialization
    into memory bvector and then logical operation

    \ingroup bvserial 
*/
template<class BV>
class operation_deserializer
{
public:
    typedef BV bvector_type;
public:
    /**
    \brief Deserialize bvector using buffer as set operation argument
    
    \param bv - target bvector
    \param buf - serialized buffer as a logical argument
    \param temp_block - temporary block to avoid re-allocations
    \param op - set algebra operation (default: OR)
    \param exit_on_one - quick exit if set operation found some result
    
    \return bitcount
    */
    static
    unsigned deserialize(bvector_type&        bv, 
                         const unsigned char* buf, 
                         bm::word_t*          temp_block,
                         set_operation        op = bm::set_OR,
                         bool                 exit_on_one = false ///<! exit early if any one are found
                         );
private:
    /** experimental 3-way deserializator TARGET = MASK (OR/AND/XOR) BUF
    \param bv_target - target bvector
    \param bv_mask - mask bvector (MASK)
    \param buf - buffer argument
    \param temp_block - operational temporary block to avoid re-allocations
    \param op - logical operation
    */
    static
    void deserialize(bvector_type&        bv_target,
                     const bvector_type&  bv_mask,
                     const unsigned char* buf, 
                     bm::word_t*          temp_block,
                     set_operation        op);

private:
    typedef 
        typename BV::blocks_manager_type               blocks_manager_type;
    typedef 
        serial_stream_iterator<bm::decoder>            serial_stream_current;
    typedef 
        serial_stream_iterator<bm::decoder_big_endian> serial_stream_be;
    typedef 
        serial_stream_iterator<bm::decoder_little_endian> serial_stream_le;

};





//---------------------------------------------------------------------

template<class BV>
serializer<BV>::serializer(const allocator_type&   alloc,
                           bm::word_t*             temp_block)
: alloc_(alloc),
  gap_serial_(false),
  byte_order_serial_(true),
  compression_level_(4)
{
    if (temp_block == 0)
    {
        temp_block_ = alloc_.alloc_bit_block();
        own_temp_block_ = true;
    }
    else
    {
        temp_block_ = temp_block;
        own_temp_block_ = false;
    }
}

template<class BV>
serializer<BV>::serializer(bm::word_t*    temp_block)
: alloc_(allocator_type()),
  gap_serial_(false),
  byte_order_serial_(true),
  compression_level_(4)
{
    if (temp_block == 0)
    {
        temp_block_ = alloc_.alloc_bit_block();
        own_temp_block_ = true;
    }
    else
    {
        temp_block_ = temp_block;
        own_temp_block_ = false;
    }
}


template<class BV>
void serializer<BV>::set_compression_level(unsigned clevel)
{
    compression_level_ = clevel;
}

template<class BV>
unsigned serializer<BV>::get_compression_level() const
{
    return compression_level_;
}

template<class BV>
serializer<BV>::~serializer()
{
    if (own_temp_block_)
        alloc_.free_bit_block(temp_block_);
}


template<class BV>
void serializer<BV>::gap_length_serialization(bool value)
{
    gap_serial_ = value;
}

template<class BV>
void serializer<BV>::byte_order_serialization(bool value)
{
    byte_order_serial_ = value;
}

template<class BV>
void serializer<BV>::encode_header(const BV& bv, bm::encoder& enc)
{
    const blocks_manager_type& bman = bv.get_blocks_manager();

    unsigned char header_flag = 0;
    if (bv.size() == bm::id_max) // no dynamic resize
        header_flag |= BM_HM_DEFAULT;
    else 
        header_flag |= BM_HM_RESIZE;

    if (!byte_order_serial_) 
        header_flag |= BM_HM_NO_BO;

    if (!gap_serial_) 
        header_flag |= BM_HM_NO_GAPL;

    enc.put_8(header_flag);

    if (byte_order_serial_)
    {
        ByteOrder bo = globals<true>::byte_order();
        enc.put_8((unsigned char)bo);
    }

    // keep GAP levels information
    if (gap_serial_)
    {
        enc.put_16(bman.glen(), bm::gap_levels);
    }

    // save size (only if bvector has been down-sized)
    if (header_flag & BM_HM_RESIZE) 
    {
        enc.put_32(bv.size());
    }
    
}

template<class BV>
void serializer<BV>::gamma_gap_block(bm::gap_word_t* gap_block, bm::encoder& enc)
{
    unsigned len = gap_length(gap_block);

    // Use Elias Gamma encoding 
    if (len > 6 && (compression_level_ > 3)) 
    {
        encoder::position_type enc_pos0 = enc.get_pos();
        {
            bit_out_type bout(enc);
            gamma_encoder_func gamma(bout);

            enc.put_8(set_block_gap_egamma);		
            enc.put_16(gap_block[0]);

            for_each_dgap(gap_block, gamma);        
        }

        // evaluate gamma coding efficiency
        encoder::position_type enc_pos1 = enc.get_pos();
        unsigned gamma_size = (unsigned)(enc_pos1 - enc_pos0);        
        if (gamma_size > (len-1)*sizeof(gap_word_t))
        {
            enc.set_pos(enc_pos0);
        }
        else
        {
            return;
        }
    }

    // save as plain GAP block 
    enc.put_8(set_block_gap);
    enc.put_16(gap_block, len-1);
}

template<class BV>
void serializer<BV>::gamma_gap_array(const bm::gap_word_t* gap_array, 
                                     unsigned              arr_len, 
                                     bm::encoder&          enc,
                                     bool                  inverted)
{
    if (compression_level_ > 3 && arr_len > 25)
    {        
        encoder::position_type enc_pos0 = enc.get_pos();
        {
            bit_out_type bout(enc);

            enc.put_8(
                inverted ? set_block_arrgap_egamma_inv 
                         : set_block_arrgap_egamma);

            bout.gamma(arr_len);

            gap_word_t prev = gap_array[0];
            bout.gamma(prev + 1);

            for (unsigned i = 1; i < arr_len; ++i)
            {
                gap_word_t curr = gap_array[i];
                bout.gamma(curr - prev);
                prev = curr;
            }
        }

        encoder::position_type enc_pos1 = enc.get_pos();
        unsigned gamma_size = (unsigned)(enc_pos1 - enc_pos0);            
        if (gamma_size > (arr_len)*sizeof(gap_word_t))
        {
            enc.set_pos(enc_pos0);
        }
        else
        {
            return;
        }
    }

    // save as an plain array
    enc.put_prefixed_array_16(inverted ? set_block_arrgap_inv : set_block_arrgap, 
                              gap_array, arr_len, true);
}


template<class BV>
void serializer<BV>::encode_gap_block(bm::gap_word_t* gap_block, bm::encoder& enc)
{
    if (compression_level_ > 2)
    {
        gap_word_t*  gap_temp_block = (gap_word_t*) temp_block_;    
        gap_word_t arr_len;

        unsigned bc = gap_bit_count_unr(gap_block);
        if (bc == 1)
        {
            arr_len = gap_convert_to_arr(gap_temp_block,
                                         gap_block,
                                         bm::gap_equiv_len-10);
            BM_ASSERT(arr_len == 1);
	        enc.put_8(set_block_bit_1bit);				
	        enc.put_16(gap_temp_block[0]);
            return;
        }

        unsigned len = gap_length(gap_block);
        bool invert, use_array;
        invert = use_array = false;
        
        if (bc < len-1) 
        {
            use_array = true;
        }
        else  // inverted array is a better alternative?
        {
            unsigned inverted_bc = bm::gap_max_bits - bc;
            if (inverted_bc < len-1)
            {
                use_array = invert = true;
            }
        }
        if (use_array)
        {
            arr_len = gap_convert_to_arr(gap_temp_block,
                                         gap_block,
                                         bm::gap_equiv_len-10,
                                         invert);
            if (arr_len)
            {
                gamma_gap_array(gap_temp_block, arr_len, enc, invert);
                return;
            }
        }
    }

    gamma_gap_block(gap_block, enc);
}

template<class BV>
void serializer<BV>::encode_bit_interval(const bm::word_t* blk, 
                                         bm::encoder&      enc,
                                         unsigned          //size_control
                                         )
{
    enc.put_8(set_block_bit_0runs);
    enc.put_8((blk[0]==0) ? 0 : 1); // encode start 

    unsigned i,j;//,k;

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk[i] == 0)
        {
            // scan fwd to find 0 island length
            for (j = i+1; j < bm::set_block_size; ++j)
            {
                if (blk[j] != 0)
                    break;
            }
            enc.put_16((gap_word_t)(j-i)); 
            i = j - 1;
        }
        else
        {
            // scan fwd to find non-0 island length
            for (j = i+1; j < bm::set_block_size; ++j)
            {
                if (blk[j] == 0)
                {
                    // look ahead to identify and ignore short 0-run
                    if (((j+1 < bm::set_block_size) && blk[j+1]) ||
                        ((j+2 < bm::set_block_size) && blk[j+2])
                       )
                    {
                        ++j; // skip zero word
                        continue;
                    }
                    break;
                }
            }
            enc.put_16((gap_word_t)(j-i)); 
            // stream all bit-words now
            BM_ASSERT(i < j);
            enc.put_32(blk + i, j - i);
            i = j - 1;
        }
    }
}

template<class BV>
void serializer<BV>::serialize(const BV& bv,
                               typename serializer<BV>::buffer& buf,
                               const statistics_type* bv_stat)
{
    statistics_type stat;
    if (!bv_stat)
    {
        bv.calc_stat(&stat);
        bv_stat = &stat;
    }
    
    buf.resize(bv_stat->max_serialize_mem);
    
    unsigned slen = this->serialize(bv, buf.data(), buf.size());
    BM_ASSERT(slen <= buf.size()); // or we have a BIG problem with prediction
    BM_ASSERT(slen);
    
    buf.resize(slen);
}


template<class BV>
unsigned serializer<BV>::serialize(const BV& bv, 
                                   unsigned char* buf, size_t buf_size)
{
    BM_ASSERT(temp_block_);
    
    const blocks_manager_type& bman = bv.get_blocks_manager();

    gap_word_t*  gap_temp_block = (gap_word_t*) temp_block_;
    
    bm::encoder enc(buf, buf_size);  // create the encoder
    encode_header(bv, enc);

    unsigned i,j;


    // save blocks.
    for (i = 0; i < bm::set_total_blocks; ++i)
    {
        bm::word_t* blk = bman.get_block(i);
        // -----------------------------------------
        // Empty or ONE block serialization

        bool flag;
        flag = bm::check_block_zero(blk, false/*shallow check*/);
        if (flag)
        {
        zero_block:
            unsigned next_nb = bman.find_next_nz_block(i+1, false);
            if (next_nb == bm::set_total_blocks) // no more blocks
            {
                enc.put_8(set_block_azero);
                return enc.size();
            }
            unsigned nb = next_nb - i;
            
            if (nb > 1 && nb < 128)
            {
                // special (but frequent) case -- GAP delta fits in 7bits
                unsigned char c = (unsigned char)((1u << 7) | nb);
                enc.put_8(c);
            }
            else 
            {
                SER_NEXT_GRP(enc, nb, set_block_1zero, 
                                      set_block_8zero, 
                                      set_block_16zero, 
                                      set_block_32zero) 
            }
            i = next_nb - 1;
            continue;
        }
        else
        {
            flag = bm::check_block_one(blk, false);
            if (flag)
            {
                // Look ahead for similar blocks
                for(j = i+1; j < bm::set_total_blocks; ++j)
                {
                   bm::word_t* blk_next = bman.get_block(j);
                   if (flag != bm::check_block_one(blk_next, false))
                       break;
                }
                if (j == bm::set_total_blocks)
                {
                    enc.put_8(set_block_aone);
                    break;
                }
                else
                {
                   unsigned nb = j - i;
                   SER_NEXT_GRP(enc, nb, set_block_1one, 
                                         set_block_8one, 
                                         set_block_16one, 
                                         set_block_32one) 
                }
                i = j - 1;
                continue;
            }
        }

        // ------------------------------
        // GAP serialization

        if (BM_IS_GAP(blk))
        {
            gap_word_t* gblk = BMGAP_PTR(blk);
            encode_gap_block(gblk, enc);
            continue;
        }
                
        // ----------------------------------------------
        // BIT BLOCK serialization

        {
        if (compression_level_ <= 1)
        {
            enc.put_prefixed_array_32(set_block_bit, blk, bm::set_block_size);
            continue;            
        }

        // compute bit-block statistics: bit-count and number of GAPS
        unsigned block_bc = 0;
        bm::id_t bit_gaps = 
            bm::bit_block_calc_count_change(blk, blk + bm::set_block_size, &block_bc);
        unsigned block_bc_inv = bm::gap_max_bits - block_bc;
        switch (block_bc)
        {
        case 1: // corner case: only 1 bit on
            {
                bm::id_t bit_idx = 0;
                bm::bit_find_in_block(blk, bit_idx, &bit_idx);
                enc.put_8(set_block_bit_1bit); enc.put_16(bm::short_t(bit_idx));
                continue;
            }
        case 0: goto zero_block; // empty block
        default:
            break;
        }
       
       
        // compute alternative representation sizes
        //
        unsigned arr_block_size = unsigned(sizeof(gap_word_t) + (block_bc * sizeof(gap_word_t)));
        unsigned arr_block_size_inv = unsigned(sizeof(gap_word_t) + (block_bc_inv * sizeof(gap_word_t)));
        unsigned gap_block_size = unsigned(sizeof(gap_word_t) + ((bit_gaps+1) * sizeof(gap_word_t)));
        unsigned interval_block_size;
        interval_block_size = bit_count_nonzero_size(blk, bm::set_block_size);
        
        bool inverted = false;

        if (arr_block_size_inv < arr_block_size &&
            arr_block_size_inv < gap_block_size &&
            arr_block_size_inv < interval_block_size)
        {
            inverted = true;
            goto bit_as_array;
        }

        // if interval representation is not a good alternative
        if ((interval_block_size > arr_block_size) || 
            (interval_block_size > gap_block_size))
        {
            if (gap_block_size < (bm::gap_equiv_len-64) &&
                gap_block_size < arr_block_size)
            {
                unsigned len = bit_convert_to_gap(gap_temp_block, 
                                                  blk, 
                                                  bm::gap_max_bits, 
                                                  bm::gap_equiv_len-64);
                if (len) // save as GAP
                {
                    gamma_gap_block(gap_temp_block, enc);
                    continue;
                }
            }
            
            if (arr_block_size < ((bm::gap_equiv_len-64) * sizeof(gap_word_t)))
            {
            bit_as_array:
                gap_word_t arr_len;
                unsigned mask = inverted ? ~0u : 0u;
                arr_len = bit_convert_to_arr(gap_temp_block, 
                                             blk, 
                                             bm::gap_max_bits, 
                                             bm::gap_equiv_len-64,
                                             mask);
                if (arr_len)
                {
                    gamma_gap_array(gap_temp_block, arr_len, enc, inverted);
                    continue;
                }
                
            }
            // full bit-block
            enc.put_prefixed_array_32(set_block_bit, blk, bm::set_block_size);
            continue;            
        }
        
        // if interval block is a winner
        // it needs to have a compelling advantage of 25% over bit block
        //
        unsigned threashold_block_size =
            bm::set_block_size * sizeof(bm::word_t);
        threashold_block_size -= threashold_block_size / 4;
            
        if (interval_block_size < arr_block_size &&
            interval_block_size < gap_block_size &&
            interval_block_size < (bm::set_block_size * sizeof(bm::word_t))
            )
        {
            encode_bit_interval(blk, enc, interval_block_size);
            continue;
        }
        
        if (gap_block_size < bm::gap_equiv_len &&
            gap_block_size < arr_block_size)
        {
            unsigned len = bit_convert_to_gap(gap_temp_block, 
                                              blk, 
                                              bm::gap_max_bits, 
                                              bm::gap_equiv_len-64);
            if (len) // save as GAP
            {
                gamma_gap_block(gap_temp_block, enc);
                continue;
            }
        }
        
         
        // if array is best
        if (arr_block_size < bm::gap_equiv_len-64)
        {
            goto bit_as_array;
        }
        // full bit-block
        enc.put_prefixed_array_32(set_block_bit, blk, bm::set_block_size);
        continue;            
        }
    }

    enc.put_8(set_block_end);

    unsigned encoded_size = enc.size();
    return encoded_size;

}



/// Bit mask flags for serialization algorithm
/// \ingroup bvserial 
enum serialization_flags {
    BM_NO_BYTE_ORDER = 1,       ///< save no byte-order info (save some space)
    BM_NO_GAP_LENGTH = (1 << 1) ///< save no GAP info (save some space)
};

/*!
   \brief Saves bitvector into memory.

   Function serializes content of the bitvector into memory.
   Serialization adaptively uses compression(variation of GAP encoding) 
   when it is benefitial. 
   
   \param bv - source bvecor
   \param buf - pointer on target memory area. No range checking in the
   function. It is responsibility of programmer to allocate sufficient 
   amount of memory using information from calc_stat function.

   \param temp_block - pointer on temporary memory block. Cannot be 0; 
   If you want to save memory across multiple bvectors
   allocate temporary block using allocate_tempblock and pass it to 
   serialize.
   (Serialize does not deallocate temp_block.)

   \param serialization_flags
   Flags controlling serilization (bit-mask) 
   (use OR-ed serialization flags)

   \ingroup bvserial 

   \return Size of serialization block.
   \sa calc_stat, serialization_flags

*/
/*!
 Serialization format:
 <pre>

 | HEADER | BLOCKS |

 Header structure:
   BYTE : Serialization header (bit mask of BM_HM_*)
   BYTE : Byte order ( 0 - Big Endian, 1 - Little Endian)
   INT16: Reserved (0)
   INT16: Reserved Flags (0)

 </pre>
*/
template<class BV>
unsigned serialize(const BV& bv, 
                   unsigned char* buf, 
                   bm::word_t*    temp_block = 0,
                   unsigned       serialization_flags = 0)
{
    bm::serializer<BV> bv_serial(bv.get_allocator(), temp_block);
    if (serialization_flags & BM_NO_BYTE_ORDER) 
        bv_serial.byte_order_serialization(false);
        
    if (serialization_flags & BM_NO_GAP_LENGTH) 
        bv_serial.gap_length_serialization(false);
    else
        bv_serial.gap_length_serialization(true);

    bv_serial.set_compression_level(4);
    
    return bv_serial.serialize(bv, buf, 0);
}

/*!
   @brief Saves bitvector into memory.
   Allocates temporary memory block for bvector.
 
   \param bv - source bvecor
   \param buf - pointer on target memory area. No range checking in the
   function. It is responsibility of programmer to allocate sufficient 
   amount of memory using information from calc_stat function.

   \param serialization_flags
   Flags controlling serilization (bit-mask) 
   (use OR-ed serialization flags)
 
   \ingroup bvserial
*/
template<class BV>
unsigned serialize(BV& bv, 
                   unsigned char* buf, 
                   unsigned  serialization_flags=0)
{
    return bm::serialize(bv, buf, 0, serialization_flags);
}



/*!
    @brief Bitvector deserialization from memory.

    @param bv - target bvector
    @param buf - pointer on memory which keeps serialized bvector
    @param temp_block - pointer on temporary block, 
            if NULL bvector allocates own.
    @return Number of bytes consumed by deserializer.

    Function desrializes bitvector from memory block containig results
    of previous serialization. Function does not remove bits 
    which are currently set. Effectively it means OR logical operation 
    between current bitset and previously serialized one.

    @ingroup bvserial
*/
template<class BV>
unsigned deserialize(BV& bv, 
                     const unsigned char* buf, 
                     bm::word_t* temp_block=0)
{
    ByteOrder bo_current = globals<true>::byte_order();

    bm::decoder dec(buf);
    unsigned char header_flag = dec.get_8();
    ByteOrder bo = bo_current;
    if (!(header_flag & BM_HM_NO_BO))
    {
        bo = (bm::ByteOrder) dec.get_8();
    }

    if (bo_current == bo)
    {
        deserializer<BV, bm::decoder> deserial;
        return deserial.deserialize(bv, buf, temp_block);
    }
    switch (bo_current) 
    {
    case BigEndian:
        {
        deserializer<BV, bm::decoder_big_endian> deserial;
        return deserial.deserialize(bv, buf, temp_block);
        }
    case LittleEndian:
        {
        deserializer<BV, bm::decoder_little_endian> deserial;
        return deserial.deserialize(bv, buf, temp_block);
        }
    default:
        BM_ASSERT(0);
    };
    return 0;
}

template<class DEC>
unsigned deseriaizer_base<DEC>::read_id_list(decoder_type&   decoder, 
		    								 unsigned        block_type, 
				   					         bm::gap_word_t* dst_arr)
{
    typedef bit_in<DEC> bit_in_type;

	gap_word_t len = 0;

    switch (block_type)
    {
    case set_block_bit_1bit:
        dst_arr[0] = decoder.get_16();
		len = 1;
		break;
    case set_block_arrgap:
    case set_block_arrgap_inv:
        len = decoder.get_16();
        decoder.get_16(dst_arr, len);
		break;
    case set_block_arrgap_egamma:
    case set_block_arrgap_egamma_inv:
        {
            bit_in_type bin(decoder);
            len = (gap_word_t)bin.gamma();
            gap_word_t prev = 0;
            for (gap_word_t k = 0; k < len; ++k)
            {
                gap_word_t bit_idx = (gap_word_t)bin.gamma();
                if (k == 0) --bit_idx; // TODO: optimization
                bit_idx = (gap_word_t)(bit_idx + prev);
                prev = bit_idx;
				dst_arr[k] = bit_idx;
            } // for
        }
        break;
    default:
        BM_ASSERT(0);
    }
	return len;
}


template<class DEC>
void deseriaizer_base<DEC>::read_gap_block(decoder_type&   decoder, 
                                           unsigned        block_type, 
                                           bm::gap_word_t* dst_block,
                                           bm::gap_word_t& gap_head)
{
    typedef bit_in<DEC> bit_in_type;

    switch (block_type)
    {
    case set_block_gap:
        {
            unsigned len = gap_length(&gap_head);
            --len;
            *dst_block = gap_head;
            decoder.get_16(dst_block+1, len - 1);
            dst_block[len] = gap_max_bits - 1;
        }
        break;
    case set_block_bit_1bit:
        {
			gap_set_all(dst_block, bm::gap_max_bits, 0);
            gap_word_t bit_idx = decoder.get_16();
			gap_add_value(dst_block, bit_idx);
        }
        break;
    case set_block_arrgap:
    case set_block_arrgap_inv:
        {
            gap_set_all(dst_block, bm::gap_max_bits, 0);
            gap_word_t len = decoder.get_16();

            for (gap_word_t k = 0; k < len; ++k)
            {
                gap_word_t bit_idx = decoder.get_16();
				gap_add_value(dst_block, bit_idx);
            } // for
        }
        break;
    case set_block_arrgap_egamma:
    case set_block_arrgap_egamma_inv:
        {
        	unsigned arr_len = read_id_list(decoder, block_type, id_array_);
            dst_block[0] = 0;
            //unsigned gap_len = 
                gap_set_array(dst_block, id_array_, arr_len);
        }
        break;
    case set_block_gap_egamma:
        {
        unsigned len = (gap_head >> 3);
        --len;
        // read gamma GAP block into a dest block
        {
            *dst_block = gap_head;
            gap_word_t* gap_data_ptr = dst_block + 1;

            bit_in_type bin(decoder);
            {
				gap_word_t v = (gap_word_t)bin.gamma();
                gap_word_t gap_sum = *gap_data_ptr = (gap_word_t)(v - 1);
                for (unsigned i = 1; i < len; ++i)
                {					
                    v = (gap_word_t)bin.gamma();
                    gap_sum = (gap_word_t)(gap_sum + v);
                    *(++gap_data_ptr) = gap_sum;
                }
                dst_block[len+1] = bm::gap_max_bits - 1;
            }
        }

        }
        break;        
    default:
        BM_ASSERT(0);
    }

    if (block_type == set_block_arrgap_egamma_inv || 
        block_type == set_block_arrgap_inv)
    {
        gap_invert(dst_block);
    }
}


template<class BV, class DEC>
void 
deserializer<BV, DEC>::deserialize_gap(unsigned char btype, decoder_type& dec, 
                                       bvector_type&  bv, blocks_manager_type& bman,
                                       unsigned i,
                                       bm::word_t* blk)
{
    //typedef bit_in<DEC> bit_in_type;
    gap_word_t gap_head; 

    switch (btype)
    {
    case set_block_gap: 
    case set_block_gapbit:
    {
        gap_head = (gap_word_t)
            (sizeof(gap_word_t) == 2 ? dec.get_16() : dec.get_32());

        unsigned len = gap_length(&gap_head);
        int level = gap_calc_level(len, bman.glen());
        --len;
        if (level == -1)  // Too big to be GAP: convert to BIT block
        {
            *gap_temp_block_ = gap_head;
            dec.get_16(gap_temp_block_+1, len - 1);
            gap_temp_block_[len] = gap_max_bits - 1;

            if (blk == 0)  // block does not exist yet
            {
                blk = bman.get_allocator().alloc_bit_block();
                bman.set_block(i, blk);
                gap_convert_to_bitset(blk, gap_temp_block_);                
            }
            else  // We have some data already here. Apply OR operation.
            {
                gap_convert_to_bitset(temp_block_, 
                                      gap_temp_block_);
                bv.combine_operation_with_block(i, 
                                                temp_block_, 
                                                0, 
                                                BM_OR);
            }
            return;
        } // level == -1

        set_gap_level(&gap_head, level);

        if (blk == 0)
        {
            BM_ASSERT(level >= 0);
            gap_word_t* gap_blk = 
              bman.get_allocator().alloc_gap_block(unsigned(level), bman.glen());
            gap_word_t* gap_blk_ptr = BMGAP_PTR(gap_blk);
            *gap_blk_ptr = gap_head;
            bm::set_gap_level(gap_blk_ptr, level);
            blk = bman.set_block(i, (bm::word_t*)BMPTR_SETBIT0(gap_blk));
            BM_ASSERT(blk == 0);
            
            dec.get_16(gap_blk + 1, len - 1);
            gap_blk[len] = bm::gap_max_bits - 1;
        }
        else // target block exists
        {
            // read GAP block into a temp memory and perform OR
            *gap_temp_block_ = gap_head;
            dec.get_16(gap_temp_block_ + 1, len - 1);
            gap_temp_block_[len] = bm::gap_max_bits - 1;
            break;
        }
        return;
    }
    case set_block_arrgap: 
    case set_block_arrgap_egamma:
        {
        	unsigned arr_len = this->read_id_list(dec, btype, this->id_array_);
            gap_temp_block_[0] = 0; // reset unused bits in gap header
            unsigned gap_len =
              gap_set_array(gap_temp_block_, this->id_array_, arr_len);

              int level = gap_calc_level(gap_len, bman.glen());
              if (level == -1)  // Too big to be GAP: convert to BIT block
              {
                  gap_convert_to_bitset(temp_block_, gap_temp_block_);
                  bv.combine_operation_with_block(i,
                                                  temp_block_,
                                                  0,
                                                  BM_OR);
                  return;
              }

            break;
        }
    case set_block_gap_egamma:            
        gap_head = (gap_word_t)
            (sizeof(gap_word_t) == 2 ? dec.get_16() : dec.get_32());
    case set_block_arrgap_egamma_inv:
    case set_block_arrgap_inv:
        this->read_gap_block(dec, btype, gap_temp_block_, gap_head);
        break;
    default:
        BM_ASSERT(0);
    }

    bv.combine_operation_with_block(i, 
                                   (bm::word_t*)gap_temp_block_, 
                                    1, 
                                    BM_OR);

}


template<class BV, class DEC>
unsigned deserializer<BV, DEC>::deserialize(bvector_type&        bv, 
                                            const unsigned char* buf,
                                            bm::word_t*          temp_block)
{
    blocks_manager_type& bman = bv.get_blocks_manager();
    if (!bman.is_init())
    {
        bman.init_tree();
    }

    bm::wordop_t* tmp_buf = 
        temp_block ? (bm::wordop_t*) temp_block 
                   : (bm::wordop_t*)bman.check_allocate_tempblock();

    temp_block_ = temp_block = (word_t*)tmp_buf;
    bm::strategy  strat = bv.get_new_blocks_strat();
    bv.set_new_blocks_strat(BM_GAP);

    decoder_type dec(buf);

    BM_SET_MMX_GUARD

    // Reading header

    unsigned char header_flag =  dec.get_8();
    if (!(header_flag & BM_HM_NO_BO))
    {
        /*ByteOrder bo = (bm::ByteOrder)*/dec.get_8();
    }

    if (header_flag & BM_HM_ID_LIST)
    {
        // special case: the next comes plain list of integers
        if (header_flag & BM_HM_RESIZE)
        {
            unsigned bv_size = dec.get_32();
            if (bv_size > bv.size())
            {
                bv.resize(bv_size);
            }
        }
        

        for (unsigned cnt = dec.get_32(); cnt; --cnt) {
            bm::id_t id = dec.get_32();
            bv.set(id);
        } // for
        // -1 for compatibility with other deserialization branches
        return dec.size()-1;
    }

    unsigned i;

    if (!(header_flag & BM_HM_NO_GAPL)) 
    {
        //gap_word_t glevels[bm::gap_levels];
        // read GAP levels information
        for (i = 0; i < bm::gap_levels; ++i)
        {
            /*glevels[i] =*/ dec.get_16();
        }
    }

    if (header_flag & (1 << 1))
    {
        unsigned bv_size = dec.get_32();
        if (bv_size > bv.size())
        {
            bv.resize(bv_size);
        }
    }

    unsigned char btype;
    unsigned nb;

    for (i = 0; i < bm::set_total_blocks; ++i)
    {
        btype = dec.get_8();
        bm::word_t* blk = bman.get_block(i);
        
        // pre-check if we have short zero-run packaging here
        //
        if (btype & (1 << 7))
        {
            nb = btype & ~(1 << 7);
            i += nb-1;
            continue;
        }        

        switch (btype)
        {
        case set_block_azero: 
        case set_block_end:
            i = bm::set_total_blocks;
            break;
        case set_block_1zero:
            continue;
        case set_block_8zero:
            nb = dec.get_8();
            i += nb-1;
            continue;
        case set_block_16zero:
            nb = dec.get_16();
            i += nb-1;
            continue;
        case set_block_32zero:
            nb = dec.get_32();
            i += nb-1;
            continue;
        case set_block_aone:
            for (;i < bm::set_total_blocks; ++i)
            {
                bman.set_block_all_set(i);
            }
            break;
        case set_block_1one:
            bman.set_block_all_set(i);
            continue;
        case set_block_8one:
            BM_SET_ONE_BLOCKS(dec.get_8());
            continue;
        case set_block_16one:
            BM_SET_ONE_BLOCKS(dec.get_16());
            continue;
        case set_block_32one:
            BM_SET_ONE_BLOCKS(dec.get_32());
            continue;
        case set_block_bit: 
        {
            if (blk == 0)
            {
                blk = bman.get_allocator().alloc_bit_block();
                bman.set_block(i, blk);
                dec.get_32(blk, bm::set_block_size);
                continue;                
            }
			
            dec.get_32(temp_block, bm::set_block_size);
            bv.combine_operation_with_block(i, 
                                            temp_block, 
                                            0, BM_OR);
			
            continue;
        }
        case set_block_bit_1bit:
        {
            unsigned bit_idx = dec.get_16();
            bit_idx += i * bm::bits_in_block; 
            bv.set_bit(bit_idx);
            continue;
        }
        case set_block_bit_0runs:
        {
            //TODO: optimization if block exists
            bit_block_set(temp_block, 0);

            unsigned char run_type = dec.get_8();
            for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
            {
                unsigned run_length = dec.get_16();
                if (run_type)
                {
                    unsigned run_end = j + run_length;
                    for (;j < run_end; ++j)
                    {
                        BM_ASSERT(j < bm::set_block_size);
                        temp_block[j] = dec.get_32();
                    }
                }
                else
                {
                    j += run_length;
                }
            } // for

            bv.combine_operation_with_block(i, 
                                            temp_block,
                                            0, BM_OR);            
            continue;
        }
        case set_block_bit_interval: 
        {
            unsigned head_idx, tail_idx;
            head_idx = dec.get_16();
            tail_idx = dec.get_16();

            if (blk == 0)
            {
                blk = bman.get_allocator().alloc_bit_block();
                bman.set_block(i, blk);
                for (unsigned k = 0; k < head_idx; ++k)
                {
                    blk[k] = 0;
                }
                dec.get_32(blk + head_idx, tail_idx - head_idx + 1);
                for (unsigned j = tail_idx + 1; j < bm::set_block_size; ++j)
                {
                    blk[j] = 0;
                }
                continue;
            }
            bit_block_set(temp_block, 0);
            dec.get_32(temp_block + head_idx, tail_idx - head_idx + 1);

            bv.combine_operation_with_block(i, 
                                            temp_block,
                                            0, BM_OR);
            continue;
        }
        case set_block_gap: 
        case set_block_gapbit:
        case set_block_arrgap:
        case set_block_gap_egamma:
        case set_block_arrgap_egamma:
        case set_block_arrgap_egamma_inv:
        case set_block_arrgap_inv:    
            deserialize_gap(btype, dec, bv, bman, i, blk);
            continue;
        case set_block_arrbit:
        {
            gap_word_t len = (gap_word_t)
                (sizeof(gap_word_t) == 2 ? dec.get_16() : dec.get_32());

            if (BM_IS_GAP(blk))
            {
                // convert from GAP cause generic bitblock is faster
                blk = bman.deoptimize_block(i);
            }
            else
            {
                if (blk == 0)  // block does not exists yet
                {
                    blk = bman.get_allocator().alloc_bit_block();
                    bm::bit_block_set(blk, 0);
                    bman.set_block(i, blk);
                }
            }

            // Get the array one by one and set the bits.
            for (unsigned k = 0; k < len; ++k)
            {
                gap_word_t bit_idx = dec.get_16();
				bm::set_bit(blk, bit_idx);
            }
            continue;
        }
        default:
            BM_ASSERT(0); // unknown block type
        } // switch
    } // for i

    bv.forget_count();
    bv.set_new_blocks_strat(strat);

    return dec.size();
}



template<class DEC>
serial_stream_iterator<DEC>::serial_stream_iterator(const unsigned char* buf)
  : decoder_(buf),
    end_of_stream_(false),
    bv_size_(0),
    state_(e_unknown),
    id_cnt_(0),
    block_idx_(0),
    mono_block_cnt_(0)
{
    ::memset(bit_func_table_, 0, sizeof(bit_func_table_));

    bit_func_table_[bm::set_AND] = 
        &serial_stream_iterator<DEC>::get_bit_block_AND;
    bit_func_table_[bm::set_ASSIGN] = 
        &serial_stream_iterator<DEC>::get_bit_block_ASSIGN;
    bit_func_table_[bm::set_OR]     = 
        &serial_stream_iterator<DEC>::get_bit_block_OR;
    bit_func_table_[bm::set_SUB] = 
        &serial_stream_iterator<DEC>::get_bit_block_SUB;
    bit_func_table_[bm::set_XOR] = 
        &serial_stream_iterator<DEC>::get_bit_block_XOR;
    bit_func_table_[bm::set_COUNT] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT;
    bit_func_table_[bm::set_COUNT_AND] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT_AND;
    bit_func_table_[bm::set_COUNT_XOR] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT_XOR;
    bit_func_table_[bm::set_COUNT_OR] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT_OR;
    bit_func_table_[bm::set_COUNT_SUB_AB] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT_SUB_AB;
    bit_func_table_[bm::set_COUNT_SUB_BA] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT_SUB_BA;
    bit_func_table_[bm::set_COUNT_A] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT_A;
    bit_func_table_[bm::set_COUNT_B] = 
        &serial_stream_iterator<DEC>::get_bit_block_COUNT;


    // reading stream header
    unsigned char header_flag =  decoder_.get_8();
    if (!(header_flag & BM_HM_NO_BO))
    {
        /*ByteOrder bo = (bm::ByteOrder)*/decoder_.get_8();
    }

    // check if bitvector comes as an inverted, sorted list of ints
    //
    if (header_flag & BM_HM_ID_LIST)
    {
        // special case: the next comes plain list of unsigned integers
        if (header_flag & BM_HM_RESIZE)
        {
            bv_size_ = decoder_.get_32();
        }

        state_ = e_list_ids;
        id_cnt_ = decoder_.get_32();
        next(); // read first id
    }
    else
    {
        if (!(header_flag & BM_HM_NO_GAPL)) 
        {
            unsigned i;
            // keep GAP levels info
            for (i = 0; i < bm::gap_levels; ++i)
            {
                glevels_[i] = decoder_.get_16();
            }
        }

        if (header_flag & (1 << 1))
        {
            bv_size_ = decoder_.get_32();
        }
        state_ = e_blocks;
    }
}

template<class DEC>
void serial_stream_iterator<DEC>::next()
{
    if (is_eof())
    {
        ++block_idx_;
        return;
    }

    switch (state_) 
    {
    case e_list_ids:
        // read inverted ids one by one
        if (id_cnt_ == 0)
        {
            end_of_stream_ = true;
            state_ = e_unknown;
            break;
        }
        last_id_ = decoder_.get_32();
        --id_cnt_;
        break;

    case e_blocks:
        if (block_idx_ == bm::set_total_blocks)
        {
            end_of_stream_ = true;
            state_ = e_unknown;
            break;
        }

        block_type_ = decoder_.get_8();

        // pre-check for 7-bit zero block
        //
        if (block_type_ & (1u << 7u))
        {
            mono_block_cnt_ = (block_type_ & ~(1u << 7u)) - 1;
            state_ = e_zero_blocks;
            break;
        }

        switch (block_type_)
        {
        case set_block_azero:
        case set_block_end:
            end_of_stream_ = true; state_ = e_unknown;
            break;
        case set_block_1zero:
            state_ = e_zero_blocks;
            mono_block_cnt_ = 0;
            break;
        case set_block_8zero:
            state_ = e_zero_blocks;
            mono_block_cnt_ = decoder_.get_8()-1;
            break;
        case set_block_16zero:
            state_ = e_zero_blocks;
            mono_block_cnt_ = decoder_.get_16()-1;
            break;
        case set_block_32zero:
            state_ = e_zero_blocks;
            mono_block_cnt_ = decoder_.get_32()-1;
            break;
        case set_block_aone:
            state_ = e_one_blocks;
            mono_block_cnt_ = bm::set_total_blocks - block_idx_;
            break;
        case set_block_1one:
            state_ = e_one_blocks;
            mono_block_cnt_ = 0;
            break;
        case set_block_8one:
            state_ = e_one_blocks;
            mono_block_cnt_ = decoder_.get_8()-1;
            break;
        case set_block_16one:
            state_ = e_one_blocks;
            mono_block_cnt_ = decoder_.get_16()-1;
            break;
        case set_block_32one:
            state_ = e_one_blocks;
            mono_block_cnt_ = decoder_.get_32()-1;
            break;

        case set_block_bit:
        case set_block_bit_interval:
        case set_block_bit_0runs:
        case set_block_arrbit:
            state_ = e_bit_block;
            break;

        case set_block_gap:
        case set_block_gap_egamma:
            gap_head_ = (gap_word_t)
                (sizeof(gap_word_t) == 2 ? 
                    decoder_.get_16() : decoder_.get_32());
        case set_block_arrgap:
        case set_block_arrgap_egamma:
        case set_block_arrgap_egamma_inv:
        case set_block_arrgap_inv:
		case set_block_bit_1bit:
            state_ = e_gap_block;
            break;        
        case set_block_gapbit:
            state_ = e_gap_block; //e_bit_block; // TODO: make a better decision here
            break;
        
        default:
            BM_ASSERT(0);
        }// switch

        break;

    case e_zero_blocks:
    case e_one_blocks:
        ++block_idx_;
        if (!mono_block_cnt_)
        {
            state_ = e_blocks; // get new block token
            break;
        }
        --mono_block_cnt_;
        break;

    case e_unknown:
    default:
        BM_ASSERT(0);
    } // switch
}

template<class DEC>
void serial_stream_iterator<DEC>::skip_mono_blocks()
{
	BM_ASSERT(state_ == e_zero_blocks || state_ == e_one_blocks);
    if (!mono_block_cnt_)
    {
		++block_idx_;
    }
	else
	{
		block_idx_ += mono_block_cnt_+1;
		mono_block_cnt_ = 0;
	}
    state_ = e_blocks;
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_ASSIGN(
                                            bm::word_t*  dst_block,
                                            bm::word_t*  /*tmp_block*/)
{
    BM_ASSERT(this->state_ == e_bit_block);
    unsigned count = 0;

    switch (this->block_type_)
    {
    case set_block_bit:
        decoder_.get_32(dst_block, bm::set_block_size);
        break;
    case set_block_bit_0runs: 
        {
        if (dst_block)
            bit_block_set(dst_block, 0);
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            if (run_type)
            {
				decoder_.get_32(dst_block ? dst_block + j : dst_block, run_length);
			}
			j += run_length;
        } // for
        }
        break;
    case set_block_bit_interval:
        {
            unsigned head_idx = decoder_.get_16();
            unsigned tail_idx = decoder_.get_16();
            if (dst_block) 
            {
                for (unsigned i = 0; i < head_idx; ++i)
                    dst_block[i] = 0;
                decoder_.get_32(dst_block + head_idx, 
                                tail_idx - head_idx + 1);
                for (unsigned j = tail_idx + 1; j < bm::set_block_size; ++j)
                    dst_block[j] = 0;
            }
            else
            {
                int pos = int(tail_idx - head_idx) + 1;
                pos *= 4u;
                decoder_.seek(pos);
            }
        }
        break;
    case set_block_arrbit:
    case set_block_bit_1bit:
        get_arr_bit(dst_block, true /*clear target*/);
        break;
    case set_block_gapbit:
        BM_ASSERT(0);
        break;
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_OR(bm::word_t*  dst_block,
                                              bm::word_t*  /*tmp_block*/)
{
    BM_ASSERT(this->state_ == e_bit_block);
    unsigned count = 0;
    switch (block_type_)
    {
    case set_block_bit:
        decoder_.get_32_OR(dst_block, bm::set_block_size);
        break;
    case set_block_bit_interval:
        {
        unsigned head_idx = decoder_.get_16();
        unsigned tail_idx = decoder_.get_16();
        for (unsigned i = head_idx; i <= tail_idx; ++i)
            dst_block[i] |= decoder_.get_32();
        }
        break;
    case set_block_bit_0runs:
        {
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            if (run_type)
            {
                unsigned run_end = j + run_length;
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    dst_block[j] |= decoder_.get_32();
                }
            }
            else
            {
                j += run_length;
            }
        } // for
        }
        break;
    case set_block_bit_1bit:
    case set_block_arrbit:
        get_arr_bit(dst_block, false /*don't clear target*/);
        break;		
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_AND(bm::word_t* BMRESTRICT dst_block,
                                               bm::word_t* BMRESTRICT tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block != tmp_block);
    unsigned count = 0;
    switch (block_type_)
    {
    case set_block_bit:
        decoder_.get_32_AND(dst_block, bm::set_block_size);
        break;
    case set_block_bit_0runs:
        {
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();

            unsigned run_end = j + run_length;
            if (run_type)
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    dst_block[j] &= decoder_.get_32();
                }
            }
            else
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    dst_block[j] = 0;
                }
            }
        } // for
        }
        break;
    case set_block_bit_interval:
        {
            unsigned head_idx = decoder_.get_16();
            unsigned tail_idx = decoder_.get_16();
            unsigned i;
            for ( i = 0; i < head_idx; ++i)
                dst_block[i] = 0;
            for ( i = head_idx; i <= tail_idx; ++i)
                dst_block[i] &= decoder_.get_32();
            for ( i = tail_idx + 1; i < bm::set_block_size; ++i)
                dst_block[i] = 0;
        }
        break;
    case set_block_bit_1bit:
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /*clear target*/);
        if (dst_block)
            bit_block_and(dst_block, tmp_block);
        break;		
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_XOR(bm::word_t*  dst_block,
                                               bm::word_t*  tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block != tmp_block);

    unsigned count = 0;
    switch (block_type_)
    {
    case set_block_bit:
        for (unsigned i = 0; i < bm::set_block_size; ++i)
            dst_block[i] ^= decoder_.get_32();
        break;
    case set_block_bit_0runs:
        {
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            if (run_type)
            {
                unsigned run_end = j + run_length;
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    dst_block[j] ^= decoder_.get_32();
                }
            }
            else
            {
                j += run_length;
            }
        } // for
        }
        break;
    case set_block_bit_interval:
        {
            unsigned head_idx = decoder_.get_16();
            unsigned tail_idx = decoder_.get_16();
            for (unsigned i = head_idx; i <= tail_idx; ++i)
                dst_block[i] ^= decoder_.get_32();
        }
        break;
    case set_block_bit_1bit:
        // TODO: optimization
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /*clear target*/);
        if (dst_block)
        {
            bit_block_xor(dst_block, tmp_block);
        }
        break;
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_SUB(bm::word_t*  dst_block,
                                               bm::word_t*  tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block != tmp_block);

    unsigned count = 0;
    switch (block_type_)
    {
    case set_block_bit:
        for (unsigned i = 0; i < bm::set_block_size; ++i)
            dst_block[i] &= ~decoder_.get_32();
        break;
    case set_block_bit_0runs:
        {
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            if (run_type)
            {
                unsigned run_end = j + run_length;
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    dst_block[j] &= ~decoder_.get_32();
                }
            }
            else
            {
                j += run_length;
            }
        } // for
        }
        break;
    case set_block_bit_interval:
        {
            unsigned head_idx = decoder_.get_16();
            unsigned tail_idx = decoder_.get_16();
            for (unsigned i = head_idx; i <= tail_idx; ++i)
                dst_block[i] &= ~decoder_.get_32();
        }
        break;
    case set_block_bit_1bit:
        // TODO: optimization
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /*clear target*/);
        if (dst_block)
            bit_block_sub(dst_block, tmp_block);
        break;
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}


template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_COUNT(bm::word_t*  /*dst_block*/,
                                                 bm::word_t*  /*tmp_block*/)
{
    BM_ASSERT(this->state_ == e_bit_block);

    unsigned count = 0;
    switch (block_type_)
    {
    case set_block_bit:
        for (unsigned i = 0; i < bm::set_block_size; ++i)
            count += word_bitcount(decoder_.get_32());
        break;
    case set_block_bit_0runs:
        {
        //count = 0;
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            if (run_type)
            {
                unsigned run_end = j + run_length;
                for (;j < run_end; ++j)
                {
                    count += word_bitcount(decoder_.get_32());
                }
            }
            else
            {
                j += run_length;
            }
        } // for
        return count;
        }
    case set_block_bit_interval:
        {
            unsigned head_idx = decoder_.get_16();
            unsigned tail_idx = decoder_.get_16();
            for (unsigned i = head_idx; i <= tail_idx; ++i)
                count += word_bitcount(decoder_.get_32());
        }
        break;
    case set_block_arrbit:
        count += get_arr_bit(0);
        break;
    case set_block_bit_1bit:
        ++count;
        decoder_.get_16();
        break;
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_COUNT_A(bm::word_t*  dst_block,
                                                   bm::word_t*  /*tmp_block*/)
{
    BM_ASSERT(this->state_ == e_bit_block);
    unsigned count = 0;
    if (dst_block)
    {
        // count the block bitcount
        count = bm::bit_block_count(dst_block);
    }

    switch (block_type_)
    {
    case set_block_bit:
        decoder_.get_32(0, bm::set_block_size);
        break;
    case set_block_bit_0runs:
        {
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            if (run_type)
            {
                unsigned run_end = j + run_length;
                for (;j < run_end; ++j)
                {
                    decoder_.get_32();
                }
            }
            else
            {
                j += run_length;
            }
        } // for
        }
        break;

    case set_block_bit_interval:
        {
            unsigned head_idx = decoder_.get_16();
            unsigned tail_idx = decoder_.get_16();
            for (unsigned i = head_idx; i <= tail_idx; ++i)
                decoder_.get_32();
        }
        break;
    case set_block_arrbit:
        get_arr_bit(0);
        break;
    case set_block_bit_1bit:
        decoder_.get_16();
        break;
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}


template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_COUNT_AND(bm::word_t*  dst_block,
                                                     bm::word_t*  tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block);

    unsigned count = 0;
    switch (block_type_)
    {
    case set_block_bit:
        for (unsigned i = 0; i < bm::set_block_size; ++i)
            count += word_bitcount(dst_block[i] & decoder_.get_32());
        break;
    case set_block_bit_0runs:
        {
        //count = 0;
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            if (run_type)
            {
                unsigned run_end = j + run_length;
                for (;j < run_end; ++j)
                {
                    count += word_bitcount(dst_block[j] & decoder_.get_32());
                }
            }
            else
            {
                j += run_length;
            }
        } // for
        return count;
        }
    case set_block_bit_interval:
        {
        unsigned head_idx = decoder_.get_16();
        unsigned tail_idx = decoder_.get_16();
        for (unsigned i = head_idx; i <= tail_idx; ++i)
            count += word_bitcount(dst_block[i] & decoder_.get_32());
        }
        break;
    case set_block_bit_1bit:
        // TODO: optimization
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /*clear target*/);
        count += 
            bit_operation_and_count(dst_block, tmp_block);
        break;
    default:
        BM_ASSERT(0);
    } // switch
    return count;
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_COUNT_OR(bm::word_t*  dst_block,
                                                    bm::word_t*  tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block);

    bitblock_sum_adapter count_adapter;
    switch (block_type_)
    {
    case set_block_bit:
        {
        bitblock_get_adapter ga(dst_block);
        bit_COUNT_OR<bm::word_t> func;
        
        bit_recomb(ga,
                   decoder_,
                   func,
                   count_adapter
                  );
        }
        break;
    case set_block_bit_0runs: 
        {
        unsigned count = 0;
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            unsigned run_end = j + run_length;
            if (run_type)
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    count += word_bitcount(dst_block[j] | decoder_.get_32());
                }
            }
            else
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    count += word_bitcount(dst_block[j]);
                }
            }
        } // for
        return count;
        }
    case set_block_bit_interval:
        {
        unsigned head_idx = decoder_.get_16();
        unsigned tail_idx = decoder_.get_16();
        unsigned count = 0;
        unsigned i;
        for (i = 0; i < head_idx; ++i)
            count += word_bitcount(dst_block[i]);
        for (i = head_idx; i <= tail_idx; ++i)
            count += word_bitcount(dst_block[i] | decoder_.get_32());
        for (i = tail_idx + 1; i < bm::set_block_size; ++i)
            count += word_bitcount(dst_block[i]);
        return count;
        }
    case set_block_bit_1bit:
        // TODO: optimization
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /* clear target*/);
        return 
            bit_operation_or_count(dst_block, tmp_block);
    default:
        BM_ASSERT(0);
    } // switch
    return count_adapter.sum();
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_COUNT_XOR(bm::word_t*  dst_block,
                                                     bm::word_t*  tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block);

    bitblock_sum_adapter count_adapter;
    switch (block_type_)
    {
    case set_block_bit:
        {
        bitblock_get_adapter ga(dst_block);
        bit_COUNT_XOR<bm::word_t> func;
        
        bit_recomb(ga,
                   decoder_,
                   func,
                   count_adapter
                  );
        }
        break;
    case set_block_bit_0runs: 
        {
        unsigned count = 0;
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            unsigned run_end = j + run_length;
            if (run_type)
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    count += word_bitcount(dst_block[j] ^ decoder_.get_32());
                }
            }
            else
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    count += word_bitcount(dst_block[j]);
                }
            }
        } // for
        return count;
        }
    case set_block_bit_interval:
        {
        unsigned head_idx = decoder_.get_16();
        unsigned tail_idx = decoder_.get_16();
        unsigned count = 0;
        unsigned i;
        for (i = 0; i < head_idx; ++i)
            count += word_bitcount(dst_block[i]);
        for (i = head_idx; i <= tail_idx; ++i)
            count += word_bitcount(dst_block[i] ^ decoder_.get_32());
        for (i = tail_idx + 1; i < bm::set_block_size; ++i)
            count += word_bitcount(dst_block[i]);
        return count;
        }
    case set_block_bit_1bit:
        // TODO: optimization
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /* clear target*/);
        return 
            bit_operation_xor_count(dst_block, tmp_block);
    default:
        BM_ASSERT(0);
    } // switch
    return count_adapter.sum();
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_COUNT_SUB_AB(bm::word_t*  dst_block,
                                                        bm::word_t*  tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block);

    bitblock_sum_adapter count_adapter;
    switch (block_type_)
    {
    case set_block_bit:
        {
        bitblock_get_adapter ga(dst_block);
        bit_COUNT_SUB_AB<bm::word_t> func;
        
        bit_recomb(ga, 
                   decoder_,
                   func,
                   count_adapter
                  );
        }
        break;
    case set_block_bit_0runs: 
        {
        unsigned count = 0;
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            unsigned run_end = j + run_length;
            if (run_type)
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    count += word_bitcount(dst_block[j] & ~decoder_.get_32());
                }
            }
            else
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    count += word_bitcount(dst_block[j]);
                }
            }
        } // for
        return count;
        }
    case set_block_bit_interval:
        {
        unsigned head_idx = decoder_.get_16();
        unsigned tail_idx = decoder_.get_16();
        unsigned count = 0;
        unsigned i;
        for (i = 0; i < head_idx; ++i)
            count += word_bitcount(dst_block[i]);
        for (i = head_idx; i <= tail_idx; ++i)
            count += word_bitcount(dst_block[i] & (~decoder_.get_32()));
        for (i = tail_idx + 1; i < bm::set_block_size; ++i)
            count += word_bitcount(dst_block[i]);
        return count;
        }
        break;
    case set_block_bit_1bit:
        //TODO: optimization
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /* clear target*/);
        return 
            bit_operation_sub_count(dst_block, tmp_block);
    default:
        BM_ASSERT(0);
    } // switch
    return count_adapter.sum();
}

template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block_COUNT_SUB_BA(bm::word_t*  dst_block,
                                                        bm::word_t*  tmp_block)
{
    BM_ASSERT(this->state_ == e_bit_block);
    BM_ASSERT(dst_block);

    bitblock_sum_adapter count_adapter;
    switch (block_type_)
    {
    case set_block_bit:
        {
        bitblock_get_adapter ga(dst_block);
        bit_COUNT_SUB_BA<bm::word_t> func;

        bit_recomb(ga,
                   decoder_,
                   func,
                   count_adapter
                  );
        }
        break;
    case set_block_bit_0runs: 
        {
        unsigned count = 0;
        unsigned char run_type = decoder_.get_8();
        for (unsigned j = 0; j < bm::set_block_size;run_type = !run_type)
        {
            unsigned run_length = decoder_.get_16();
            unsigned run_end = j + run_length;
            if (run_type)
            {
                for (;j < run_end; ++j)
                {
                    BM_ASSERT(j < bm::set_block_size);
                    count += word_bitcount(decoder_.get_32() & (~dst_block[j]));
                }
            }
            else
            {
                j += run_length;
            }
        } // for
        return count;
        }

    case set_block_bit_interval:
        {
        unsigned head_idx = decoder_.get_16();
        unsigned tail_idx = decoder_.get_16();
        unsigned count = 0;
        unsigned i;
        for (i = head_idx; i <= tail_idx; ++i)
            count += word_bitcount(decoder_.get_32() & (~dst_block[i]));
        return count;
        }
        break;
    case set_block_bit_1bit:
        // TODO: optimization
    case set_block_arrbit:
        get_arr_bit(tmp_block, true /* clear target*/);
        return 
            bit_operation_sub_count(tmp_block, dst_block);
    default:
        BM_ASSERT(0);
    } // switch
    return count_adapter.sum();
}



template<class DEC>
unsigned serial_stream_iterator<DEC>::get_arr_bit(bm::word_t* dst_block, 
                                                  bool        clear_target)
{
    BM_ASSERT(this->block_type_ == set_block_arrbit || 
              this->block_type_ == set_block_bit_1bit);
    
    gap_word_t len = decoder_.get_16(); // array length / 1bit_idx
    if (dst_block)
    {
        if (clear_target)
            bit_block_set(dst_block, 0);

        if (this->block_type_ == set_block_bit_1bit)
        {
            // len contains idx of 1 bit set
            set_bit(dst_block, len);
            return 1;
        }

        for (unsigned k = 0; k < len; ++k)
        {
            gap_word_t bit_idx = decoder_.get_16();
            set_bit(dst_block, bit_idx);
        }
    }
    else
    {
        if (this->block_type_ == set_block_bit_1bit)
        {
            return 1; // nothing to do: len var already consumed 16bits
        }
        // fwd the decocing stream
        decoder_.seek(len * 2);
    }
    return len;
}

template<class DEC>
unsigned serial_stream_iterator<DEC>::get_bit()
{
    BM_ASSERT(this->block_type_ == set_block_bit_1bit);
    ++(this->block_idx_);
    this->state_ = e_blocks;

	return decoder_.get_16(); // 1bit_idx	
}

template<class DEC>
void 
serial_stream_iterator<DEC>::get_gap_block(bm::gap_word_t* dst_block)
{
    BM_ASSERT(this->state_ == e_gap_block || 
              this->block_type_ == set_block_bit_1bit);
    BM_ASSERT(dst_block);

    this->read_gap_block(this->decoder_,
                   this->block_type_,
                   dst_block,
                   this->gap_head_);

    ++(this->block_idx_);
    this->state_ = e_blocks;
}


template<class DEC>
unsigned 
serial_stream_iterator<DEC>::get_bit_block(bm::word_t*    dst_block,
                                           bm::word_t*    tmp_block,
                                           set_operation  op)
{
    BM_ASSERT(this->state_ == e_bit_block);
    get_bit_func_type bit_func = bit_func_table_[op];
    BM_ASSERT(bit_func);
    unsigned cnt = ((*this).*(bit_func))(dst_block, tmp_block);
    this->state_ = e_blocks;
    ++(this->block_idx_);
    return cnt;
}



template<class BV>
unsigned operation_deserializer<BV>::deserialize(
                                        bvector_type&        bv, 
                                        const unsigned char* buf, 
                                        bm::word_t*          temp_block,
                                        set_operation        op,
                                        bool                 exit_on_one
                                        )
{
    ByteOrder bo_current = globals<true>::byte_order();

    bm::decoder dec(buf);
    unsigned char header_flag = dec.get_8();
    ByteOrder bo = bo_current;
    if (!(header_flag & BM_HM_NO_BO))
    {
        bo = (bm::ByteOrder) dec.get_8();
    }

    blocks_manager_type& bman = bv.get_blocks_manager();
    bit_block_guard<blocks_manager_type> bg(bman);
    if (temp_block == 0)
    {
        temp_block = bg.allocate();
    }

    if (bo_current == bo)
    {
        serial_stream_current ss(buf);
        return 
            iterator_deserializer<BV, serial_stream_current>::
                deserialize(bv, ss, temp_block, op, exit_on_one);
    }
    switch (bo_current) 
    {
    case BigEndian:
        {
        serial_stream_be ss(buf);
        return 
            iterator_deserializer<BV, serial_stream_be>::
                deserialize(bv, ss, temp_block, op, exit_on_one);
        }
    case LittleEndian:
        {
        serial_stream_le ss(buf);
        return 
            iterator_deserializer<BV, serial_stream_le>::
                deserialize(bv, ss, temp_block, op, exit_on_one);
        }
    default:
        BM_ASSERT(0);
    };
    return 0;
}


template<class BV>
void operation_deserializer<BV>::deserialize(
                     bvector_type&        bv_target,
                     const bvector_type&  bv_mask,
                     const unsigned char* buf, 
                     bm::word_t*          temp_block,
                     set_operation        op)
{
    ByteOrder bo_current = globals<true>::byte_order();

    bm::decoder dec(buf);
    unsigned char header_flag = dec.get_8();
    ByteOrder bo = bo_current;
    if (!(header_flag & BM_HM_NO_BO))
    {
        bo = (bm::ByteOrder) dec.get_8();
    }

    blocks_manager_type& bman = bv_target.get_blocks_manager();
    if (!bman.is_init())
    {
        bman.init_tree();
    }
    bit_block_guard<blocks_manager_type> bg(bman);
    if (temp_block == 0)
    {
        temp_block = bg.allocate();
    }

    if (bo_current == bo)
    {
        serial_stream_current ss(buf);
        iterator_deserializer<BV, serial_stream_current>::
            deserialize(bv_target, bv_mask, ss, temp_block, op);
        return;
    }
    switch (bo_current) 
    {
    case BigEndian:
        {
            serial_stream_be ss(buf);
            iterator_deserializer<BV, serial_stream_be>::
                deserialize(bv_target, bv_mask, ss, temp_block, op);
        }
    case LittleEndian:
        {
            serial_stream_le ss(buf);
            iterator_deserializer<BV, serial_stream_le>::
                deserialize(bv_target, bv_mask, ss, temp_block, op);
        }
    default:
        BM_ASSERT(0);
    };
}


template<class BV, class SerialIterator>
void iterator_deserializer<BV, SerialIterator>::load_id_list(
                                            bvector_type&         bv, 
                                            serial_iterator_type& sit,
                                            unsigned              id_count,
                                            bool                  set_clear)
{
    const unsigned win_size = 64;
    bm::id_t id_buffer[win_size+1];

    if (set_clear)  // set bits
    {
        for (unsigned i = 0; i <= id_count;)
        {
            unsigned j;
            for (j = 0; j < win_size && i <= id_count; ++j, ++i) 
            {
                id_buffer[j] = sit.get_id();
                sit.next();
            } // for j
            bm::combine_or(bv, id_buffer, id_buffer + j);
        } // for i
    } 
    else // clear bits
    {
        for (unsigned i = 0; i <= id_count;)
        {
            unsigned j;
            for (j = 0; j < win_size && i <= id_count; ++j, ++i) 
            {
                id_buffer[j] = sit.get_id();
                sit.next();
            } // for j
            bm::combine_sub(bv, id_buffer, id_buffer + j);
        } // for i
    }
}

template<class BV, class SerialIterator>
unsigned 
iterator_deserializer<BV, SerialIterator>::finalize_target_vector(
                                                blocks_manager_type& bman,
                                                set_operation        op,
                                                unsigned             bv_block_idx)
{
    unsigned count = 0;
    switch (op)
    {
    case set_OR:    case set_SUB:     case set_XOR:
    case set_COUNT: case set_COUNT_B: case set_COUNT_AND:
    case set_COUNT_SUB_BA:
        // nothing to do
        break;
    case set_AND: case set_ASSIGN:
        // clear the rest of the target vector
        {
            unsigned i, j;
            bman.get_block_coord(bv_block_idx, i, j);
            bm::word_t*** blk_root = bman.top_blocks_root();
            unsigned top_size = bman.top_block_size();
            for (; i < top_size; ++i)
            {
                bm::word_t** blk_blk = blk_root[i];
                if (blk_blk == 0) 
                {
                    bv_block_idx+=bm::set_array_size-j;
                    j = 0;
                    continue;
                }
                for (;j < bm::set_array_size; ++j, ++bv_block_idx)
                {
                    //if (blk_blk[j])
                        bman.zero_block(bv_block_idx);
                } // for j
                j = 0;
            } // for i

        }
        break;
    case set_COUNT_A: case set_COUNT_OR: case set_COUNT_XOR:
    case set_COUNT_SUB_AB:
        // count bits in the target vector
        {
            unsigned i, j;
            bman.get_block_coord(bv_block_idx, i, j);
            bm::word_t*** blk_root = bman.top_blocks_root();
            unsigned top_size = bman.top_block_size();
            for (;i < top_size; ++i)
            {
                bm::word_t** blk_blk = blk_root[i];
                if (blk_blk == 0) 
                {
                    bv_block_idx+=bm::set_array_size-j;
                    j = 0;
                    continue;
                }
                for (;j < bm::set_array_size; ++j, ++bv_block_idx)
                {
                    if (blk_blk[j])
                        count += bman.block_bitcount(blk_blk[j]);
                } // for j
                j = 0;
            } // for i

        }
        break;
    case set_END:
    default:
        BM_ASSERT(0);
    }
    return count;
}

template<class BV, class SerialIterator>
unsigned 
iterator_deserializer<BV, SerialIterator>::process_id_list(
                                    bvector_type&         bv, 
                                    serial_iterator_type& sit,
                                    set_operation         op)
{
    unsigned count = 0;
    unsigned id_count = sit.get_id_count();
    bool set_clear = true;
    switch (op)
    {
    case set_AND:
        {
            // TODO: use some more optimal algorithm without temp vector
            BV bv_tmp(BM_GAP);
            load_id_list(bv_tmp, sit, id_count, true);
            bv &= bv_tmp;
        }
        break;
    case set_ASSIGN:
        bv.clear(true);
        // intentional case fall through here (not a bug)
    case set_OR:
        set_clear = true;
        // fall to SUB
    case set_SUB:
        load_id_list(bv, sit, id_count, set_clear);
        break;
    case set_XOR:
        for (unsigned i = 0; i < id_count; ++i)
        {
            bm::id_t id = sit.get_id();
            bv[id] ^= true;
            sit.next();
        } // for
        break;
    case set_COUNT: case set_COUNT_B:
        for (unsigned i = 0; i < id_count; ++i)
        {
            /* bm::id_t id = */ sit.get_id();
            ++count;
            sit.next();
        } // for
        break;
    case set_COUNT_A:
        return bv.count();
    case set_COUNT_AND:
        for (unsigned i = 0; i < id_count; ++i)
        {
            bm::id_t id = sit.get_id();
            count += bv.get_bit(id);
            sit.next();
        } // for
        break;
    case set_COUNT_XOR:
        {
            // TODO: get rid of the temp vector
            BV bv_tmp(BM_GAP);
            load_id_list(bv_tmp, sit, id_count, true);
            count += bm::count_xor(bv, bv_tmp);
        }
        break;
    case set_COUNT_OR:
        {
            // TODO: get rid of the temp. vector
            BV bv_tmp(BM_GAP);
            load_id_list(bv_tmp, sit, id_count, true);
            count += bm::count_or(bv, bv_tmp);
        }
        break;
    case set_COUNT_SUB_AB:
        {
            // TODO: get rid of the temp. vector
            BV bv_tmp(bv);
            load_id_list(bv_tmp, sit, id_count, false);
            count += bv_tmp.count();
        }
        break;
    case set_COUNT_SUB_BA:
        {
            BV bv_tmp(BM_GAP);
            load_id_list(bv_tmp, sit, id_count, true);
            count += bm::count_sub(bv_tmp, bv);        
        }
        break;
    case set_END:
    default:
        BM_ASSERT(0);
    } // switch

    return count;
}


template<class BV, class SerialIterator>
void iterator_deserializer<BV, SerialIterator>::deserialize(
                     bvector_type&         bv_target,
                     const bvector_type&   bv_mask,
                     serial_iterator_type& sit, 
                     bm::word_t*           temp_block,
                     set_operation         op)
{
    BM_ASSERT(temp_block);
    BM_ASSERT(op == bm::set_AND ||
              op == bm::set_OR || op == bm::set_XOR || op == bm::set_SUB);

    gap_word_t   gap_temp_block[bm::gap_equiv_len*3];
    gap_temp_block[0] = 0;

    bv_target.clear(true); // clear and free the memory

    const blocks_manager_type& bman_mask = bv_mask.get_blocks_manager();
          blocks_manager_type& bman_target = bv_target.get_blocks_manager();
    
    if (!bman_target.is_init())
    {
        bman_target.init_tree();
    }

    unsigned bv_size = sit.bv_size();
    if (bv_mask.size() > bv_size) 
    {
        bv_size = bv_mask.size();    
    }
    if (bv_target.size() < bv_size)
    {
        bv_target.resize(bv_size);
    }

    unsigned top_blocks = bman_mask.effective_top_block_size();

    BM_SET_MMX_GUARD

    typename serial_iterator_type::iterator_state state;
    state = sit.get_state();
    if (state == serial_iterator_type::e_list_ids)
    {
        bv_target = bv_mask;
        process_id_list(bv_target, sit, op);
        return;
    }

    unsigned bv_block_idx = 0;
    for (;1;)
    {
		bv_block_idx = sit.block_idx();
        // early exit check to avoid over-scan
        {
            unsigned tb_idx = bv_block_idx >> bm::set_array_shift; // current top block
            if (tb_idx > top_blocks)
            {
                if (op == bm::set_AND)
                    break;
                if (sit.is_eof())
                    break;
            }
        }
        
        if (sit.is_eof()) // argument stream ended
        {
            if (op == bm::set_AND)
                break;
            // (OR, SUB, XOR) need to scan fwd until mask vector ends
            state = serial_iterator_type::e_zero_blocks;
        }
        else 
        {
            state = sit.state();
        }

        switch (state)
        {
        case serial_iterator_type::e_blocks:
            sit.next();
            continue;
        case serial_iterator_type::e_bit_block:
            {
		        bm::set_operation sop = op;
                const bm::word_t* blk_mask = bman_mask.get_block(bv_block_idx);
                bm::word_t* blk_target = 0;
                if (!blk_mask) 
                {
                    switch (op)
                    {
                    case set_AND: case set_SUB: 
                        // first arg is 0, so the operation gives us zero
                        // all we need to do is to seek the input stream
                        sop = set_ASSIGN;
                        break;
                    case set_OR: case set_XOR:
                        blk_target = bman_target.make_bit_block(bv_block_idx);
                        break;
                    case set_ASSIGN:
                    case set_COUNT:
                    case set_COUNT_AND:
                    case set_COUNT_XOR:
                    case set_COUNT_OR:
                    case set_COUNT_SUB_AB:
                    case set_COUNT_SUB_BA:
                    case set_COUNT_A:
                    case set_COUNT_B:
                    case set_END:
                    default:
                        BM_ASSERT(0);
                    }
                }
                else // block exists
                {
                    int is_gap = BM_IS_GAP(blk_mask);
                    blk_target = bman_target.copy_bit_block(bv_block_idx, blk_mask, is_gap);
                }

                // 2 bit blocks recombination
                sit.get_bit_block(blk_target, temp_block, sop);
            }
            break;

        case serial_iterator_type::e_zero_blocks:
            {
				if (op == set_AND)
				{
					sit.skip_mono_blocks();
					break;
				}
                sit.next();
   			    // set_SUB: set_OR: set_XOR: 
				bman_target.copy_block(bv_block_idx, bman_mask);
            }
            break;

        case serial_iterator_type::e_one_blocks:
            {
                BM_ASSERT(bv_block_idx == sit.block_idx());
                const bm::word_t* blk_mask = bman_mask.get_block(bv_block_idx);
                sit.next();

                switch (op)
                {
                case set_OR:
                    bman_target.set_block_all_set(bv_block_idx);
                    break;
                case set_SUB:
                    break;
                case set_AND:
					bman_target.copy_block(bv_block_idx, bman_mask);
                    break;
                case set_XOR:
                    if (blk_mask)
                    {
                        int is_gap = BM_IS_GAP(blk_mask);
                        bm::word_t* blk_target = 
                            bman_target.copy_bit_block(bv_block_idx, blk_mask, is_gap);
                        bit_block_xor(blk_target, FULL_BLOCK_REAL_ADDR);
                    }
                    else
                    {
                        // 0 XOR 1 = 1
                        bman_target.set_block_all_set(bv_block_idx);
                    }
                    break;
                default:
                    BM_ASSERT(0);
                } // switch

            }
            break;

        case serial_iterator_type::e_gap_block:
            {
				// Single bit-in-block optimization				
				if (sit.get_block_type() == set_block_bit_1bit)
				{
					if (op == set_AND)
					{
						unsigned bit_idx = sit.get_bit();
						unsigned bn = (bv_block_idx << bm::set_block_shift) | bit_idx;
						bool bval_mask = bv_mask.test(bn);
						bv_target.set_bit(bn, bval_mask);						
						break;
					}
				}
				
                const bm::word_t* blk_mask = bman_mask.get_block(bv_block_idx);

                sit.get_gap_block(gap_temp_block);
                // gap_word_t   gap_head = gap_temp_block[0];

                unsigned len = gap_length(gap_temp_block);
                int level = gap_calc_level(len, bman_target.glen());
                --len;

                
                if (!blk_mask)
                {
					if (op == set_OR || op == set_XOR)
					{
                        bman_target.set_gap_block(bv_block_idx, gap_temp_block, level);
					}
                }
                else  // mask block exists
                {
                    bm::operation bop = bm::setop2op(op);
                    bman_target.copy_block(bv_block_idx, bman_mask);
                    bv_target.combine_operation_with_block(
                                        bv_block_idx, 
                                        (bm::word_t*)gap_temp_block, 
                                        1,  // GAP
                                        bop);
                }
                
            }
            break;

        default:
            BM_ASSERT(0);
        } // switch


    } // for (deserialization)

	bv_target.forget_count();
}



template<class BV, class SerialIterator>
unsigned 
iterator_deserializer<BV, SerialIterator>::deserialize(
                                       bvector_type&         bv, 
                                       serial_iterator_type& sit, 
                                       bm::word_t*           temp_block,
                                       set_operation         op,
                                       bool                  exit_on_one)
{
    BM_ASSERT(temp_block);

    unsigned count = 0;
    gap_word_t   gap_temp_block[bm::gap_equiv_len*3];
    gap_temp_block[0] = 0;

    blocks_manager_type& bman = bv.get_blocks_manager();
    if (!bman.is_init())
    {
        bman.init_tree();
    }

    bv.forget_count();
    if (sit.bv_size() && (sit.bv_size() > bv.size())) 
    {
        bv.resize(sit.bv_size());
    }

    BM_SET_MMX_GUARD

    typename serial_iterator_type::iterator_state state;
    state = sit.get_state();
    if (state == serial_iterator_type::e_list_ids)
    {
        count = process_id_list(bv, sit, op);
        return count;
    }

    unsigned bv_block_idx = 0;

    for (;1;)
    {
        bm::set_operation sop = op;
        if (sit.is_eof()) // argument stream ended
        {
            count += finalize_target_vector(bman, op, bv_block_idx);
            return count;
        }

        state = sit.state();
        switch (state)
        {
        case serial_iterator_type::e_blocks:
            sit.next();
            continue;
        case serial_iterator_type::e_bit_block:
            {
            BM_ASSERT(sit.block_idx() == bv_block_idx);
            bm::word_t* blk = bman.get_block(bv_block_idx);

            if (!blk) 
            {
                switch (op)
                {
                case set_AND:          case set_SUB: case set_COUNT_AND:
                case set_COUNT_SUB_AB: case set_COUNT_A:
                    // one arg is 0, so the operation gives us zero
                    // all we need to do is to seek the input stream
                    sop = set_ASSIGN;
                    break;

                case set_OR: case set_XOR: case set_ASSIGN:
                    blk = bman.make_bit_block(bv_block_idx);
                    break;

                case set_COUNT:        case set_COUNT_XOR: case set_COUNT_OR:
                case set_COUNT_SUB_BA: case set_COUNT_B:
                    // first arg is not required (should work as is)
                    sop = set_COUNT;
                    break;

                case set_END:
                default:
                    BM_ASSERT(0);
                }
            }
            else // block exists
            {
                int is_gap = BM_IS_GAP(blk);
                if (is_gap || IS_FULL_BLOCK(blk))
                {
                    if (IS_FULL_BLOCK(blk) && is_const_set_operation(op))
                    {
                    }
                    else
                    {
                        // TODO: make sure const operations do not 
                        // deoptimize GAP blocks
                        blk = bman.deoptimize_block(bv_block_idx);
                    }
                }
            }

            // 2 bit blocks recombination
            unsigned c = sit.get_bit_block(blk, temp_block, sop);
            count += c;
			if (exit_on_one && count) // early exit
				return count;

            }
            break;

        case serial_iterator_type::e_zero_blocks:
            {
            BM_ASSERT(bv_block_idx == sit.block_idx());
            bm::word_t* blk = bman.get_block(bv_block_idx);
            sit.next();

            if (blk)
            {
                switch (op)
                {
                case set_AND: case set_ASSIGN:
                    // the result is 0
                    //blk =
                    bman.zero_block(bv_block_idx);
                    break;

                case set_SUB: case set_COUNT_AND:    case set_OR:
                case set_XOR: case set_COUNT_SUB_BA: case set_COUNT_B:
                    // nothing to do
                    break;
                
                case set_COUNT_SUB_AB: case set_COUNT_A: case set_COUNT_OR:
                case set_COUNT:        case set_COUNT_XOR:
                    // valid bit block recombined with 0 block
                    // results with same block data
                    // all we need is to bitcount bv block
                    {
                    count += blk ? bman.block_bitcount(blk) : 0;
					if (exit_on_one && count) // early exit
						return count;
                    }
                    break;
                case set_END:
                default:
                    BM_ASSERT(0);
                } // switch op
            } // if blk
            }
            break;

        case serial_iterator_type::e_one_blocks:
            {
            BM_ASSERT(bv_block_idx == sit.block_idx());
            bm::word_t* blk = bman.get_block(bv_block_idx);
            sit.next();

            switch (op)
            {
            case set_OR: case set_ASSIGN:
                bman.set_block_all_set(bv_block_idx);
                break;
            case set_COUNT_OR: case set_COUNT_B: case set_COUNT:
                // block becomes all set
                count += bm::bits_in_block;
                break;
            case set_SUB:
                //blk =
                bman.zero_block(bv_block_idx);
                break;
            case set_COUNT_SUB_AB: case set_AND:
                // nothing to do
                break;
            case set_COUNT_AND: case set_COUNT_A:
                count += blk ? bman.block_bitcount(blk) : 0;
                break;
            default:
                if (blk)
                {
                    switch (op)
                    {
                    case set_XOR:
                        blk = bman.deoptimize_block(bv_block_idx);
                        bit_block_xor(blk, FULL_BLOCK_REAL_ADDR);
                        break;
                    case set_COUNT_XOR:
                        {
                        count += 
                            combine_count_operation_with_block(
                                                blk,
                                                FULL_BLOCK_REAL_ADDR,
                                                bm::COUNT_XOR);
                        }
                        break;
                    case set_COUNT_SUB_BA:
                        {
                        count += 
                            combine_count_operation_with_block(
                                                blk,
                                                FULL_BLOCK_REAL_ADDR,
                                                bm::COUNT_SUB_BA);
                        }
                        break;
                    default:
                        BM_ASSERT(0);
                    } // switch
                }
                else // blk == 0 
                {
                    switch (op)
                    {
                    case set_XOR:
                        // 0 XOR 1 = 1
                        bman.set_block_all_set(bv_block_idx);
                        break;
                    case set_COUNT_XOR:
                        count += bm::bits_in_block;
                        break;
                    case set_COUNT_SUB_BA:
                        // 1 - 0 = 1
                        count += bm::bits_in_block;
                        break;
                    default:
                        break;
                    } // switch
                } // else
            } // switch
            if (exit_on_one && count) // early exit
				   return count;

            }
            break;

        case serial_iterator_type::e_gap_block:
            {
            BM_ASSERT(bv_block_idx == sit.block_idx());
            bm::word_t* blk = bman.get_block(bv_block_idx);

            sit.get_gap_block(gap_temp_block);

            unsigned len = gap_length(gap_temp_block);
            int level = gap_calc_level(len, bman.glen());
            --len;

            bool const_op = is_const_set_operation(op);
            if (const_op)
            {
                distance_metric metric = operation2metric(op);
                bm::word_t* gptr = (bm::word_t*)gap_temp_block;
                BMSET_PTRGAP(gptr);
                unsigned c = 
                    combine_count_operation_with_block(
                                        blk,
                                        gptr,
                                        metric);
                count += c;
                if (exit_on_one && count) // early exit
				    return count;

            }
            else // non-const set operation
            {
                if ((sop == set_ASSIGN) && blk) // target block override
                {
                    //blk =
                    bman.zero_block(bv_block_idx);
                    sop = set_OR;
                }
                if (blk == 0) // target block not found
                {
                    switch (sop)
                    {
                    case set_AND: case set_SUB:
                        break;
                    case set_OR: case set_XOR: case set_ASSIGN:
                        bman.set_gap_block(
                            bv_block_idx, gap_temp_block, level);
                        break;
                    
                    default:
                        BM_ASSERT(0);
                    } // switch
                }
                else  // target block exists
                {
                    bm::operation bop = bm::setop2op(op);
                    if (level == -1) // too big to GAP
                    {
                        gap_convert_to_bitset(temp_block, gap_temp_block);
                        bv.combine_operation_with_block(bv_block_idx, 
                                                        temp_block, 
                                                        0, // BIT
                                                        bop);
                    }
                    else // GAP fits
                    {
                        set_gap_level(gap_temp_block, level);
                        bv.combine_operation_with_block(
                                                bv_block_idx, 
                                                (bm::word_t*)gap_temp_block, 
                                                1,  // GAP
                                                bop);
                    }
                }
                if (exit_on_one) 
                {
                    blk = bman.get_block(bv_block_idx);
                    if (blk) 
                    {
                        bool z = bm::check_block_zero(blk, true/*deep check*/);
                        if (!z) 
                            return 1;
                    } 
                } // if exit_on_one

            } // if else non-const op
            }
            break;

        default:
            BM_ASSERT(0);
        } // switch

        ++bv_block_idx;

    } // for (deserialization)

    return count;
}




} // namespace bm

#include "bmundef.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif


#endif
