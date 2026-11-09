// Freeze released BM v9.3.1 opcode IDs; validate the new shared namespace.
#include "bm.h"
#include "bmserial.h"
#include <cassert>
#include <iostream>

constexpr unsigned char released_opcodes[] = {
    bm::set_block_end, bm::set_block_1zero, bm::set_block_1one,
    bm::set_block_8zero, bm::set_block_8one, bm::set_block_16zero,
    bm::set_block_16one, bm::set_block_32zero, bm::set_block_32one,
    bm::set_block_azero, bm::set_block_aone, bm::set_block_bit,
    bm::set_block_sgapbit, bm::set_block_sgapgap, bm::set_block_gap,
    bm::set_block_gapbit, bm::set_block_arrbit, bm::set_block_bit_interval,
    bm::set_block_arrgap, bm::set_block_bit_1bit, bm::set_block_gap_egamma,
    bm::set_block_arrgap_egamma, bm::set_block_bit_0runs, bm::set_block_arrgap_egamma_inv,
    bm::set_block_arrgap_inv, bm::set_block_64zero, bm::set_block_64one,
    bm::set_block_gap_bienc, bm::set_block_arrgap_bienc, bm::set_block_arrgap_bienc_inv,
    bm::set_block_arrbit_inv, bm::set_block_arr_bienc, bm::set_block_arr_bienc_inv,
    bm::set_block_bitgap_bienc, bm::set_block_bit_digest0, bm::set_block_ref_eq,
    bm::set_block_xor_ref8, bm::set_block_xor_ref16, bm::set_block_xor_ref32,
    bm::set_block_xor_gap_ref8, bm::set_block_xor_gap_ref16, bm::set_block_xor_gap_ref32,
    bm::set_block_xor_chain, bm::set_block_gap_bienc_v2, bm::set_block_arrgap_bienc_v2,
    bm::set_block_arrgap_bienc_inv_v2, bm::set_block_bitgap_bienc_v2, bm::set_nb_bookmark16,
    bm::set_nb_bookmark24, bm::set_nb_bookmark32, bm::set_nb_sync_mark8,
    bm::set_nb_sync_mark16, bm::set_nb_sync_mark24, bm::set_nb_sync_mark32,
    bm::set_nb_sync_mark48, bm::set_nb_sync_mark64, bm::set_sblock_bienc,
    bm::set_block_arr_bienc_8bh, bm::set_block_xor_ref8_um, bm::set_block_xor_ref16_um,
    bm::set_block_xor_ref32_um, bm::set_block_gap_bienc_v3, bm::set_block_gap_bienc_v3s,
    bm::set_block_arr_bienc_v3, bm::set_block_arr_bienc_inv_v3, bm::set_block_arr_bienc_v3s,
    bm::set_block_arr_bienc_inv_v3s, bm::set_block_gap_egamma_v3, bm::set_sblock_bienc_v3,
    bm::set_sblock_bienc_gaps_v3,
};
constexpr bool released_assignments_unchanged()
{
    for (unsigned i = 0; i < sizeof(released_opcodes); ++i)
        if (released_opcodes[i] != i) return false;
    return true;
}
static_assert(sizeof(released_opcodes) == 70 && released_assignments_unchanged(),
              "Released opcode IDs 0..69 must remain unchanged");
static_assert(bm::serialization_new_opcode_numbers_valid(), "Invalid new opcode registry");
int main()
{
    bool used[128] = {};
    for (unsigned code : released_opcodes) used[code] = true;
    for (unsigned code : bm::serialization_new_opcodes) { assert(!used[code]); used[code] = true; }
    for (unsigned code : {70u,71u,82u,83u,86u,87u,88u,89u,91u,93u,95u}) assert(!used[code]);
    for (unsigned code = 99; code < 128; ++code) assert(!used[code]);
    unsigned vacant = 0; for (bool occupied : used) vacant += !occupied;
    assert(vacant == 40);
    std::cout << "Opcode registry: 70 released IDs unchanged, 18 retained IDs unique, 40 vacant PASS\n";
}
