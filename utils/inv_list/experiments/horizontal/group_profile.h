// Measurement-only support for the Gov2 isolated-group survey.
// Borrow completed level-10 BLOBs to attribute both alternatives by opcode.
#ifndef INV_LIST_GROUP_PROFILE_H
#define INV_LIST_GROUP_PROFILE_H

#include <iostream>
#include "bmdbg.h"

namespace
{
/// Write code, name, count and record bytes for each emitted opcode in profile.
void write_records(const bm::serialization_opcode_profile& profile)
{
    bool first = true;
    for (unsigned i = 0; i < 256; ++i)
    {
        const auto& record = profile.opcode[i];
        if (record.count)
        {
            std::cout << (first ? "" : ",") << i << ':' << bm::serialization_opcode_name(i)
                      << ':' << record.count << ':' << record.bytes;
            first = false;
        }
    } // for i
}

/// Measure group at its original ordinal and superblock index, using both
/// serializers and reusable scratch. ids contains n sorted group-local IDs.
template<class BV>
void profile_group(const BV& group, unsigned ordinal, unsigned sb,
                   const unsigned* ids, unsigned n, bm::serializer<BV>& normal,
                   bm::serializer<BV>& direct, typename bm::serializer<BV>::buffer& buffer,
                   bm::serialization_opcode_collector<BV>& collector)
{
    bm::serialization_opcode_profile base, alt;
    normal.serialize(group, buffer);
    collector.add_blob(buffer.data(), buffer.size(), base);
    direct.serialize(group, buffer);
    collector.add_blob(buffer.data(), buffer.size(), alt);
    if (!base.reconciles() || !alt.reconciles())
    {
        throw std::logic_error("Incomplete group byte accounting");
    }
    unsigned grouped = 0;
    for (unsigned code : {bm::set_sblock_bienc, bm::set_sblock_bienc_v3,
                          bm::set_sblock_bienc_gaps_v3, bm::set_sblock_bienc_wave_runs,
                          bm::set_sblock_bienc_adaptive,
                          bm::set_sblock_bienc_large})
    {
        if (alt.opcode[code].count)
        {
            throw std::logic_error("Superblock record in constituent control");
        }
        if (base.opcode[code].count)
        {
            grouped = code;
        }
    } // for code
    typename BV::statistics stats;
    group.calc_stat(&stats);
    unsigned occupied = 1, runs = 1, slot_count = 1, max_slot_count = 1;
    for (unsigned i = 1; i < n; ++i)
    {
        runs += ids[i] != ids[i - 1] + 1;
        if (ids[i] / bm::bits_in_block != ids[i - 1] / bm::bits_in_block)
        {
            ++occupied;
            slot_count = 0;
        }
        max_slot_count = std::max(max_slot_count, ++slot_count);
    } // for i
    std::cout << ordinal << '\t' << sb << '\t' << n << '\t' << occupied << '\t'
              << stats.bit_blocks << '\t' << stats.gap_blocks << '\t' << runs << '\t'
              << max_slot_count << '\t' << grouped << '\t' << base.bytes << '\t' << alt.bytes
              << '\t' << base.header_bytes + base.id_list_bytes
              << '\t' << alt.header_bytes + alt.id_list_bytes << '\t';
    write_records(base);
    std::cout << '\t';
    write_records(alt);
    std::cout << '\n';
}
} // namespace
#endif
