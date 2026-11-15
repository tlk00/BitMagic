// Measurement-only L8 wave features and exact padded record costs.
#include <fstream>
#include <iostream>
#include <array>
#include <vector>
#include "bm.h"
#include "bmserial.h"

using BV = bm::bvector<>;

/// Borrow the optimized source while its serializer calls the score observer.
struct Observation
{
    const BV* source;
    unsigned ordinal;
};

/// Export cheap source features and exact labels for one eligible superblock.
static void observe(void* context, const bm::sblock_wave_score_observation& score)
{
    auto& row = *static_cast<Observation*>(context);
    bm::bv_sub_survey s;
    s.init(nullptr, score.sb);
    row.source->get_blocks_manager().is_sparse_sblock(score.sb, 0, s);
    unsigned first = 256, last = 0, runs = 0, peak = 0, transitions = 0;
    std::array<unsigned, 4> levels{};
    for (unsigned j = 0; j < 256; ++j)
    {
        if (!s.bc_arr[j]) continue;
        first = std::min(first, j);
        last = j;
        runs += !j || !s.bc_arr[j-1];
        peak = std::max(peak, s.bc_arr[j]);
        if (BM_IS_GAP(s.blocks[j]))
            ++levels[bm::gap_level(BMGAP_PTR(s.blocks[j]))];
        else
            transitions += s.gc_arr[j];
    } // for j
    std::cout << row.ordinal << ' ' << score.sb << ' ' << score.n << ' '
              << score.blocks << ' ' << runs << ' ' << last-first+1 << ' '
              << peak << ' ' << s.gap_blocks << ' ' << s.gap_len_sum << ' ' << transitions;
    for (auto count : levels) std::cout << ' ' << count;
    std::cout << ' ' << score.incumbent << ' ' << std::min(score.incumbent,score.best_flat.bytes)
              << ' ' << std::min(score.incumbent,score.best_hier.bytes) << '\n';
}

/// Sample every 55th posting list across a raw Gov2 collection, including ordinal zero.
/// Arguments are the input collection and separate FULL-group output path.
int main(int argc, char** argv)
{
    if (argc != 3) return 2;
    std::ifstream in(argv[1], std::ios::binary);
    std::ofstream full(argv[2]);
    if (!in || !full) return 1;
    full << "vector sb full residual\n";
    std::cout << "vector sb n blocks runs span peak gap_blocks gap_words bit_changes gap0 gap1 gap2 gap3 rival flat hier\n";
    bm::serializer<BV> ser;
    ser.set_compression_level(8);
    ser.set_sblock_wave_experiment(1);
    Observation row{};
    ser.set_wave_score_observer(observe, &row);
    bm::serializer<BV>::buffer buf;
    BM_DECLARE_TEMP_BLOCK(tmp)
    unsigned n;
    for (unsigned i = 0; in.peek() != EOF; ++i)
    {
        if (!in.read(reinterpret_cast<char*>(&n), 4)) return 1;
        if (i % 55)
        {
            in.seekg(std::streamoff(n)*4, std::ios::cur);
            continue;
        }
        std::vector<unsigned> ids(n);
        if (!in.read(reinterpret_cast<char*>(ids.data()), std::streamsize(n)*4)) return 1;
        BV bv;
        if (n) bv.set(ids.data(), n, bm::BM_SORTED);
        bv.optimize(tmp);
        auto& manager = bv.get_blocks_manager();
        for (unsigned sb = 0; sb < manager.top_block_size(); ++sb)
        {
            bm::bv_sub_survey s{};
            s.init(nullptr, sb);
            auto root = manager.top_blocks_root();
            if (root && reinterpret_cast<const bm::word_t*>(root[sb]) == FULL_BLOCK_FAKE_ADDR)
                full << i << ' ' << sb << " 256 0\n";
            else
            {
                manager.is_sparse_sblock(sb, 0, s);
                if (s.full_blocks)
                    full << i << ' ' << sb << ' ' << s.full_blocks << ' '
                         << s.bv_count - 65536u*s.full_blocks << '\n';
            }
        } // for sb
        row = {&bv, i};
        ser.serialize(bv, buf);
    } // for i
    return in.bad() || !full || !std::cout;
}
