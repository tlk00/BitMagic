// Focused opt-in coverage for independent block selection.
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#ifdef TEST_BM64
#include "bm64.h"
#else
#include "bm.h"
#endif
#include "bmserial.h"
static void check(bool b) { if (!b) throw std::runtime_error("selection regression"); }
int main()
{
    using BV = bm::bvector<>;
    std::mt19937 rng(143);
    for (unsigned pattern = 0; pattern < 7; ++pattern)
    {
        BV v;
        if (pattern == 1) v.set_range(0, 65535);
        if (pattern == 2) { v.set(7); v.set(65540); v.optimize(); }
        if (pattern >= 3)
            for (unsigned i = 0; i < (pattern == 6 ? 140000u : 65536u); ++i)
                if (pattern == 3 ? i % 97 < 44 : pattern == 4 ? rng() % 100 < 85 : rng() % 100 < 35)
                    v.set(i);
        size_t sizes[9];
        bm::id64_t saved = 0;
        for (unsigned mode = 0; mode < 9; ++mode)
        {
            bm::serializer<BV> s; s.set_compression_level(6);
            if (mode == 1)
            {
                bm::periodic_options o;
                o.min_step = o.max_step = 1; o.order = 1; o.sample = 1;
                o.min_length = 4; o.max_holes = 1;
                o.max_chains = o.formula_chains = 256; o.trial_budget = 2304;
                s.set_periodic_options(o);
            }
            if (mode == 2) { bm::bwt_options o; o.diagnostics = true; s.set_bwt_full_options(o); }
            if (mode == 3) s.set_block_selection_options();
            if (mode == 4) s.set_block_selection_compose_options();
            if (mode == 5) s.set_block_selection_multistep_options(64, 2, 4);
            if (mode == 6) s.set_block_selection_step2_options();
            if (mode >= 3 && mode <= 6)
                s.set_block_selection_bwt_min_gain(0); // unrestricted controls
            if (mode == 7)
            {
                s.set_block_selection_step2_options();
                s.set_block_selection_bwt_min_gain(32);
            }
            if (mode == 8) s.set_block_selection_step2_options(); // default 16
            bm::serializer<BV>::buffer buf;
            s.serialize(v, buf); sizes[mode] = buf.size();
            BV decoded; bm::deserialize(decoded, buf.buf()); check(decoded == v);
            if (mode == 3)
            {
                const auto& stats = s.get_block_selection_statistics(); saved = stats.saved;
                for (size_t block_i = 0; block_i < stats.blocks.size(); ++block_i)
                {
                    const auto& b = stats.blocks[block_i];
                    check(b.selected == std::min(b.baseline, std::min(b.periodic, b.bwt)));
                    check(b.kind == (b.selected == b.baseline ? 0u : b.selected == b.periodic ? 1u : 2u));
                } // for block_i
            }
            if (mode == 4)
            {
                const auto& stats = s.get_block_selection_statistics();
                check(stats.composed);
                for (size_t block_i = 0; block_i < stats.blocks.size(); ++block_i)
                {
                    const auto& b = stats.blocks[block_i];
                    check(b.previous == std::min(b.baseline,
                                                 std::min(b.periodic, b.bwt)));
                    check(b.selected == std::min(b.previous, b.composed));
                    check(b.kind == (b.selected == b.baseline ? 0u :
                              b.selected == b.periodic ? 1u :
                              b.selected == b.bwt ? 2u : 3u));
                } // for block_i
                check(sizes[0] - sizes[mode] == stats.saved);
            }
            if (mode == 5)
            {
                const auto& stats = s.get_block_selection_statistics();
                check(stats.multistep && stats.composed);
                for (size_t block_i = 0; block_i < stats.blocks.size(); ++block_i)
                {
                    const auto& b = stats.blocks[block_i];
                    check(b.fourway == std::min(b.previous, b.composed));
                    check(b.selected == std::min(b.fourway, b.multistep));
                    check(b.kind == 4 ? b.multistep < b.fourway :
                                              b.selected == b.fourway);
                    check(b.multistep_last_step <= 64);
                    check(b.multistep_stop <= 3);
                    check(b.multistep_groups <= 64);
                } // for block_i
                check(sizes[0] - sizes[mode] == stats.saved);
            }
            if (mode == 6)
            {
                const auto& stats = s.get_block_selection_statistics();
                check(stats.step2 && stats.composed && !stats.multistep);
                for (size_t block_i = 0; block_i < stats.blocks.size(); ++block_i)
                {
                    const auto& b = stats.blocks[block_i];
                    check(b.fourway == std::min(b.previous, b.composed));
                    check(b.selected == std::min(b.fourway, b.step2));
                    check(b.kind == 4 ? b.step2 < b.fourway :
                                              b.selected == b.fourway);
                    check(!b.step2_chains || b.step2_formula + b.step2_residual + 2 == b.step2);
                } // for block_i
                check(sizes[0] - sizes[mode] == stats.saved);
            }
            if (mode == 7 || mode == 8)
            {
                const auto& stats = s.get_block_selection_statistics();
                const unsigned min_gain = mode == 7 ? 32u : 16u;
                check(stats.bwt_min_gain == min_gain);
                bm::id64_t rejected = 0, lost = 0;
                for (size_t block_i = 0; block_i < stats.blocks.size(); ++block_i)
                {
                    const auto& b = stats.blocks[block_i];
                    check(b.prethreshold_selected == std::min(b.fourway, b.step2));
                    if (b.prethreshold_kind >= 2 &&
                        b.periodic - b.prethreshold_selected < min_gain)
                    {
                        check(b.selected == b.periodic);
                        check(b.kind == (b.periodic < b.baseline ? 1u : 0u));
                        ++rejected;
                        lost += b.selected - b.prethreshold_selected;
                    }
                    else
                    {
                        check(b.selected == b.prethreshold_selected);
                        check(b.kind == b.prethreshold_kind);
                    }
                } // for block_i
                check(stats.threshold_rejected == rejected);
                check(stats.threshold_lost_bytes == lost);
                check(sizes[0] - sizes[mode] == stats.saved);
            }
        }
        check(sizes[3] <= sizes[1] && sizes[3] <= sizes[2]);
        check(sizes[0] - sizes[3] == saved);
        check(sizes[4] <= sizes[3]);
        check(sizes[5] <= sizes[4]);
        check(sizes[6] <= sizes[4]);
        check(sizes[7] >= sizes[6] && sizes[7] <= sizes[0]);
        check(sizes[8] >= sizes[6] && sizes[8] <= sizes[7]);
    }
    std::cout << "Block selection regression PASS\n";
}
