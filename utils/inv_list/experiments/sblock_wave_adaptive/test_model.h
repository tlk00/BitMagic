// Independent pre-cleanup trial encoder, retained only to verify dry scores.
#ifndef TEST_SBLOCK_WAVE_MODEL_H
#define TEST_SBLOCK_WAVE_MODEL_H
namespace control_wave_test {
inline void require(bool v) { if(!v) throw std::logic_error("Invalid test model"); }
struct trial_model {
    unsigned wave_size=1024,wave_count=16384;
    std::vector<unsigned> occupied, empty, values;
    unsigned correction = 0, upper = 0;

    bool prepare(const unsigned* original, unsigned n, bm::sblock_wave_record& r)
    {
        r = bm::sblock_wave_record(); occupied.clear(); empty.clear(); values.resize(n);
        require(n >= 2);
        for (unsigned i = 0; i < n; ++i)
        {
            unsigned wave = original[i] / wave_size;
            if (occupied.empty() || occupied.back() != wave) occupied.push_back(wave);
            values[i] = unsigned(occupied.size() - 1) * wave_size + original[i] % wave_size;
        }
        if (occupied.size() == wave_count) return false;
        for (unsigned wave = 0, j = 0; wave < wave_count; ++wave)
            if (j < occupied.size() && occupied[j] == wave) ++j;
            else empty.push_back(wave);
        r.eligible = true;
        correction = values[1] - values[0] - 1;
        for (unsigned i = 2; i < n; ++i)
            if (values[i] - values[i-1] - 1 < correction)
                correction = values[i] - values[i-1] - 1;
        upper = unsigned(occupied.size()) * wave_size - 1 - (n - 1) * correction;
        for (unsigned i = 0; i < n; ++i) values[i] -= i * correction;
        require(values.back() <= upper);
        return true;
    }

    // One accumulator spans header, map, bounds and position payload.
    template<class BitOut, class Position>
    void encode(BitOut& out, Position position, unsigned sb, unsigned mode,
                bm::sblock_wave_cost& cost) const
    {
        const unsigned n = unsigned(values.size());
        const auto& map_values = (mode & 2) ? empty : occupied;
        out.put_bits(mode, 2); // polarity and explicit-boundary flag
        unsigned sb_mode = sb < 8 ? 0 : sb <= 255 ? 1 : sb <= 65535 ? 2 : 3;
        out.put_bits(sb_mode, 2);
        switch (sb_mode)
        {
        case 0: out.gamma(sb + 1); break;
        case 1: out.put_bits(sb, 8); break;
        case 2: out.put_16_no(sb); break;
        default: out.put_32_no(sb); break;
        }
        out.gamma(n); out.gamma(unsigned(map_values.size())); out.gamma(correction + 1);
        cost.header = position(); // includes placeholder opcode
        out.bic_encode_u32_cm(map_values.data(), unsigned(map_values.size()), 0, wave_count - 1);
        cost.map = position() - cost.header;
        if (mode & 1)
        {
            const bm::id64_t before = position();
            out.gamma(values.front() + 1);
            out.gamma(upper - values.back() + 1);
            cost.header += position() - before;
        }
        const bm::id64_t before = position();
        if (mode & 1)
        {
            if (n > 2) out.bic_encode_u32_cm(values.data() + 1, n - 2,
                                           values.front() + 1, values.back() - 1);
        }
        else out.bic_encode_u32_cm(values.data(), n, 0, upper);
        cost.bic = position() - before;
    }

};
}
#endif
