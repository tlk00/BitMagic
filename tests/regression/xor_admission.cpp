/** \file xor_admission.cpp
    Standalone regression for internal level-dependent XOR admission.
    Build with the repository src include directory.
*/
#include <iostream>
#include <stdexcept>
#include "bm.h"
#include "bmserial.h"

/// Check internal overrides and level policy using medium-density XOR residuals.
int main()
{
    bm::bvector<> source, reference;
    unsigned state = 123456789;
    auto random = [&]()
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    };
    for (unsigned i = 0; i < 65536; ++i)
    {
        bool ref = (random() & 1) != 0;
        bool residual = random() % 100 < 38;
        if (ref)
        {
            reference.set(i);
        }
        if (ref != residual)
        {
            source.set(i);
        }
    } // for i
    bm::bv_ref_vector<bm::bvector<>> refs;
    refs.add(&source, 0);
    refs.add(&reference, 1);
    bm::serializer<bm::bvector<>> writer;
    auto match = [&](unsigned level, unsigned limit)
    {
        writer.set_compression_level(level);
        bm::xor_sim_params params;
        params.admission_max_metric = limit;
        bm::xor_sim_model<bm::bvector<>> model;
        if (!writer.compute_sim_model(model, refs, params))
        {
            throw std::runtime_error("Model not built");
        }
        return model.matr.get(0, 0).match != bm::e_no_xor_match;
    };
    if (match(6, 0) || match(7, 21845) || match(10, 21845))
    {
        throw std::runtime_error("Legacy policy unexpectedly admitted residual");
    }
    for (unsigned level = 7; level <= 10; ++level)
    {
        if (!match(level, 0) || !match(level, 29789))
        {
            throw std::runtime_error("Relaxed policy rejected residual");
        }
    } // for level
    if (!match(6, 29789) || match(6, 0))
    {
        throw std::runtime_error("Override or level reset failed");
    }
    std::cout << "Policy and override checks passed\n";
}
