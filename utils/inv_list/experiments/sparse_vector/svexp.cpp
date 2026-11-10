// Gov2 sparse-vector collection experiment. Build with make svexp-neon.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "bm.h"
#include "bmserial.h"
#include "bmfio.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmthreadpool.h"
#include "bmundef.h"

namespace {
using sv_type = bm::sparse_vector<unsigned, bm::bvector<> >;
using clock_type = std::chrono::steady_clock;
using pool_type = bm::thread_pool<bm::task_descr*, bm::spin_lock<bm::pad0_struct> >;
using batch_type = bm::task_batch<bm::bvector<>::allocator_type>;

struct options {
    std::string action, input, stream, index, representation, codec;
    unsigned threads = 10, batch = 512, progress = 1000;
    uint64_t max_vectors = 0;
};

uint64_t number(const char* str) {
    std::string s(str);
    if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("Invalid unsigned number: " + s);
    size_t used = 0;
    const auto n = std::stoull(s, &used);
    if (used != s.size()) throw std::runtime_error("Invalid unsigned number: " + s);
    return n;
}

unsigned small_number(const char* str) {
    const uint64_t value = number(str);
    if (value > std::numeric_limits<unsigned>::max())
        throw std::runtime_error("Number exceeds unsigned range");
    return unsigned(value);
}

options parse(int argc, char** argv) {
    if (argc < 2) throw std::runtime_error("Usage: svexp encode|verify|decode --stream FILE --index FILE --representation raw|gaps --codec plain|xor|select [--input FILE] [--threads N] [--batch N] [--max-vectors N] [--progress N]");
    options o; o.action = argv[1];
    for (int i = 2; i < argc; i += 2) {
        if (i + 1 >= argc) throw std::runtime_error("Missing option value");
        const std::string key(argv[i]), value(argv[i+1]);
        if (key == "--input") o.input = value;
        else if (key == "--stream") o.stream = value;
        else if (key == "--index") o.index = value;
        else if (key == "--representation") o.representation = value;
        else if (key == "--codec") o.codec = value;
        else if (key == "--threads") o.threads = small_number(argv[i+1]);
        else if (key == "--batch") o.batch = small_number(argv[i+1]);
        else if (key == "--max-vectors") o.max_vectors = number(argv[i+1]);
        else if (key == "--progress") o.progress = small_number(argv[i+1]);
        else throw std::runtime_error("Unknown option: " + key);
    }
    if ((o.action != "encode" && o.action != "verify" && o.action != "decode") ||
        o.stream.empty() || o.index.empty() ||
        (o.action != "decode" && o.input.empty()) ||
        (o.representation != "raw" && o.representation != "gaps") ||
        (o.codec != "plain" && o.codec != "xor" && o.codec != "select") ||
        !o.threads || o.threads > 256 || !o.batch || o.batch > 512)
        throw std::runtime_error("Invalid or incomplete svexp options");
    return o;
}

bool read_vector(std::ifstream& file, std::vector<unsigned>& data) {
    unsigned n = 0;
    file.read(reinterpret_cast<char*>(&n), sizeof(n));
    if (!file || !n || n > 100000000u) return false;
    data.resize(n);
    file.read(reinterpret_cast<char*>(data.data()), std::streamsize(size_t(n)*sizeof(unsigned)));
    if (!file) throw std::runtime_error("Truncated Gov2 vector");
    for (size_t i = 1; i < data.size(); ++i)
        if (data[i] <= data[i-1]) throw std::runtime_error("Gov2 vector is not strictly increasing");
    return true;
}

unsigned min_gap(const std::vector<unsigned>& input) {
    if (input.size() < 2) return 0;
    unsigned result = std::numeric_limits<unsigned>::max();
    for (size_t i = 1; i < input.size(); ++i)
        result = std::min(result, input[i] - input[i-1]);
    return result;
}

void make_sparse(const std::vector<unsigned>& input, bool gaps, sv_type& sv) {
    auto out = sv.get_back_inserter();
    if (gaps) {
        const unsigned minimum = min_gap(input);
        out.add(minimum);
        out.add(input.front());
        for (size_t i = 1; i < input.size(); ++i)
            out.add(input[i] - input[i-1] - minimum);
    } else {
        for (unsigned value : input) out.add(value);
    }
    out.flush();
    BM_DECLARE_TEMP_BLOCK(tb)
    sv.optimize(tb);
}

void check_sparse(const sv_type& sv, const std::vector<unsigned>& input, bool gaps) {
    if (sv.size() != input.size() + (gaps ? 1u : 0u))
        throw std::runtime_error("Sparse vector length mismatch: got " +
                                 std::to_string(sv.size()) + " expected " +
                                 std::to_string(input.size() + (gaps ? 1u : 0u)));
    if (gaps) {
        const unsigned minimum = sv.get(0);
        uint64_t value = sv.get(1);
        if (value != input.front()) throw std::runtime_error("Sparse first ID mismatch");
        for (size_t i = 1; i < input.size(); ++i) {
            value += uint64_t(minimum) + sv.get(bm::id_t(i+1));
            if (value != input[i]) throw std::runtime_error("Sparse gap reconstruction mismatch");
        }
    } else {
        for (size_t i = 0; i < input.size(); ++i)
            if (sv.get(bm::id_t(i)) != input[i])
                throw std::runtime_error("Sparse raw ID mismatch");
    }
}

struct job {
    std::vector<unsigned> input;
    std::string bytes, error;
    uint64_t begin = 0, end = 0;
    uint64_t input_end = 0;
    unsigned selected[5] = {};
    uint64_t selected_saved = 0;
    uint64_t xor_blocks = 0;
    double seconds = 0;
};

std::vector<uint64_t> read_index(const options& o, uint64_t length) {
    std::ifstream file(o.index);
    std::string magic, repr, codec;
    uint64_t bytes = 0, count = 0;
    if (!(file >> magic >> repr >> codec >> bytes >> count) ||
        magic != "BM_SV_GOV2_INDEX_V1" || repr != o.representation ||
        codec != o.codec || bytes != length || !count || count > 100000000)
        throw std::runtime_error("Invalid or mismatched sparse-vector index");
    std::vector<uint64_t> offsets;
    offsets.reserve(size_t(count+1));
    for (uint64_t i = 0; i <= count; ++i) {
        uint64_t offset;
        if (!(file >> offset) || offset > length || (i && offset <= offsets.back()) ||
            (!i && offset)) throw std::runtime_error("Invalid sparse-vector offset");
        offsets.push_back(offset);
    }
    std::string extra;
    if (offsets.back() != length || file >> extra ||
        (o.max_vectors && count != o.max_vectors))
        throw std::runtime_error("Sparse-vector index count or end mismatch");
    return offsets;
}

void execute(const options& o) {
    const bool encode = o.action == "encode", verify = o.action == "verify";
    const bool gaps = o.representation == "gaps";
    std::ifstream input;
    uint64_t input_bytes = 0;
    if (encode || verify) {
        input.open(o.input, std::ios::binary);
        if (!input) throw std::runtime_error("Cannot open Gov2 input");
        input.seekg(0, std::ios::end);
        input_bytes = uint64_t(input.tellg());
        input.seekg(0);
    }
    std::ofstream output;
    std::ifstream stored;
    std::vector<uint64_t> offsets(1, 0);
    if (encode) {
        output.open(o.stream, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("Cannot open sparse-vector output");
    } else {
        stored.open(o.stream, std::ios::binary);
        if (!stored) throw std::runtime_error("Cannot open sparse-vector input");
        stored.seekg(0, std::ios::end);
        const auto size = stored.tellg();
        if (size <= 0) throw std::runtime_error("Empty sparse-vector collection");
        offsets = read_index(o, uint64_t(size));
    }
    pool_type pool;
    pool.start(o.threads);
    batch_type batch;
    batch.get_task_vector().reserve(o.batch);
    std::vector<std::unique_ptr<job> > jobs;
    jobs.reserve(o.batch);
    const auto start = clock_type::now();
    uint64_t prepared = 0, completed = 0, total_ints = 0, output_bytes = 0;
    uint64_t input_consumed = 0;
    uint64_t selected[5] = {}, selected_saved = 0, xor_blocks = 0;
    size_t pending_input = 0;
    double worker_seconds = 0;
    auto drain = [&]() {
        for (size_t i = 0; i < jobs.size(); ++i) {
            bm::thread_pool_executor<pool_type>::wait_for_batch_done(pool, batch, unsigned(i), unsigned(i));
            job& j = *jobs[i];
            if (!j.error.empty()) throw std::runtime_error(j.error);
            if (encode) {
                if (j.bytes.empty()) throw std::runtime_error("Empty sparse-vector record");
                output.write(j.bytes.data(), std::streamsize(j.bytes.size()));
                if (!output) throw std::runtime_error("Error writing sparse-vector collection");
                output_bytes += j.bytes.size();
                offsets.push_back(output_bytes);
            } else output_bytes = j.end;
            total_ints += j.input.size();
            if (encode || verify) input_consumed = j.input_end;
            worker_seconds += j.seconds;
            selected_saved += j.selected_saved;
            xor_blocks += j.xor_blocks;
            for (unsigned k = 0; k < 5; ++k) selected[k] += j.selected[k];
            ++completed;
            if (o.progress && (completed % o.progress == 0 ||
                               (o.max_vectors && completed == o.max_vectors) ||
                               (!encode && completed == offsets.size()-1))) {
                const double elapsed = std::chrono::duration<double>(clock_type::now()-start).count();
                std::cerr << "Progress " << o.action << " vectors=" << completed
                          << " input_bytes=" << input_consumed
                          << " output_bytes=" << output_bytes << " elapsed_s=" << std::fixed
                          << std::setprecision(2) << elapsed << " vectors_per_s="
                          << (elapsed > 0 ? double(completed)/elapsed : 0) << '\n';
            }
        }
        jobs.clear(); batch.get_task_vector().clear(); pending_input = 0;
    };
    try {
        while (true) {
            if (o.max_vectors && prepared == o.max_vectors) break;
            if (!encode && prepared == offsets.size()-1) break;
            std::unique_ptr<job> next(new job);
            if (encode || verify) {
                if (input.peek() == EOF) {
                    if (encode && !o.max_vectors) break;
                    throw std::runtime_error("Gov2 input ended before expected vector count");
                }
                if (!read_vector(input, next->input))
                    throw std::runtime_error("Invalid Gov2 record");
                next->input_end = uint64_t(input.tellg());
            }
            if (!encode) {
                next->begin = offsets[size_t(prepared)];
                next->end = offsets[size_t(prepared+1)];
            }
            const size_t bytes = next->input.size()*sizeof(unsigned);
            if (!jobs.empty() && (jobs.size() == o.batch ||
                                  pending_input + bytes > 512ull*1024*1024)) drain();
            pending_input += bytes;
            jobs.push_back(std::move(next));
            job* target = jobs.back().get();
            bm::task_function_t task([target, &o, encode, verify, gaps](void*) -> int {
                try {
                    const auto begin = clock_type::now();
                    if (encode) {
                        sv_type sv;
                        make_sparse(target->input, gaps, sv);
                        bm::sparse_vector_serializer<sv_type> writer;
                        writer.disable_xor_compression();
                        auto& bvs = writer.get_bv_serializer();
                        bvs.set_compression_level(6);
                        if (o.codec == "xor") writer.enable_xor_compression();
                        if (o.codec == "select") {
                            bvs.set_block_selection_step2_options();
                            bvs.set_block_selection_bwt_min_gain(16);
                            bvs.set_periodic_step1_interleaved(4);
                        }
                        std::ostringstream stream(std::ios::out | std::ios::binary);
                        bm::streams_encoder channel(stream);
                        if (!writer.serialize(sv, channel) || !channel.finish())
                            throw std::runtime_error("Sparse-vector serialization failed");
                        target->bytes = stream.str();
                        if (o.codec == "xor") {
                            const auto* counts = bvs.get_compression_stat();
                            target->xor_blocks = counts[bm::set_block_ref_eq] +
                                                 counts[bm::set_block_xor_ref32] +
                                                 counts[bm::set_block_xor_chain];
                        }
                        if (o.codec == "select") {
                            const auto& s = bvs.get_block_selection_statistics();
                            for (unsigned k = 0; k < 5; ++k) target->selected[k] = unsigned(s.chosen[k]);
                            target->selected_saved = s.saved;
                        }
                    } else {
                        std::ifstream stream(o.stream, std::ios::binary);
                        if (!stream) throw std::runtime_error("Cannot open indexed sparse-vector stream");
                        stream.seekg(std::streamoff(target->begin));
                        sv_type sv;
                        bm::streams_sparse_vector_deserializer<sv_type> reader;
                        if (!reader.deserialize(sv, stream) ||
                            stream.tellg() != std::streamoff(target->end))
                            throw std::runtime_error("Sparse-vector decode or record boundary mismatch");
                        if (verify) check_sparse(sv, target->input, gaps);
                    }
                    target->seconds = std::chrono::duration<double>(clock_type::now()-begin).count();
                    return 0;
                } catch (const std::exception& e) { target->error = e.what(); return 1; }
                  catch (...) { target->error = "Unknown sparse-vector worker failure"; return 1; }
            });
            batch.add(task, 0);
            pool.get_job_queue().push(batch.get_task(unsigned(jobs.size()-1)));
            ++prepared;
        }
        drain();
    } catch (...) {
        pool.set_stop_mode(pool_type::stop_when_done);
        pool.join();
        throw;
    }
    pool.set_stop_mode(pool_type::stop_when_done);
    pool.join();
    if (!completed) throw std::runtime_error("No sparse vectors processed");
    if (encode) {
        output.close();
        if (!output) throw std::runtime_error("Sparse-vector output close failed");
        std::ofstream index(o.index, std::ios::trunc);
        if (!index) throw std::runtime_error("Cannot write sparse-vector index");
        index << "BM_SV_GOV2_INDEX_V1 " << o.representation << ' ' << o.codec << ' '
              << output_bytes << ' ' << completed << '\n';
        for (auto offset : offsets) index << offset << '\n';
        index.close();
        if (!index) throw std::runtime_error("Sparse-vector index close failed");
    } else if (verify && !o.max_vectors && input_consumed != input_bytes)
        throw std::runtime_error("Gov2 input has extra vectors");
    const double elapsed = std::chrono::duration<double>(clock_type::now()-start).count();
    std::cout << "Result action=" << o.action << " representation=" << o.representation
              << " codec=" << o.codec << " vectors=" << completed
              << " original_ints=" << total_ints << " stream_bytes=" << output_bytes
              << " elapsed_s=" << elapsed << " worker_s=" << worker_seconds;
    if (o.codec == "select" && encode)
        std::cout << " ordinary=" << selected[0] << " periodic=" << selected[1]
                  << " bwt=" << selected[2] << " composed=" << selected[3]
                  << " step2=" << selected[4] << " selected_saved=" << selected_saved;
    if (o.codec == "xor" && encode)
        std::cout << " xor_blocks=" << xor_blocks;
    std::cout << '\n';
}
}

int main(int argc, char** argv) {
    try { execute(parse(argc, argv)); return 0; }
    catch (const std::exception& e) { std::cerr << "svexp: " << e.what() << '\n'; return 1; }
}
