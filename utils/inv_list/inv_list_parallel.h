#ifndef INV_LIST_PARALLEL_H_INCLUDED
#define INV_LIST_PARALLEL_H_INCLUDED

// Ordered, bounded producer/worker/writer path shared by compression policies.
// Included by inv_list.cpp after its collection I/O helpers and options.
#include <memory>
#include <sstream>
#include <limits>
#include <cerrno>
#include <cstring>
#include "inv_list_progress.h"
#include "inv_list_reports.h"
#include "inv_list_read_buffer.h"
#include "inv_list_worker_error.h"

struct plain_compression_policy
{
    struct statistics {};
    static const char* name() { return c_level >= 7 ? "compression-level-preset" : "baseline"; }
    static void configure(bm::serializer<bm::bvector<> >&) {}
    static statistics snapshot(const bm::serializer<bm::bvector<> >&) { return {}; }
    static void add(statistics&, const statistics&) {}
    static void progress(std::ostream&, const statistics&) {}

};

struct block_selection_step2_policy
{
    typedef bm::serializer<bm::bvector<> >::block_selection_statistics_type statistics;
    static const char* name() { return "fourway-or-BWT-periodic2"; }
    static void configure(bm::serializer<bm::bvector<> >& s)
        { s.set_block_selection_step2_options();
          if (bwt_min_gain_specified) s.set_block_selection_bwt_min_gain(bwt_min_gain);
          if (periodic_step1_interleaved) s.set_periodic_step1_interleaved(periodic_step1_interleaved); }
    static statistics snapshot(const bm::serializer<bm::bvector<> >& s)
        { return s.get_block_selection_statistics(); }
    static void add(statistics& d, const statistics& s) { bm::block_selection_add(d, s); }
    static void progress(std::ostream& out, const statistics& s)
    {
        inv_report_block_selection_step2_policy(out, s);
    }

};

struct gap_xor_compression_policy : plain_compression_policy
{
    static const char* name() { return "step2-with-GAP-XOR"; }
    static void configure(bm::serializer<bm::bvector<> >& s)
    {
        block_selection_step2_policy::configure(s);
        s.set_gap_xor_options(true);
    }
};

struct gap_route_compression_policy
{
    typedef bm::gap_route_statistics statistics;
    static const char* name() { return "step2-GAP-route-survey"; }
    static void configure(bm::serializer<bm::bvector<> >& s)
    {
        block_selection_step2_policy::configure(s);
        s.set_gap_route_diagnostics(true);
    }
    static statistics snapshot(const bm::serializer<bm::bvector<> >& s)
        { return s.get_gap_route_statistics(); }
    static void add(statistics& d, const statistics& s)
        { bm::gap_route_add(d, s); }
    static void progress(std::ostream& out, const statistics& s)
    {
        inv_report_gap_route_compression_policy(out, s);
    }

};

struct sblock_adaptive_compression_policy : plain_compression_policy
{
    static const char* name() { return "step2-adaptive-superblock"; }
    static void configure(bm::serializer<bm::bvector<> >& s)
    {
        block_selection_step2_policy::configure(s);
        s.set_sblock_adaptive_options(true);
    }
};

struct sblock_adaptive_gap_xor_policy : sblock_adaptive_compression_policy
{
    static const char* name() { return "step2-adaptive-superblock-with-GAP-XOR"; }
    static void configure(bm::serializer<bm::bvector<> >& s)
    {
        sblock_adaptive_compression_policy::configure(s);
        s.set_gap_xor_options(true);
    }
};

template<class Policy>
struct parallel_compression_job
{
    std::vector<unsigned> input;
    bm::serializer<bm::bvector<> >::buffer blob;
    typename Policy::statistics stats;
    bm::id64_t min_delta_count = 0;
    bm::id64_t vector_index = 0;
    size_t input_count = 0;
    std::streamoff input_end = 0;
    double seconds = 0;
    std::exception_ptr error; ///< Worker exception, published by task completion.
};

template<class Policy>
static void compress_inv_dump_file_parallel(const std::string& fname,
                                             const std::string& out_fname)
{
    typedef bm::bvector<> bv_type;
    typedef bm::thread_pool<bm::task_descr*,
                            bm::spin_lock<bm::pad0_struct> > pool_type;
    typedef bm::task_batch<bv_type::allocator_type> batch_type;
    typedef std::chrono::steady_clock clock_type;
    const size_t max_input_bytes = 512ull * 1024 * 1024;

    cout << "Reading input collection: " << fname << endl;
    cout << "Writing to BV collection: " << out_fname << endl;
    cout << "Compression level: " << c_level << endl;
    if (compression_diagnostics) cout << "Parallel compression: experiment=" << Policy::name()
         << " threads=" << compression_threads
         << " batch_vectors=" << compression_batch_vectors
         << " input_window_mib=512" << endl;
    bm::chrono_taker tt1(std::cout, "1. Convert collection", 1, &timing_map);
    std::ifstream fin(fname.c_str(), std::ios::in | std::ios::binary);
    if (!fin) throw std::runtime_error("Cannot open input file");
    fin.seekg(0, std::ios::end);
    const std::streamoff fsize = fin.tellg();
    fin.seekg(0, std::ios::beg);
    std::ofstream out(out_fname.c_str(), std::ios::out | std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open bvect out file");

    std::vector<std::unique_ptr<parallel_compression_job<Policy> > > jobs;
    jobs.reserve(compression_batch_vectors);
    batch_type batch;
    batch.get_task_vector().reserve(compression_batch_vectors);
    pool_type pool;
    pool.start(compression_threads);
    bm::id64_t total_ints = 0, total_min_delta = 0;
    bm::id64_t prepared = 0, emitted = 0, output_bytes = 0;
    std::vector<bm::id64_t> blob_offsets;
    if (!bv_index_out_file.empty()) blob_offsets.push_back(0);
    typename Policy::statistics totals;
    bm::serialization_opcode_profile opcode_profile;
    std::unique_ptr<bm::serialization_opcode_collector<bv_type> > opcode_collector;
    if (profile_opcodes) {
        opcode_collector.reset(new bm::serialization_opcode_collector<bv_type>);
        opcode_collector->set_collect_ids();
    }
    std::ofstream survey;
    if (!opcode_survey_file.empty()) {
        survey.open(opcode_survey_file);
        if (!survey) throw std::runtime_error("Cannot open opcode survey output");
        survey << "vector_index,source_offset,vector_ids,vector_bytes,opcode,records,bytes,ids\n";
    }
    size_t pending_input_bytes = 0, peak_input_bytes = 0;
    double worker_seconds = 0, head_wait_seconds = 0;
    const auto start = clock_type::now();

    inv_list_progress_timer progress_timer(progress_vectors != 0 && !progress_silent, start);
    std::streamoff input_consumed = 0;
    double last_vector_seconds = 0;
    size_t last_vector_bytes = 0;
    auto report_progress = [&](bool count_due)
    {
        const auto now = clock_type::now();
        if (!progress_timer.due(count_due, now))
        {
            return;
        }
        const double elapsed = std::chrono::duration<double>(now - start).count();
        const inv_progress_values values{emitted, prepared,
            inv_list_completed_tasks(batch, prepared), input_consumed,
            output_bytes, elapsed, last_vector_seconds, last_vector_bytes};
        std::ostringstream progress;
        inv_report_compress_progress(progress, values);
        Policy::progress(progress, totals);
        if (opcode_collector)
        {
            progress << " profiled_records=" << opcode_profile.records;
        }
        std::cerr << progress.str() << std::endl;
        progress_timer.reported(now);
    };

    auto drain = [&]()
    {
        if (jobs.empty()) return;
        for (size_t i = 0; i < jobs.size(); ++i)
        {
            const auto wait_begin = compression_diagnostics ? clock_type::now() : clock_type::time_point{};
            inv_list_wait_for_task(pool, batch, unsigned(i), progress_vectors != 0 && !progress_silent,
                                   [&]() { report_progress(false); });
            if (compression_diagnostics) head_wait_seconds += std::chrono::duration<double>(clock_type::now() - wait_begin).count();
            parallel_compression_job<Policy>& job = *jobs[i];
            const int task_code = batch.get_task(unsigned(i))->err_code;
            if (job.error || task_code)
            {
                std::ostringstream context;
                context << "Compression worker failed: input=" << fname
                        << " output=" << out_fname << " vector_index=" << job.vector_index;
                inv_list_check_worker(job.error, task_code, context.str());
            }
            if (!job.blob.size()) throw std::runtime_error("Empty worker compression BLOB");
            errno = 0;
            out.write(reinterpret_cast<const char*>(job.blob.data()),
                      std::streamsize(job.blob.size()));
            const int write_errno = errno;
            if (!out)
            {
                std::ostringstream message;
                message << "Error writing BV collection: path=" << out_fname
                        << " vector_index=" << job.vector_index
                        << " completed_vectors=" << emitted
                        << " output_offset=" << output_bytes
                        << " requested_bytes=" << job.blob.size()
                        << " stream_state=" << unsigned(out.rdstate())
                        << " errno=" << write_errno;
                if (write_errno) message << " (" << std::strerror(write_errno) << ')';
                throw std::runtime_error(message.str());
            }
            if (opcode_collector) {
                if (opcode_survey_file.empty())
                    opcode_collector->add_blob(job.blob.data(), job.blob.size(), opcode_profile);
                else {
                    bm::serialization_opcode_profile local;
                    opcode_collector->add_blob(job.blob.data(), job.blob.size(), local);
                    if (!local.reconciles() || !local.cardinalities_reconcile() || local.ids != job.input_count)
                        throw std::logic_error("Per-vector survey accounting failed");
                    const auto offset=job.input_end-std::streamoff(job.input_count*4+4);
                    auto row = [&](int code, bm::id64_t count, bm::id64_t bytes, bm::id64_t ids) {
                        survey << job.vector_index << ',' << offset << ',' << job.input_count
                               << ',' << job.blob.size() << ',' << code << ',' << count
                               << ',' << bytes << ',' << ids << '\n';
                    };
                    row(-1,0,local.header_bytes,0);
                    if (local.id_list_bytes) row(-2,0,local.id_list_bytes,local.id_list_ids);
                    for (unsigned code=0;code<256;++code) {
                        const auto& e=local.opcode[code];
                        if(e.count) row(int(code),e.count,e.bytes,e.ids);
                    }
                    if (!survey) throw std::runtime_error("Error writing opcode survey");
                    opcode_profile.add(local);
                }
            }
            output_bytes += job.blob.size();
            if (!bv_index_out_file.empty()) blob_offsets.push_back(output_bytes);
            total_ints += job.input_count;
            total_min_delta += job.min_delta_count;
            if (compression_diagnostics) worker_seconds += job.seconds;
            Policy::add(totals, job.stats);
            ++emitted;
            input_consumed = job.input_end;
            last_vector_seconds = job.seconds;
            last_vector_bytes = job.blob.size();
            const bool count_due = progress_vectors &&
                (emitted % progress_vectors == 0 ||
                 (max_vectors && emitted == max_vectors) || job.input_end == fsize);
            report_progress(count_due);
        }
        jobs.clear();
        batch.get_task_vector().clear();
        pending_input_bytes = 0;
    };

    try
    {
        while (!max_vectors || prepared < max_vectors)
        {
            std::unique_ptr<parallel_compression_job<Policy> > job(
                new parallel_compression_job<Policy>);
            if (io_read_u32_coll(fin, job->input))
                throw std::runtime_error("Error reading input file");
            unsigned min_delta = 0, min_delta_cnt = 0;
            if (validate_inp_vec(job->input, min_delta, min_delta_cnt) ||
                !min_delta || !min_delta_cnt)
                throw std::runtime_error("Input vector validation failed");
            job->min_delta_count = min_delta_cnt;
            job->vector_index = prepared;
            job->input_count = job->input.size();
            job->input_end = fin.tellg();
            const bool at_eof = job->input_end == fsize;
            const size_t input_bytes = job->input.size() * sizeof(unsigned) + sizeof(unsigned);
            if (!jobs.empty() &&
                (jobs.size() == compression_batch_vectors ||
                 pending_input_bytes + input_bytes > max_input_bytes))
                drain();
            pending_input_bytes += input_bytes;
            if (pending_input_bytes > peak_input_bytes) peak_input_bytes = pending_input_bytes;
            jobs.push_back(std::move(job));
            parallel_compression_job<Policy>* target = jobs.back().get();
            bm::task_function_t task([target](void*) -> int
            {
                try
                {
                    const auto begin = clock_type::now();
                    BM_DECLARE_TEMP_BLOCK(tb)
                    bv_type bv;
                    bv.set(target->input.data(), bv_type::size_type(target->input.size()), bm::BM_SORTED);
                    bv.optimize(tb);
                    bm::serializer<bv_type> serializer;
                    configure_inv_serializer(serializer, inv_serializer_options(), Policy::configure);
                    serializer.set_periodic_vector_index(target->vector_index);
                    serializer.set_bwt_vector_index(target->vector_index);
                    serializer.serialize(bv, target->blob);
                    if (!opcode_survey_file.empty()) {
                        bm::serialization_profile_streambuf view(target->blob.data(),target->blob.size());
                        std::istream stream(&view); bm::streams_decoder input(stream);
                        bm::streams_deserializer<bv_type> reader; bv_type restored;
                        if (!reader.deserialize(restored,input) || input.tell()!=target->blob.size() || restored!=bv)
                            throw std::runtime_error("Survey raw-vector round trip failed");
                    }

                    target->stats = Policy::snapshot(serializer);
                    target->seconds = std::chrono::duration<double>(clock_type::now() - begin).count();
                    target->input.clear();
                    target->input.shrink_to_fit();
                    return 0;
                }
                catch (...)
                {
                    target->error = std::current_exception();
                    return 1;
                }
            });
            batch.add(task, 0);
            pool.get_job_queue().push(batch.get_task(unsigned(jobs.size() - 1)));
            ++prepared;
            report_progress(false);
            if (at_eof) break;
        }
        drain();
    }
    catch (...)
    {
        pool.set_stop_mode(pool_type::stop_when_done);
        pool.join();
        throw;
    }
    pool.set_stop_mode(pool_type::stop_when_done);
    pool.join();
    if (!opcode_survey_file.empty()) {
        survey.close();
        if (!survey) throw std::runtime_error("Error closing opcode survey");
        cout << "Opcode survey verified_vectors=" << emitted << " rows_file=" << opcode_survey_file << endl;
    }
    errno = 0;
    out.close();
    const int close_errno = errno;
    if (!out)
    {
        std::ostringstream message;
        message << "Error closing BV collection: path=" << out_fname
                << " completed_vectors=" << emitted
                << " output_bytes=" << output_bytes
                << " stream_state=" << unsigned(out.rdstate())
                << " errno=" << close_errno;
        if (close_errno) message << " (" << std::strerror(close_errno) << ')';
        throw std::runtime_error(message.str());
    }
    if (!bv_index_out_file.empty())
    {
        std::ofstream index(bv_index_out_file.c_str(), std::ios::out | std::ios::trunc);
        if (!index) throw std::runtime_error("Cannot open BLOB index output");
        index << "BM_INV_BLOB_INDEX_V1 " << output_bytes << ' ' << emitted << '\n';
        for (bm::id64_t offset : blob_offsets) index << offset << '\n';
        index.close();
        if (!index) throw std::runtime_error("Error writing BLOB index");
    }

    if (opcode_collector)
    {
        if (!opcode_profile.reconciles() || !opcode_profile.cardinalities_reconcile() ||
            opcode_profile.ids != total_ints || opcode_profile.bytes != output_bytes ||
            opcode_profile.blobs != emitted)
            throw std::runtime_error("Opcode profile does not match collection bytes");
        bm::print_serialization_opcode_profile(cout, opcode_profile);
    }

    if (prepared!=emitted) throw std::logic_error("Prepared/emitted vector count mismatch");
    inv_report_compression_summary(cout, compression_diagnostics, prepared, emitted, peak_input_bytes, worker_seconds, head_wait_seconds, fin.tellg(), fsize, total_ints, total_min_delta, output_bytes);

}

// The sidecar marks BLOB boundaries, which cannot be inferred from a byte
// offset without parsing the preceding self-delimiting BLOB.
static std::vector<bm::id64_t> read_inv_blob_offsets(
    const std::string& path, bm::id64_t file_bytes)
{
    std::ifstream input(path.c_str());
    std::string magic;
    bm::id64_t indexed_bytes = 0, count = 0;
    if (!(input >> magic >> indexed_bytes >> count) ||
        magic != "BM_INV_BLOB_INDEX_V1" || indexed_bytes != file_bytes ||
        !count || count > 100000000)
        throw std::runtime_error("Invalid or stale BLOB index header");
    std::vector<bm::id64_t> offsets;
    offsets.reserve(size_t(count + 1));
    for (bm::id64_t i = 0; i <= count; ++i)
    {
        bm::id64_t offset;
        if (!(input >> offset) || offset > file_bytes ||
            (i && offset <= offsets.back()) || (!i && offset))
            throw std::runtime_error("Invalid BLOB index offset");
        offsets.push_back(offset);
    }
    std::string extra;
    if (offsets.back() != file_bytes || (input >> extra))
        throw std::runtime_error("BLOB index does not cover exact file");
    return offsets;
}

struct parallel_read_job
{
    std::vector<unsigned> input;
    bm::id64_t begin = 0, end = 0;
    const unsigned char* blob = nullptr; ///< Borrowed batch slice, valid until all tasks join.
    std::streamoff input_end = 0;
    size_t input_count = 0;
    double seconds = 0;
    std::exception_ptr error; ///< Worker exception, published by task completion.
    bm::deserializer<bm::bvector<>, bm::decoder>::sblock_records_vect_type sblock_records;
};

static unsigned sblock_bic_length_bin(unsigned length)
{
    if (!length || length > 65535)
        throw std::runtime_error("Invalid BIC superblock length");
    if (length <= 16) return length - 1;
    unsigned exponent = 0;
    for (unsigned n = length; n >>= 1; ) ++exponent;
    const unsigned base = 1u << exponent;
    const unsigned quarter = (length - base) * 4 / base;
    return 16 + (exponent - 4) * 4 + quarter;
}

static void read_inv_dump_file_parallel(const std::string& input_name,
                                        const std::string& bv_name, bool verify)
{
    typedef bm::bvector<> bv_type;
    typedef bm::thread_pool<bm::task_descr*,
                            bm::spin_lock<bm::pad0_struct> > pool_type;
    typedef bm::task_batch<bv_type::allocator_type> batch_type;
    typedef std::chrono::steady_clock clock_type;
    bm::chrono_taker timer(std::cout, verify ? "2. Verify collection" :
                                      "3. Decode collection", 1, &timing_map);
    std::ifstream bv_file(bv_name.c_str(), std::ios::binary);
    if (!bv_file) throw std::runtime_error("Cannot open bvect dump file");
    bv_file.seekg(0, std::ios::end);
    const std::streamoff file_bytes = bv_file.tellg();
    if (file_bytes <= 0) throw std::runtime_error("Empty or unreadable BV collection");
    bv_file.seekg(0);
    inv_list_read_buffer encoded_buffer(bv_name, read_mmap);
    const auto offsets = read_inv_blob_offsets(bv_index_in_file,
                                                bm::id64_t(file_bytes));
    const bm::id64_t count = offsets.size() - 1;
    if (verify && max_vectors && max_vectors != count)
        throw std::runtime_error("BLOB index count differs from -max-vectors");
    std::ifstream fin;
    std::streamoff input_size = 0;
    if (verify)
    {
        fin.open(input_name.c_str(), std::ios::binary);
        if (!fin) throw std::runtime_error("Cannot open input file");
        fin.seekg(0, std::ios::end);
        input_size = fin.tellg(); fin.seekg(0);
        if (input_size < 0) throw std::runtime_error("Cannot size input file");
    }
    cout << "Parallel " << (verify ? "verification" : "decode")
         << ": threads=" << compression_threads
         << " batch_vectors=" << compression_batch_vectors
         << " indexed_vectors=" << count
         << " read_memory_mib=" << read_memory_bytes / (1024 * 1024)
         << " read_backend=" << (read_mmap ? "mmap" : "buffered") << endl;
    std::vector<std::unique_ptr<parallel_read_job> > jobs;
    jobs.reserve(compression_batch_vectors);
    batch_type batch;
    batch.get_task_vector().reserve(compression_batch_vectors);
    pool_type pool;
    pool.start(compression_threads);
    const auto start = clock_type::now();
    bm::id64_t prepared = 0, emitted = 0, total_ints = 0;
    std::vector<bm::id64_t> sblock_exact(sblock_len_stats ? 65536 : 0, 0);
    bm::id64_t sblock_bin_count[64] = {}, sblock_bin_bytes[64] = {};
    bm::id64_t sblock_count = 0, sblock_bytes = 0, sblock_length_sum = 0;
    unsigned sblock_min = 65536, sblock_max = 0;
    size_t staged_bytes = 0, peak_staged_bytes = 0;
    double worker_seconds = 0;
    std::streamoff input_consumed = 0;
    inv_list_progress_timer progress_timer(progress_vectors != 0 && !progress_silent, start);
    bm::id64_t encoded_consumed = 0;
    auto report_progress = [&](bool count_due)
    {
        const auto now = clock_type::now();
        if (!progress_timer.due(count_due, now))
        {
            return;
        }
        const double elapsed = std::chrono::duration<double>(now - start).count();
        const inv_progress_values values{emitted, prepared,
            batch.size() ? inv_list_completed_tasks(batch, prepared) : emitted, input_consumed,
            encoded_consumed, elapsed, 0, 0};
        std::ostringstream progress;
        inv_report_read_progress(progress, values, verify);
        std::cerr << progress.str() << std::endl;
        progress_timer.reported(now);
    };
    auto drain = [&]()
    {
        for (size_t i = 0; i < jobs.size(); ++i)
        {
            inv_list_wait_for_task(pool, batch, unsigned(i), progress_vectors != 0 && !progress_silent,
                                   [&]() { report_progress(false); });
            const parallel_read_job& job = *jobs[i];
            const int task_code = batch.get_task(unsigned(i))->err_code;
            if (job.error || task_code)
            {
                std::ostringstream context;
                context << (verify ? "Verification" : "Decode")
                        << " worker failed: input=" << bv_name
                        << " raw_input=" << input_name
                        << " vector_index=" << emitted
                        << " blob_begin=" << job.begin << " blob_end=" << job.end;
                inv_list_check_worker(job.error, task_code, context.str());
            }
            total_ints += job.input_count;
            if (sblock_len_stats)
                for (const auto& rec : job.sblock_records)
                {
                    const unsigned bin = sblock_bic_length_bin(rec.length);
                    ++sblock_exact[rec.length];
                    ++sblock_bin_count[bin];
                    sblock_bin_bytes[bin] += rec.bytes;
                    ++sblock_count;
                    sblock_bytes += rec.bytes;
                    sblock_length_sum += rec.length;
                    sblock_min = std::min(sblock_min, rec.length);
                    sblock_max = std::max(sblock_max, rec.length);
                }
            if (compression_diagnostics) worker_seconds += job.seconds;
            input_consumed = job.input_end;
            ++emitted;
            encoded_consumed = job.end;
            const bool count_due = progress_vectors &&
                (emitted % progress_vectors == 0 || emitted == count);
            report_progress(count_due);
        }
        jobs.clear(); batch.get_task_vector().clear(); staged_bytes = 0;
    };
    try
    {
        unsigned pending_length = 0;
        while (prepared < count)
        {
            const bm::id64_t batch_begin = offsets[size_t(prepared)];
            bm::id64_t batch_end = batch_begin;
            bm::id64_t raw_bytes = 0;
            while (prepared < count && jobs.size() < compression_batch_vectors)
            {
                const bm::id64_t end = offsets[size_t(prepared + 1)];
                if (verify && !pending_length)
                {
                    fin.read(reinterpret_cast<char*>(&pending_length), sizeof(pending_length));
                    if (!fin || !pending_length)
                    {
                        throw std::runtime_error("Invalid raw vector length");
                    }
                }
                const bm::id64_t next_raw = bm::id64_t(pending_length) * sizeof(unsigned);
                if (next_raw > read_memory_bytes - raw_bytes ||
                    end - batch_begin > read_memory_bytes - raw_bytes - next_raw)
                {
                    if (jobs.empty())
                    {
                        throw std::runtime_error("Vector exceeds read batch budget; increase -read-memory-mib");
                    }
                    break;
                }
                std::unique_ptr<parallel_read_job> job(new parallel_read_job);
                job->begin = offsets[size_t(prepared)];
                job->end = end;
                if (verify)
                {
                    job->input.resize(pending_length);
                    fin.read(reinterpret_cast<char*>(job->input.data()), std::streamsize(next_raw));
                    if (!fin)
                    {
                        throw std::runtime_error("Truncated raw input vector");
                    }
                    pending_length = 0;
                    job->input_count = job->input.size();
                    job->input_end = fin.tellg();
                }
                raw_bytes += next_raw;
                batch_end = end;
                jobs.push_back(std::move(job));
                ++prepared;
                report_progress(false);
            } // while prepared
            staged_bytes = size_t(raw_bytes + batch_end - batch_begin);
            peak_staged_bytes = std::max(peak_staged_bytes, staged_bytes);
            encoded_buffer.load(bv_file, batch_begin, batch_end - batch_begin);
            report_progress(false);
            for (size_t i = 0; i < jobs.size(); ++i)
            {
                parallel_read_job* target = jobs[i].get();
                target->blob = encoded_buffer.data() + size_t(target->begin - batch_begin);
                bm::task_function_t task([target, verify](void*) -> int
                {
                    try
                    {
                        const auto begin = clock_type::now();
                        const size_t blob_size = size_t(target->end - target->begin);
                        bm::serialization_profile_streambuf view(target->blob, blob_size);
                        std::istream stream(&view);
                        bm::streams_decoder channel(stream);
                        bm::streams_deserializer<bv_type> reader;
                        if (sblock_len_stats)
                        {
                            reader.set_sblock_bic_records(&target->sblock_records);
                        }
                        bv_type bv;
                        if (!reader.deserialize(bv, channel) || channel.tell() != blob_size)
                        {
                            throw std::runtime_error("Indexed BLOB decode or length mismatch");
                        }
                        if (verify && compare_vect(target->input, bv))
                        {
                            throw std::runtime_error("Vector comparison failed");
                        }
                        target->seconds = std::chrono::duration<double>(clock_type::now() - begin).count();
                        return 0;
                    }
                    catch (...)
                    {
                        target->error = std::current_exception();
                        return 1;
                    }
                });
                batch.add(task, 0);
                pool.get_job_queue().push(batch.get_task(unsigned(i)));
            } // for i
            drain();
            encoded_buffer.clear();
        } // while prepared
    }
    catch (...)
    {
        pool.set_stop_mode(pool_type::stop_when_done);
        pool.join();
        throw;
    }
    pool.set_stop_mode(pool_type::stop_when_done);
    pool.join();
    if (verify && input_consumed != input_size && !max_vectors)
        throw std::runtime_error("Input collection has extra vectors");
    inv_report_read_summary(cout, verify, prepared, emitted, worker_seconds, peak_staged_bytes, sblock_len_stats, sblock_exact, sblock_bin_count, sblock_bin_bytes, sblock_count, sblock_bytes, sblock_length_sum, sblock_min, sblock_max, input_consumed, input_size, total_ints);

}

#endif
