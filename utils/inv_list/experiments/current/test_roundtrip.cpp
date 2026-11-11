/** Induced in-memory round-trip failures using real serializer subclasses.
    Build with -I../../src -I. -pthread from utils/inv_list. Crypto is disabled
    to exercise the dependency-free round-trip path. No production fault flags.
*/
#define BM_NO_SHA256
#define main inv_list_application_main
#include "../../inv_list.cpp"
#undef main
#include <cassert>

/** @internal Emit a valid BLOB with altered contents to exercise exact equality. */
class changed_set_serializer : public bm::serializer<bm::bvector<>>
{
public:
    /** Serialize an altered copy of bv into borrowed output of size capacity. */
    size_type serialize(const bm::bvector<>& bv, unsigned char* output, size_t capacity)
    {
        bm::bvector<> changed(bv);
        changed.clear_bit(3);
        changed.set(5);
        return bm::serializer<bm::bvector<>>::serialize(changed, output, capacity);
    }
};

/** @internal Append an indexed trailing byte to exercise exact BLOB consumption. */
class extra_byte_serializer : public bm::serializer<bm::bvector<>>
{
public:
    /** Serialize bv plus one extra byte into borrowed output of size capacity. */
    size_type serialize(const bm::bvector<>& bv, unsigned char* output, size_t capacity)
    {
        const size_type bytes = bm::serializer<bm::bvector<>>::serialize(bv, output, capacity);
        assert(bytes < capacity);
        output[bytes] = 0;
        return bytes + 1;
    }
};

/** Assert that a worker's encoder failure survives future.get with level context.
    @tparam Serializer Real encoder or a fault-inducing subclass.
    @param expected Diagnostic substring expected from the induced failure.
*/
template<class Serializer>
void check_failure(const char* expected)
{
    inv_hash::Options options;
    options.budget = 8 * 1024 * 1024;
    options.threads = 1;
    options.start_level = options.end_level = 6;
    auto worker = std::async(std::launch::async, [&]()
    {
        inv_hash::BV bv;
        bv.set(1);
        bv.set(3);
        bv.set(7);
        inv_hash::Metrics metric;
        Serializer serializer;
        bm::streams_deserializer<inv_hash::BV> reader;
        std::unique_ptr<unsigned char[]> buffer;
        size_t capacity = 0;
        inv_hash::roundtrip_vector(bv, options, 17, metric, serializer, reader, buffer, capacity);
    });
    bool failed = false;
    try
    {
        worker.get();
    }
    catch (const std::runtime_error& error)
    {
        const std::string text(error.what());
        assert(text.find("roundtrip level=6") != std::string::npos);
        assert(text.find(expected) != std::string::npos);
        failed = true;
    }
    assert(failed);
}

/** Check semantic corruption, extra bytes and failure propagation from a worker. */
int main()
{
    check_failure<changed_set_serializer>("Rebuilt vector comparison failed");
    check_failure<extra_byte_serializer>("BLOB decode/length mismatch");
    return 0;
}
