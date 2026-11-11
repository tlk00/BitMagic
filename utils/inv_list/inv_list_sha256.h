#ifndef INV_LIST_SHA256_H_INCLUDED
#define INV_LIST_SHA256_H_INCLUDED
/* Copyright(c) 2026 Anatoliy Kuznetsov. Apache License, Version 2.0. */
/** @file inv_list_sha256.h
    @brief Application SHA-256 adapter for bm::bvector_hash and cache provenance.
    Uses CommonCrypto, Windows CNG or OpenSSL. BM_NO_SHA256 retains exact-only
    builds without crypto dependencies. No crypto implementation enters bm core.
*/
#include <cstddef>
#include <limits>
#include <new>
#include <stdexcept>
#if !defined(BM_NO_SHA256)
# if defined(__APPLE__)
#  include <CommonCrypto/CommonDigest.h>
# elif defined(_WIN32)
#  ifndef NOMINMAX
#   define NOMINMAX
#  endif
#  include <windows.h>
#  include <bcrypt.h>
#  if defined(_MSC_VER)
#   pragma comment(lib, "bcrypt.lib")
#  endif
# else
#  include <openssl/evp.h>
# endif
#endif

#if defined(__APPLE__) && defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
/** @internal Owning SHA-256 functor. One instance is used by one worker at a time. */
class inv_list_sha256
{
public:
    static const unsigned digest_size = 32; ///< Fixed SHA-256 output bytes.
    /** Allocate provider state and initialize a fresh stream; throws on failure. */
    inv_list_sha256();
    /** Release provider resources. */
    ~inv_list_sha256();
    inv_list_sha256(const inv_list_sha256&) = delete;
    inv_list_sha256& operator=(const inv_list_sha256&) = delete;
    /** @return Whether the executable was built with a SHA-256 provider. */
    static bool available();
    /** Reset the stream while retaining provider resources; errors propagate. */
    void reset();
    /** Append size bytes borrowed from data; data may be null only for size zero. */
    void update(const void* data, size_t size);
    /** Finalize into caller-owned digest_size bytes at output, once per reset. */
    void final(unsigned char* output);
private:
#if !defined(BM_NO_SHA256)
# if defined(__APPLE__)
    CC_SHA256_CTX context_; ///< CommonCrypto state.
# elif defined(_WIN32)
    BCRYPT_ALG_HANDLE algorithm_ = nullptr; ///< CNG algorithm provider.
    BCRYPT_HASH_HANDLE context_ = nullptr; ///< CNG hash object.
# else
    EVP_MD_CTX* context_ = nullptr; ///< OpenSSL streaming context.
# endif
#endif
};

inline
bool inv_list_sha256::available()
{
#ifdef BM_NO_SHA256
    return false;
#else
    return true;
#endif
}

inline
inv_list_sha256::inv_list_sha256()
{
#if !defined(BM_NO_SHA256) && defined(_WIN32)
    if (BCryptOpenAlgorithmProvider(&algorithm_, BCRYPT_SHA256_ALGORITHM, nullptr, 0) < 0)
    {
        throw std::runtime_error("Cannot open CNG SHA-256 provider");
    }
#elif !defined(BM_NO_SHA256) && !defined(__APPLE__)
    context_ = EVP_MD_CTX_new();
    if (!context_)
    {
        throw std::bad_alloc();
    }
#endif
    try
    {
        reset();
    }
    catch (...)
    {
#if !defined(BM_NO_SHA256) && defined(_WIN32)
        BCryptCloseAlgorithmProvider(algorithm_, 0);
#elif !defined(BM_NO_SHA256) && !defined(__APPLE__)
        EVP_MD_CTX_free(context_);
#endif
        throw;
    }
}

inline
inv_list_sha256::~inv_list_sha256()
{
#if !defined(BM_NO_SHA256) && defined(_WIN32)
    if (context_)
    {
        BCryptDestroyHash(context_);
    }
    BCryptCloseAlgorithmProvider(algorithm_, 0);
#elif !defined(BM_NO_SHA256) && !defined(__APPLE__)
    EVP_MD_CTX_free(context_);
#endif
}

inline
void inv_list_sha256::reset()
{
#if defined(BM_NO_SHA256)
    throw std::runtime_error("SHA-256 disabled in this build (BM_NO_SHA256); use exact verification or rebuild with a crypto provider");
#elif defined(__APPLE__)
    if (CC_SHA256_Init(&context_) != 1)
    {
        throw std::runtime_error("SHA-256 reset failed");
    }
#elif defined(_WIN32)
    if (context_)
    {
        BCryptDestroyHash(context_);
        context_ = nullptr;
    }
    if (BCryptCreateHash(algorithm_, &context_, nullptr, 0, nullptr, 0, 0) < 0)
    {
        throw std::runtime_error("CNG SHA-256 reset failed");
    }
#else
    if (EVP_DigestInit_ex(context_, EVP_sha256(), nullptr) != 1)
    {
        throw std::runtime_error("OpenSSL SHA-256 reset failed");
    }
#endif
}

inline
void inv_list_sha256::update(const void* data, size_t size)
{
#if defined(BM_NO_SHA256)
    (void)data;
    (void)size;
    reset();
#elif defined(__APPLE__) || defined(_WIN32)
    const auto* bytes = static_cast<const unsigned char*>(data);
    while (size)
    {
        const size_t maximum = (std::numeric_limits<unsigned long>::max)() < 0xfffffffful ?
            (std::numeric_limits<unsigned long>::max)() : 0xfffffffful;
        const unsigned long take = static_cast<unsigned long>(size < maximum ? size : maximum);
# if defined(__APPLE__)
        if (CC_SHA256_Update(&context_, bytes, CC_LONG(take)) != 1)
# else
        if (BCryptHashData(context_, const_cast<PUCHAR>(bytes), ULONG(take), 0) < 0)
# endif
        {
            throw std::runtime_error("SHA-256 update failed");
        }
        bytes += take;
        size -= take;
    } // while size
#else
    if (EVP_DigestUpdate(context_, data, size) != 1)
    {
        throw std::runtime_error("OpenSSL SHA-256 update failed");
    }
#endif
}

inline
void inv_list_sha256::final(unsigned char* output)
{
#if defined(BM_NO_SHA256)
    (void)output;
    reset();
#elif defined(__APPLE__)
    if (CC_SHA256_Final(output, &context_) != 1)
    {
        throw std::runtime_error("SHA-256 final failed");
    }
#elif defined(_WIN32)
    if (BCryptFinishHash(context_, output, digest_size, 0) < 0)
    {
        throw std::runtime_error("CNG SHA-256 final failed");
    }
#else
    unsigned size = 0;
    if (EVP_DigestFinal_ex(context_, output, &size) != 1 || size != digest_size)
    {
        throw std::runtime_error("OpenSSL SHA-256 final failed");
    }
#endif
}
#if defined(__APPLE__) && defined(__clang__)
#pragma clang diagnostic pop
#endif
#endif
