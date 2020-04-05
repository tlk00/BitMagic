#ifndef BM_DNA_FINGER_H__INCLUDED__
#define BM_DNA_FINGER_H__INCLUDED__

/**
    DNA finterprint scanner utility class.
*/

/**
    Utility for keeping all DNA finger print vectors
*/
template <typename BV>
class DNA_FingerprintScanner
{
public:
    enum { eA = 0, eC, eG, eT, eN, eEnd };

    typedef BV bvector_type;
    typedef typename bvector_type::size_type    size_type;
    typedef bm::aggregator<bvector_type>        aggregator_type;

public:
    DNA_FingerprintScanner() {}

    /// Build fingerprint bit-vectors using bulk insert iterator and parallel
    /// processing
    ///
    void BuildParallel(const std::vector<char>& sequence, unsigned threads)
    {
        struct Func
        {
            DNA_FingerprintScanner*      target_idx;
            const std::vector<char>*     src_sequence;

            Func(DNA_FingerprintScanner* idx, const std::vector<char>& src)
                : target_idx(idx), src_sequence(&src) {}

            void operator() (size_t from, size_t to)
            {
                const std::vector<char>& sequence = *src_sequence;
                bvector_type bvA, bvT, bvG, bvC, bvN;

                {
                    bm::bvector<>::bulk_insert_iterator iA(bvA, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iC(bvC, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iG(bvG, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iT(bvT, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iN(bvN, bm::BM_SORTED);
                    for (size_t i = from; i < sequence.size() && (i < to); ++i)
                    {
                        unsigned pos = unsigned(i);
                        switch (sequence[i])
                        {
                        case 'A':
                            iA = pos;
                            break;
                        case 'C':
                            iC = pos;
                            break;
                        case 'G':
                            iG = pos;
                            break;
                        case 'T':
                            iT = pos;
                            break;
                        case 'N':
                            iN = pos;
                            break;
                        default:
                            break;
                        }
                    } // for
                    // Bulk insert iterator keeps an buffer, which has to be
                    // flushed, before all bits appear in the target vector
                    //
                    iA.flush(); iC.flush(); iT.flush(); iG.flush(); iN.flush();
                }
                // merge results of parallel processing back to index
                target_idx->MergeVector('A', bvA);
                target_idx->MergeVector('T', bvT);
                target_idx->MergeVector('G', bvG);
                target_idx->MergeVector('C', bvC);
                target_idx->MergeVector('N', bvN);
            }
        };

        if (threads <= 1)
        {
            BuildBulk(sequence);
            return;
        }

        // Create parallel async tasks running on a range of source sequence
        //
        std::vector<std::future<void> > futures;
        futures.reserve(8);
        unsigned range = unsigned(sequence.size() / threads);

        for (unsigned k = 0; k < sequence.size(); k += range)
        {
            futures.emplace_back(std::async(std::launch::async,
                                 Func(this, sequence),  k, k + range));
        }

        // wait for all tasks
        for (auto& e : futures)
        {
            e.wait();
        }
    }

    /// Thread sync bit-vector merge
    ///
    void MergeVector(char letter, bm::bvector<>& bv)
    {
        static std::mutex mtx_A;
        static std::mutex mtx_T;
        static std::mutex mtx_G;
        static std::mutex mtx_C;
        static std::mutex mtx_N;

        switch (letter)
        {
        case 'A':
        {
            std::lock_guard<std::mutex> guard(mtx_A);
            m_FPrintBV[eA].merge(bv);
        }
        break;
        case 'C':
        {
            std::lock_guard<std::mutex> guard(mtx_C);
            m_FPrintBV[eC].merge(bv);
        }
        break;
        case 'G':
        {
            std::lock_guard<std::mutex> guard(mtx_G);
            m_FPrintBV[eG].merge(bv);
        }
        break;
        case 'T':
        {
            std::lock_guard<std::mutex> guard(mtx_T);
            m_FPrintBV[eT].merge(bv);
        }
        break;
        case 'N':
        {
            std::lock_guard<std::mutex> guard(mtx_N);
            m_FPrintBV[eN].merge(bv);
        }
        break;
        default:
            break;
        }

    }

    /// Build index using bulk insert iterator
    ///
    void BuildBulk(const std::vector<char>& sequence)
    {
        typedef typename bvector_type::bulk_insert_iterator bulk_inserter_type;

        bulk_inserter_type iA(m_FPrintBV[eA], bm::BM_SORTED);
        bulk_inserter_type iC(m_FPrintBV[eC], bm::BM_SORTED);
        bulk_inserter_type iG(m_FPrintBV[eG], bm::BM_SORTED);
        bulk_inserter_type iT(m_FPrintBV[eT], bm::BM_SORTED);
        bulk_inserter_type iN(m_FPrintBV[eN], bm::BM_SORTED);

        for (size_t i = 0; i < sequence.size(); ++i)
        {
            unsigned pos = unsigned(i);
            switch (sequence[i])
            {
            case 'A':
                iA = pos;
                break;
            case 'C':
                iC = pos;
                break;
            case 'G':
                iG = pos;
                break;
            case 'T':
                iT = pos;
                break;
            case 'N':
                iN = pos;
                break;
            default:
                break;
            }
        } // for i

        // flush inserters explicitly to avoid exceptions from destructors
        iA.flush(); iC.flush(); iG.flush(); iT.flush(); iN.flush();

    }


    /// Return fingerprint bit-vector
    const bvector_type& GetVector(char letter) const
    {
        switch (letter)
        {
        case 'A':
            return m_FPrintBV[eA];
        case 'C':
            return m_FPrintBV[eC];
        case 'G':
            return m_FPrintBV[eG];
        case 'T':
            return m_FPrintBV[eT];
        case 'N':
            return m_FPrintBV[eN];
        default:
            break;
        }
        throw std::runtime_error("Error. Invalid letter!");
    }


    /// This method uses cache blocked aggregator with fused SHIFT+AND kernel
    ///
    void Find(const std::string& word, std::vector<size_type>& res)
    {
        res.resize(0);

        if (word.empty())
            return;
        // first we setup aggregator, add a group of vectors to be processed
        m_Agg.reset();
        m_Agg.set_compute_count(false);
        for (size_t i = 0; i < word.size(); ++i)
        {
            const bm::bvector<>& bv_mask = GetVector(word[i]);
            m_Agg.add(&bv_mask);
        }

        // now run the whole algorithm to get benefits of cache blocking
        //
        bm::bvector<> bv;
        m_Agg.combine_shift_right_and(bv);

        // translate results from bvector of word ends to result
        unsigned ws = unsigned(word.size()) - 1;
        TranslateResults(bv, ws, res);

    };

    /// This method uses cache blocked aggregator with fused SHIFT+AND kernel
    ///
    /// @return search result count
    size_type FindCount(const std::string& word)
    {
        if (word.empty())
            return 0;
        // first we setup aggregator, add a group of vectors to be processed
        m_Agg.reset();
        m_Agg.set_compute_count(true);

        for (size_t i = 0; i < word.size(); ++i)
        {
            const bm::bvector<>& bv_mask = GetVector(word[i]);
            m_Agg.add(&bv_mask);
        }

        // now run the whole algorithm to get benefits of cache blocking
        //
        bm::bvector<> bv;
        m_Agg.combine_shift_right_and(bv);
        return m_Agg.count();
    };



    /// Find a set of words in one pass using pipeline
    /// of aggregators
    ///
    void FindCollection(const std::vector<std::tuple<std::string,int> >& words,
                        std::vector<std::vector<size_type>>& hits)
    {
        std::vector<std::unique_ptr<aggregator_type> > agg_pipeline;
        unsigned ws = 0;

        for (const auto& w : words)
        {
            std::unique_ptr<aggregator_type> agg_ptr(new aggregator_type());
            agg_ptr->set_operation(aggregator_type::BM_SHIFT_R_AND);

            const std::string& word = std::get<0>(w);
            for (size_t i = 0; i < word.size(); ++i)
            {
                const bm::bvector<>& bv_mask = GetVector(word[i]);
                agg_ptr->add(&bv_mask);
            }

            agg_pipeline.emplace_back(agg_ptr.release());
            ws = unsigned(word.size()) - 1;
        }

        // run the pipeline
        bm::aggregator_pipeline_execute<aggregator_type,
           std::vector<std::unique_ptr<aggregator_type> >::iterator>(
                                agg_pipeline.begin(), agg_pipeline.end()
                                );

        // convert the results
        for (size_t i = 0; i < agg_pipeline.size(); ++i)
        {
            const aggregator_type* agg_ptr = agg_pipeline[i].get();
            auto bv = agg_ptr->get_target();
            std::vector<size_type> res;
            res.reserve(12000);
            TranslateResults(*bv, ws, res);
            hits.emplace_back(res);
        }
    }

protected:

    /// Translate search results vector using (word size) left shift
    ///
    void TranslateResults(const bm::bvector<>& bv,
                          unsigned left_shift,
                          std::vector<size_type>& res)
    {
        typename bvector_type::enumerator en = bv.first();
        for (;en.valid(); ++en)
        {
            auto pos = *en;
            res.push_back(pos - left_shift);
        } // for en
    }

private:
    bvector_type     m_FPrintBV[eEnd];
    aggregator_type  m_Agg;
};


#endif

