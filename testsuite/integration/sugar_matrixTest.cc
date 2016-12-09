/* \file sugar_matrixTest.cc
 * \brief Integration and unit tests of data matrix generation
 * \ingroup unit
 */

#include <config.h>
#include <boost/test/unit_test.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <gsl/gsl_matrix_char.h>

//This is an involved integration test of
//haplotype and genotype matrices, and it takes
//some time to run.
//We simulate a population of N=1000 diploids to equilibrium.
//We then construct a haplotype and genotype matrix from a large
//sample of diploids.  We tests that row and sum columns from the two
//matrix types agree.  Once those agreement tests pass, we compare the
//haplotype matrix to the row/sum columns of a sample based on 
//KTfwd::sample_separate (fwdpp/sugar/sampling.hpp) for the same diploids.
//These last checks ensure that two independent pieces of code, written 
//at different times, give the same results.
//
//Then, to really make sure, we do the above tests 1,000 times, evolving
//our population 100 generations in between each test.
BOOST_AUTO_TEST_CASE(singlepop_hapmatrix_exhaustive)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto generation = simulate_singlepop2(pop, rng, 0, 10000);
    std::vector<std::size_t> indlist;
    // Sample a LOT of individuals
    for (std::size_t i = 100; i < 750; i += 5)
        indlist.push_back(i);
    for (unsigned ntests = 0; ntests < 1000; ++ntests)
        {
            auto keys = mutation_keys(pop, indlist, true, true);
            // keys come out unsorted, so we have to sort for this test:
            std::sort(keys.first.begin(), keys.first.end(),
                      [&pop](const std::pair<std::size_t, KTfwd::uint_t> &a,
                             const std::pair<std::size_t, KTfwd::uint_t> &b) {
                          return pop.mutations[a.first].pos
                                 < pop.mutations[b.first].pos;
                      });
            std::sort(keys.second.begin(), keys.second.end(),
                      [&pop](const std::pair<std::size_t, KTfwd::uint_t> &a,
                             const std::pair<std::size_t, KTfwd::uint_t> &b) {
                          return pop.mutations[a.first].pos
                                 < pop.mutations[b.first].pos;
                      });
            // In order to be comparable to a sample taken from the population,
            // we must remove sites fixed in these inidivudals
            keys.first.erase(
                std::remove_if(
                    keys.first.begin(), keys.first.end(),
                    [&indlist](
                        const std::pair<std::size_t, KTfwd::uint_t> &a) {
                        return a.second == 2 * indlist.size();
                    }),
                keys.first.end());
            keys.second.erase(
                std::remove_if(
                    keys.second.begin(), keys.second.end(),
                    [&indlist](
                        const std::pair<std::size_t, KTfwd::uint_t> &a) {
                        return a.second == 2 * indlist.size();
                    }),
                keys.second.end());
            // Create data matrices and do basic sanity checks
            auto m = haplotype_matrix(pop, indlist, keys.first, keys.second);
            BOOST_REQUIRE_EQUAL(m.nrow, 2 * indlist.size());
            BOOST_REQUIRE_EQUAL(m.neutral.size(), m.nrow * keys.first.size());
            BOOST_REQUIRE_EQUAL(m.selected.size(),
                                m.nrow * keys.second.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_popfreq.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(),
                                m.selected_positions.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_popfreq.size());
			
			//Note that we can convert the data to vectors of other types
			//as needed
			std::vector<double> neutral_as_double(m.neutral.begin(),m.neutral.end());
			for(std::size_t i=0;i<m.neutral.size();++i)
			{
				BOOST_REQUIRE_EQUAL(m.neutral[i],neutral_as_double[i]);
			}

            auto gm = genotype_matrix(pop, indlist, keys.first, keys.second);
            BOOST_REQUIRE_EQUAL(gm.nrow, indlist.size());
            BOOST_REQUIRE_EQUAL(gm.neutral.size(),
                                gm.nrow * keys.first.size());
            BOOST_REQUIRE_EQUAL(gm.selected.size(),
                                gm.nrow * keys.second.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(),
                                gm.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(), gm.neutral_popfreq.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(),
                                gm.selected_positions.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(),
                                gm.selected_popfreq.size());

            // Now, compare to an independent calculation of the genotypes
            // for same individuals
            auto s = KTfwd::sample_separate(pop, indlist, true);
            // Check same total # of mutations
            BOOST_REQUIRE_EQUAL(s.first.size(), m.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(s.second.size(), m.selected_positions.size());
            BOOST_REQUIRE_EQUAL(s.first.size(), gm.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(s.second.size(), gm.selected_positions.size());

            std::vector<double> pos;
            for (auto &&i : s.first)
                pos.push_back(i.first);
            // Make sure all positions come out ok
            for (std::size_t i = 0; i < pos.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(pos[i], m.neutral_positions[i]);
                    BOOST_REQUIRE_EQUAL(pos[i], gm.neutral_positions[i]);
                }
            pos.clear();
            for (auto &&i : s.second)
                pos.push_back(i.first);
            for (std::size_t i = 0; i < pos.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(pos[i], m.selected_positions[i]);
                    BOOST_REQUIRE_EQUAL(pos[i], gm.selected_positions[i]);
                }
            auto sums = KTfwd::col_sums(m);
            auto gm_sums = KTfwd::col_sums(gm);
            BOOST_REQUIRE_EQUAL(sums.first.size(), m.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(sums.second.size(),
                                m.selected_positions.size());
            // Check that haplotype and genotype matrices
            // have same row sums in same order
            for (std::size_t i = 0; i < sums.first.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(sums.first[i], gm_sums.first[i]);
                }
            for (std::size_t i = 0; i < sums.second.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(sums.second[i], gm_sums.second[i]);
                }
            // Make sure column sums check out for neutral sites.
            // Only compare to haplotype matrix. If we've made it this far,
            // we know m == gm in terms of column sums.
            for (std::size_t c = 0; c < sums.first.size(); ++c)
                {
                    unsigned sum2 = 0;
                    for (std::size_t i = 0; i < m.nrow; ++i)
                        {
                            sum2 += (s.first[c].second[i] == '1') ? 1 : 0;
                        }
                    BOOST_REQUIRE_EQUAL(sums.first[c], sum2);
                }
            // Make sure column sums check out for selected sites
            for (std::size_t c = 0; c < sums.second.size(); ++c)
                {
                    unsigned sum2 = 0;
                    for (std::size_t i = 0; i < m.nrow; ++i)
                        {
                            sum2 += (s.second[c].second[i] == '1') ? 1 : 0;
                        }
                    BOOST_REQUIRE_EQUAL(sums.second[c], sum2);
                }
            // Now, row sums.
            sums = KTfwd::row_sums(m);
            gm_sums = KTfwd::row_sums(gm);
            // Compare row sums from haplotype and genotype matrix.
            for (std::size_t hi = 0, gi = 0; hi < sums.first.size();
                 hi += 2, gi += 1)
                {
                    // Neutral sites
                    BOOST_REQUIRE_EQUAL(sums.first[hi] + sums.first[hi + 1],
                                        gm_sums.first[gi]);
                }
            for (std::size_t hi = 0, gi = 0; hi < sums.second.size();
                 hi += 2, gi += 1)
                {
                    // Selected sites
                    BOOST_REQUIRE_EQUAL(sums.second[hi] + sums.second[hi + 1],
                                        gm_sums.second[gi]);
                }
            // Check neutral sites in haplotype matrix to the sample
            for (std::size_t r = 0; r < sums.first.size(); ++r)
                {
                    unsigned sum2 = 0;
                    for (std::size_t i = 0; i < m.neutral_positions.size();
                         ++i)
                        {
                            sum2 += (s.first[i].second[r] == '1') ? 1 : 0;
                        }
                    BOOST_REQUIRE_EQUAL(sums.first[r], sum2);
                }
            // Check selected sites
            for (std::size_t r = 0; r < sums.second.size(); ++r)
                {
                    unsigned sum2 = 0;
                    for (std::size_t i = 0; i < m.selected_positions.size();
                         ++i)
                        {
                            sum2 += (s.second[i].second[r] == '1') ? 1 : 0;
                        }
                    BOOST_REQUIRE_EQUAL(sums.second[r], sum2);
                }
            /*
if (!m.selected_positions.empty())
{
v = gsl_matrix_short_const_view_array(
m.selected.data(), m.nrow,
m.selected_positions.size());
BOOST_REQUIRE_EQUAL(v.matrix.size1, 2 * indlist.size());
// Make sure column sums check out
for (std::size_t c = 0; c < v.matrix.size2; ++c)
{
auto col_view
    = gsl_matrix_short_const_column(&v.matrix, c);
unsigned sum = 0, sum2 = 0;
for (std::size_t i = 0; i < v.matrix.size1; ++i)
    {
        sum += gsl_vector_short_get(
            &col_view.vector, i);
        sum2 += (s.second[c].second[i] == '1') ? 1
                                               : 0;
    }
BOOST_REQUIRE_EQUAL(sum, sum2);
}
// Now, row sums.
for (std::size_t r = 0; r < v.matrix.size1; ++r)
{
auto row_view
    = gsl_matrix_short_const_row(&v.matrix, r);
unsigned sum = 0, sum2 = 0;
for (std::size_t i = 0; i < v.matrix.size2; ++i)
    {
        sum += gsl_vector_short_get(
            &row_view.vector, i);
        sum2 += (s.second[i].second[r] == '1') ? 1
                                               : 0;
    }
BOOST_REQUIRE_EQUAL(sum, sum2);
}
}
*/
            // Look at properties of the GENOTYPE matrix
            m = genotype_matrix(pop, indlist, keys.first, keys.second);
            // Evolve pop 100 more generations
            generation = simulate_singlepop2(pop, rng, generation, 100);
        }
}

BOOST_AUTO_TEST_CASE(singlepop_hapmatrix_GSL_behavior)
// What happens if try to create empty matrix?
// This can happen in production code, and
// the default behavior is that GSL raises
// a signal. There are two solutions:
// 1. Don't make empty views
// 2. Turn off error handler, check return value,
// and then turn error handler back on.
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000); // a monomorphic population
    std::vector<std::size_t> indlist;
    // Sample a LOT of individuals
    for (std::size_t i = 100; i < 750; i += 5)
        indlist.push_back(i);
    auto keys = mutation_keys(pop, indlist, true, true);
    auto m = haplotype_matrix(pop, indlist, keys.first, keys.second);
    BOOST_REQUIRE(keys.first.empty());
    auto error_handler = gsl_set_error_handler_off();
    auto v = gsl_matrix_char_const_view_array(m.neutral.data(), m.nrow,
                                               m.neutral_positions.size());
    BOOST_REQUIRE_EQUAL(v.matrix.size1, 0);
    BOOST_REQUIRE_EQUAL(v.matrix.size2, 0);
    gsl_set_error_handler(error_handler);
}
