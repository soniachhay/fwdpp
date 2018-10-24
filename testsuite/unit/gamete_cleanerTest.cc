/*!
  \file gamete_cleanerTest.cc
  \ingroup unit
  \brief Testing fwdpp::fwdpp_internal::gamete_cleaner
*/
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/debug.hpp>
#include "../fixtures/fwdpp_fixtures.hpp"

BOOST_FIXTURE_TEST_SUITE(gamete_cleanerTest,
                         standard_empty_single_deme_fixture)

BOOST_AUTO_TEST_CASE(test_remove_all)
// TODO: Refactor this test to use a pop type
{
    diploids.resize(1000, std::make_pair(0, 0));
    mutations.emplace_back(mtype(0, 0, 0));
    mutations.emplace_back(mtype(1, 0, 0));
    mcounts.emplace_back(2000); // fixed
    mcounts.emplace_back(1);    // singleton
    gametes.emplace_back(gcont_t::value_type(0));
    gametes[0].mutations.push_back(0);
    gametes[0].n = 1999;
    gametes.push_back(gcont_t::value_type(1));
    gametes[1].mutations.push_back(1);
    int s = 0;
    for (auto& g : gametes)
        {
            s += g.n;
        }
    BOOST_REQUIRE_EQUAL(s, 2000);
    fwdpp::fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts, 2000,
                                          std::true_type());
    BOOST_REQUIRE_EQUAL(gametes[0].mutations.size(), 0);
}

BOOST_AUTO_TEST_CASE(test_remove_nothing)
// TODO: Refactor this test to use a pop type
{
    diploids.resize(1000, std::make_pair(0, 0));
    mutations.emplace_back(mtype(0, 0, 0));
    mutations.emplace_back(mtype(1, 0, 0));
    mcounts.emplace_back(2000); // fixed
    mcounts.emplace_back(1);    // singleton
    gametes.emplace_back(gcont_t::value_type(0));
    gametes[0].mutations.push_back(0);
    gametes[0].n = 1999;
    gametes.push_back(gcont_t::value_type(1));
    gametes[1].mutations.push_back(1);
    int s = 0;
    for (auto& g : gametes)
        {
            s += g.n;
        }
    BOOST_REQUIRE_EQUAL(s, 2000);
    fwdpp::fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts, 2000,
                                          fwdpp::remove_nothing());
    BOOST_REQUIRE_EQUAL(gametes[0].mutations.size(), 1);
}

BOOST_AUTO_TEST_CASE(test_remove_neutral)
// TODO: Refactor this test to use a pop type
{
    diploids.resize(1000, std::make_pair(0, 0));
    mutations.emplace_back(mtype(0, 0, 0));
    mutations.emplace_back(mtype(1, -0.1, 0)); // not neutral
    mcounts.emplace_back(2000);                // fixed
    mcounts.emplace_back(2000);                // fixed
    gametes.emplace_back(gcont_t::value_type(0));
    gametes[0].mutations.push_back(0);
    gametes[0].smutations.push_back(1);
    gametes[0].n = 2000;
    BOOST_REQUIRE_EQUAL(gametes[0].mutations.size(), 1);
    BOOST_REQUIRE_EQUAL(gametes[0].smutations.size(), 1);
    int s = 0;
    for (auto& g : gametes)
        {
            s += g.n;
        }
    BOOST_REQUIRE_EQUAL(s, 2000);
    fwdpp::fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts, 2000,
                                          fwdpp::remove_neutral());
    BOOST_REQUIRE_EQUAL(gametes[0].mutations.size(), 0);
    BOOST_REQUIRE_EQUAL(gametes[0].smutations.size(), 1);
}

BOOST_AUTO_TEST_SUITE_END()
