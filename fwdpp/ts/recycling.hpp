#ifndef FWDPP_TS_RECYCLING_HPP
#define FWDPP_TS_RECYCLING_HPP

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <limits>
#include <cstdint>
#include <queue>
#include <fwdpp/forward_types.hpp>

namespace fwdpp
{
    namespace ts
    {
        std::queue<std::size_t>
        make_mut_queue(
            const std::vector<std::uint32_t> &mcounts,
            const std::vector<std::uint32_t> &counts_from_preserved_nodes)
        /// \brief Make a mutation recycling queue for simulations with tree sequences
        /// \param mcounts Contribution of extant nodes to mutation counts
        /// \param counts_from_preserved_nodes Contribution of extinct nodes to mutation counts
        ///
        /// \returns std::queue<std::size_t>
        {
            std::queue<std::size_t> mutation_recycling_bin;
            for (std::size_t i = 0; i < mcounts.size(); ++i)
                {
                    if (mcounts[i] + counts_from_preserved_nodes[i] == 0)
                        {
                            mutation_recycling_bin.push(i);
                        }
                }
            return mutation_recycling_bin;
        }

        namespace detail
        {
            template <typename mcont_t, typename lookup_table>
            void
            process_mutation_index(mcont_t &mutations, lookup_table &lookup,
                                   const std::size_t i)
            /// Implementation detail for flag_mutations_for_recycling
            {
                auto itr = lookup.equal_range(mutations[i].pos);
                mutations[i].pos = std::numeric_limits<double>::max();
                while (itr.first != itr.second)
                    {
                        if (itr.first->second == i)
                            {
                                lookup.erase(itr.first);
                                break;
                            }
                        ++itr.first;
                    }
            }

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t & /*fixations*/,
                              std::vector<uint_t> & /*fixation_times*/,
                              lookup_table &mut_lookup,
                              const uint_t /*generation*/, const std::size_t i,
                              std::false_type, std::false_type)
            // Remove all fixations.
            // Do not record fixations
            {
                process_mutation_index(mutations, mut_lookup, i);
                mcounts[i] = 0;
            }

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t &fixations,
                              std::vector<uint_t> &fixation_times,
                              lookup_table &mut_lookup,
                              const uint_t generation, const std::size_t i,
                              std::false_type, std::true_type)
            // Remove all fixations
            // Record fixations
            {
                fixations.push_back(mutations[i]);
                fixation_times.push_back(generation);
                process_mutation_index(mutations, mut_lookup, i);
                mcounts[i] = 0;
            }

            struct equal_range_comparison
            {
                template <typename mutation_type>
                inline bool
                operator()(const mutation_type &m, const double p) const
                {
                    return m.pos < p;
                }
                template <typename mutation_type>
                inline bool
                operator()(const double p, const mutation_type &m) const
                {
                    return p < m.pos;
                }
            };

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t & /*fixations*/,
                              std::vector<uint_t> & /*fixation_times*/,
                              lookup_table &mut_lookup,
                              const uint_t /*generation*/, const std::size_t i,
                              std::true_type, std::false_type)
            // Only remove neutral fixations
            // Do not record fixations
            {
                if (mutations[i].neutral)
                    {
                        process_mutation_index(mutations, mut_lookup, i);
                        mcounts[i] = 0;
                    }
            }

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t &fixations,
                              std::vector<uint_t> &fixation_times,
                              lookup_table &mut_lookup,
                              const uint_t generation, const std::size_t i,
                              std::true_type, std::true_type)
            // Only remove neutral fixations.
            // Record fixations.
            // The recording is in order of mutation position.
            {
                if (mutations[i].neutral)
                    {
                        // Record the fixation here,
                        // as we will recycle this variant
                        auto loc = std::lower_bound(
                            fixations.begin(), fixations.end(),
                            mutations[i].pos,
                            [](const typename mcont_t::value_type &m,
                               const double val) { return m.pos < val; });
                        auto d = std::distance(fixations.begin(), loc);
                        fixations.insert(loc, mutations[i]);
                        fixation_times.insert(fixation_times.begin() + d,
                                              generation);
                        process_mutation_index(mutations, mut_lookup, i);
                        mcounts[i] = 0;
                    }
                else
                    {
                        // The logic here is tougher.  We are preserving
                        // selected mutations, and thus we only
                        // want to record their fixations once. Thus,
                        // we have to make sure that this fixation
                        // hasn't been recorded at all.
                        auto loc_range = std::equal_range(
                            fixations.begin(), fixations.end(),
                            mutations[i].pos,
                            detail::equal_range_comparison());
                        if (std::find(loc_range.first, loc_range.second,
                                      mutations[i])
                            == loc_range.second)
                            {
                                fixations.insert(loc_range.first,
                                                 mutations[i]);
                                auto d = std::distance(fixations.begin(),
                                                       loc_range.first);
                                fixation_times.insert(
                                    fixation_times.begin() + d, generation);
                            }
                    }
            }
        } // namespace detail

        template <typename poptype, typename mutation_count_container,
                  typename preserve_selected_fixations,
                  typename record_fixations>
        void
        flag_mutations_for_recycling(
            poptype &pop,
            mutation_count_container &mcounts_from_preserved_nodes,
            const uint_t twoN, const uint_t generation,
            const preserve_selected_fixations preserve,
            const record_fixations record)
        /*! Mark mutations for recycling.
         *
         * This function should be called immediately after simplification and mutation counting.
         *
         * \param mutations A mutation container
         * \param mcounts A container stating how many times each element in \a mutations is present in the
         * currently-alive population
         * \param mcounts_from_preserved_nodes A container recording the counts of each mutation in ancient samples
         * \param lookup The population lookup table for mutations
         * \param twoN Twice the current population size
         * \param preserve_selected_fixations If true, do not mark selected fixations for recycling.
         *
         * A mutation is marked for recycling if one of the following conditions holds:
         * 1. The sum of \a mcounts and \a mcounts_from_preserved_nodes is zero.
         * 2. \a mcounts == \a twoN, mcounts_from_preserved_nodes is zero, 
         * and \a preserve_selected_fixations is false or the mutation is neutral.
         *
         * Condition 1 refers to extinct mutations and condition 2 refers to fixations.
         *
         * When \a preserve_selected_fixations is true, selected fixations are retained
         * in the population.  We do this because simulations of phenotypes (as opposed
         * to relative fitness) require tracking the contribution of fixation to 
         * trait values.
         *
         * All variants matching the above criteria have their record in \a mcounts
         * set to zero and their position is removed from \a lookup.
         *
         * This function does not record fixations/fixation times.  The reason is that
         * fixation times are only accurate if simplification happens very often.  A future
         * release will overload this function to handle that case, or you may write your own.
         *
         * Note that \a mutations is passed in non-const! There are guaranteed to be
         * no changes to the size of the container.  However, mutations marked for recycling
         * will have there positions change to numeric_limits<double>::max().
         *
         * \version 0.7.0 Added to library
         * \version 0.7.1 Updated to change recycled mutation positions to max value of a double.
         * Updated to take a population type as an argument.
         *
         * \todo Improve treatment of fixations by allowing for variants fixed in alive AND 
         * ancient samples to be flagged.
         * \todo Return a recycling queue?
         */
        {
            for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (mcounts_from_preserved_nodes[i] == 0)
                        {
                            if (pop.mcounts[i] > twoN)
                                {
                                    throw std::runtime_error(
                                        "mutation count out of range");
                                }
                            if (pop.mcounts[i] == twoN)
                                {
                                    detail::process_fixations(
                                        pop.mutations, pop.mcounts,
                                        pop.fixations, pop.fixation_times,
                                        pop.mut_lookup, generation, i,
                                        preserve, record);
                                }
                            else if (pop.mcounts[i] == 0)
                                {
                                    detail::process_mutation_index(
                                        pop.mutations, pop.mut_lookup, i);
                                }
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
