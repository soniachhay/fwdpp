#ifndef FWDPP_POPBASE_HPP
#define FWDPP_POPBASE_HPP

#include <type_traits>
#include <vector>
#include <exception>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace fwdpp
{
    namespace poptypes
    {
        template <typename mutation_type, typename mcont, typename gcont,
                  typename mvector, typename ftvector, typename lookup_table_type>
        class popbase
        /*!
          \brief Base class for population objects
          \note Added in fwdpp 0.5.0.  Changed in 0.6.0 to be independent of ploidy.

          \version 0.9.0 Remove construction from containers of low-level types
         */
        {
            static_assert(typename fwdpp::traits::is_haploid_genome<
                              typename gcont::value_type>::type(),
                          "gcont::value_type must be a haploid_genome type");
            static_assert(
                typename fwdpp::traits::is_mutation<typename mcont::value_type>::type(),
                "mcont::value_type must be a mutation type");

          public:
            virtual ~popbase() = default;
            popbase(popbase &&) = default;
            popbase(const popbase &) = default;
            popbase &operator=(popbase &&) = default;
            popbase &operator=(const popbase &) = default;
            //! Mutation type
            using mutation_t = mutation_type;
            //! Gamete type
            using haploid_genome_t = typename gcont::value_type;
            //! Mutation vec type
            using mcont_t = mcont;
            //! Mutation count vector type
            using mcount_t = std::vector<uint_t>;
            //! Gamete vec type
            using gcont_t = gcont;
            //! Lookup table type for recording mutation positions, etc.
            using lookup_table_t = lookup_table_type;
            //! container type for fixations
            using mvector_t = mvector;
            //! container type for fixation times
            using ftvector_t = ftvector;

            mcont_t mutations;
            /*!
              Used to keep track of mutation frequencies.

              Should have memory reserved externally,
              based on some good guess.
            */

            //! Contains number of times each mutation exists
            mcount_t mcounts;
            /// Contains the contribution of ancient samples
            /// to mutation counts.  This is place within
            /// population classes out of convienience.
            /// Be aware of the following issue regarding
            /// serialization: mcounts/mcounts_from_preserved_nodes
            /// are NOT serialized in fwdpp::serialize_population.
            /// Further, when pops are deserialized, mcounts is
            /// filled by a call to fwdpp_internal::process_haploid_genomes.
            /// This behavior is generally incorrect for simulations
            /// involving tree sequence recording, and a call to
            /// fwdpp::ts::count_mutations will be needed upon
            /// deserialization.
            /// \version 0.7.3 Added to library
            mcount_t mcounts_from_preserved_nodes;
            //! Container of haploid_genomes
            gcont haploid_genomes;

            /*!
              Containers that can be used as intermediates during
              the generation of new haploid_genomes.

              The requirement to declare these was introduced in fwdpp 0.3.3.

              In previous versions of the library, vectors like this had to be
              allocated
              for every crossover event for every generation.  The result was
              an excessive
              number of requests for memory allocation.

              Now, we create the vector once per replicate.  Further, we will
              reserve memory
              here, to minimize reallocs, etc., within fwdpp.

              Internally, fwdpp's job is to make sure that this vector is
              appropriately
              and efficiently cleared, but only when needed.

            */
            typename haploid_genome_t::mutation_container neutral, selected;

            /*!
              \brief Can be used to track positions of segregating mutations.
              \note Must have interface like std::map or std::unordered_set
            */
            lookup_table_type mut_lookup;
            //! Vector of mutation_t to track fixations
            mvector fixations;
            /*! \brief vector<uint_t> records times when mutation_ts
              were added to mut_lookup
            */
            ftvector fixation_times;

            //! Constructor
            explicit popbase(
                const uint_t &initial_haploid_genome_count,
                typename haploid_genome_t::mutation_container::size_type reserve_size
                = 100)
                : // No muts in the population
                  mutations(mcont_t()), mcounts(mcount_t()),
                  mcounts_from_preserved_nodes(mcount_t()),
                  haploid_genomes(
                      gcont(1, haploid_genome_t(initial_haploid_genome_count))),
                  neutral(typename haploid_genome_t::mutation_container()),
                  selected(typename haploid_genome_t::mutation_container()),
                  mut_lookup(lookup_table_type()), fixations(mvector()),
                  fixation_times(ftvector())
            {
                // This is a good number for reserving,
                // allowing for extra allocations when recycling is doing its
                // thing
                haploid_genomes.reserve(2 * initial_haploid_genome_count);
                // Reserve memory
                neutral.reserve(reserve_size);
                selected.reserve(reserve_size);
            }

            bool
            is_equal(const popbase &rhs) const
            {
                return this->mutations == rhs.mutations && this->mcounts == rhs.mcounts
                       && this->mcounts_from_preserved_nodes
                              == rhs.mcounts_from_preserved_nodes
                       && this->haploid_genomes == rhs.haploid_genomes
                       && this->fixations == rhs.fixations
                       && this->fixation_times == rhs.fixation_times;
            }

            //! Empty all the containers
            void
            clear_containers()
            {
                mutations.clear();
                mcounts.clear();
                mcounts_from_preserved_nodes.clear();
                haploid_genomes.clear();
                mut_lookup.clear();
                fixations.clear();
                fixation_times.clear();
            }
        };
    } // namespace poptypes
} // namespace fwdpp
#endif
