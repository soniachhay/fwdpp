#ifndef FWDPP_TS_SIMPLIFICATION_SIMPLIFICATION_HPP__
#define FWDPP_TS_SIMPLIFICATION_SIMPLIFICATION_HPP__

#include <limits>
#include <vector>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <fwdpp/ts/definitions.hpp>

namespace fwdpp
{
    namespace ts
    {
        namespace simplification
        {
            struct segment
            {
                double left, right;
                TS_NODE_INT node;
                segment(double l, double r, TS_NODE_INT n) : left{l}, right{r}, node{n}
                {
                    if (right <= left)
                        {
                            throw std::invalid_argument("right must be > left");
                        }
                }
            };

            struct mutation_node_map_entry
            {
                TS_NODE_INT node;
                std::size_t site, location;
                mutation_node_map_entry(TS_NODE_INT n, std::size_t s, std::size_t l)
                    : node(n), site(s), location(l)
                {
                }
            };

            class segment_overlapper
            /// This class is an iterable object
            /// over [left, right) -> segment
            /// mappings, where the segments
            /// are the genomic intervals in
            /// child nodes overlapping with
            /// the current parent in the
            /// current genomic interval.
            {
              private:
                std::vector<segment> segment_queue, overlapping;
                std::vector<segment>::const_iterator sbeg, send;
                std::vector<segment>::iterator overlapping_end;
                double _left, _right;

                inline double
                set_partition()
                {
                    double tright = std::numeric_limits<double>::max();
                    auto b = std::begin(overlapping);
                    for (auto i = std::begin(overlapping); i < overlapping_end; ++i)
                        {
                            if (i->right > _left)
                                {
                                    *b = *i;
                                    tright = std::min(tright, b->right);
                                    ++b;
                                }
                        }
                    overlapping_end = b;
                    return tright;
                }

              public:
                segment_overlapper()
                    : segment_queue{}, overlapping{}, sbeg(std::begin(segment_queue)),
                      send(std::end(segment_queue)),
                      overlapping_end(std::end(overlapping)), _left(0),
                      _right(std::numeric_limits<double>::max())
                {
                }

                void
                init()
                {
                    sbeg = std::begin(segment_queue);
                    // The - 1 for send assumes a "cap"/sentinel value.
                    send = std::end(segment_queue) - 1;
                    overlapping.clear();
                    overlapping_end = std::end(overlapping);
                    _left = 0.0;
                    _right = std::numeric_limits<double>::max();
                }

                bool
                operator()()
                {
                    bool rv = 0;
                    if (sbeg < send)
                        {
                            _left = _right;
                            auto tright = set_partition();
                            if (num_overlaps() == 0)
                                {
                                    _left = sbeg->left;
                                }
                            while (sbeg < send && sbeg->left == _left)
                                {
                                    tright = std::min(tright, sbeg->right);
                                    overlapping_end
                                        = overlapping.insert(overlapping_end, *sbeg) + 1;
                                    ++sbeg;
                                }
                            _right = std::min(sbeg->left, tright);
                            rv = true;
                        }
                    else
                        {
                            _left = _right;
                            _right = std::numeric_limits<double>::max();
                            auto tright = set_partition();
                            if (num_overlaps() > 0)
                                {
                                    _right = tright;
                                    rv = true;
                                }
                        }
                    return rv;
                }

                inline double
                left() const
                {
                    return _left;
                }

                inline double
                right() const
                {
                    return _right;
                }

                void
                clear_queue()
                {
                    segment_queue.clear();
                }

                template <typename... T>
                inline void
                enqueue(T&&... args)
                {
                    segment_queue.emplace_back(std::forward<T>(args)...);
                }

                void
                finalize_queue(double maxlen)
                {
                    std::sort(std::begin(segment_queue), std::end(segment_queue),
                              [](const segment& a, const segment& b) {
                                  return a.left < b.left;
                              });
                    // Add sentinel
                    segment_queue.emplace_back(maxlen, maxlen + 1.0, TS_NULL_NODE);
                }

                std::int64_t
                num_overlaps()
                {
                    return std::distance(std::begin(overlapping), overlapping_end);
                }

                inline const segment&
                overlap_front() const
                {
                    return overlapping.front();
                }

                std::vector<segment>::const_iterator
                begin() const
                {
                    return std::begin(overlapping);
                }

                std::vector<segment>::const_iterator
                end() const
                {
                    return overlapping_end;
                }
            };

            class ancestry_list
            {
              private:
                void
                resize_and_fill(std::vector<std::int32_t>& v, std::size_t n)
                {
                    v.resize(n);
                    std::fill(begin(v), end(v), -1);
                }

              public:
                std::vector<segment> segments;
                std::vector<std::int32_t> head, next;

                ancestry_list() : segments(), head(), next()
                {
                }

                void
                init(std::size_t n)
                {
                    segments.clear();
                    resize_and_fill(head, n);
                    resize_and_fill(next, n);
                }

                std::int32_t
                get_list_tail(std::size_t i) const
                {
                    if (i >= head.size())
                        {
                            throw std::runtime_error("index out of range");
                        }
                    auto f = head[i];
                    while (f != -1 && next[f] != -1)
                        {
                            f = next[f];
                        }
                    return f;
                }

                void
                add_record(std::size_t i, double l, double r, TS_NODE_INT n)
                {
                    if (i >= head.size())
                        {
                            throw std::runtime_error("index out of range");
                        }
                    segments.emplace_back(l, r, n);
                    if (head[i] == -1)
                        {
                            head[i] = segments.size() - 1;
                            if (segments.size() >= next.size())
                                {
                                    next.push_back(-1);
                                }
                        }
                    else
                        {
                            next.push_back(-1);
                            auto l = get_list_tail(i);
                            next[l] = segments.size() - 1;
                        }
                }

                void
                nullify_list(std::size_t i)
                {
                    if (i >= head.size())
                        {
                            throw std::runtime_error("index out of range");
                        }
                    head[i] = -1;
                }
            };

            template <typename TableCollectionType> struct simplifier_internal_state
            /// Holds data needed during tree sequence simplification
            /// \version Added in 0.9
            {
                using table_type = TableCollectionType;
                using edge_t = typename TableCollectionType::edge_t;
                using node_t = typename TableCollectionType::node_t;
                typename table_type::edge_table new_edge_table;
                typename table_type::edge_table temp_edge_buffer;
                typename table_type::node_table new_node_table;
                typename table_type::site_table new_site_table;
                ancestry_list ancestry;
                segment_overlapper overlapper;
                // NOTE: the whole idea of mutation map could
                // go away?  Should benchmark (later) with
                // high-mutation rate simulations.
                std::vector<mutation_node_map_entry> mutation_map;

                simplifier_internal_state()
                    : new_edge_table{}, temp_edge_buffer{}, new_node_table{},
                      new_site_table{}, ancestry{}, overlapper{}, mutation_map{}
                {
                }

                void
                clear()
                {
                    new_edge_table.clear();
                    new_node_table.clear();
                    temp_edge_buffer.clear();
                    new_site_table.clear();
                }
            };

            template <typename TableCollectionType>
            inline simplifier_internal_state<TableCollectionType>
            make_simplifier_internal_state(const TableCollectionType&)
            /// Convenience function to return a simplifier_internal_state
            {
                return simplifier_internal_state<TableCollectionType>();
            }

            template <typename SimplifierState>
            inline void
            buffer_edge(SimplifierState& state, const double left, const double right,
                        const TS_NODE_INT parent, const TS_NODE_INT child)
            {
                auto itr = std::find_if(
                    state.temp_edge_buffer.rbegin(), state.temp_edge_buffer.rend(),
                    [child](const typename SimplifierState::edge_t& e) {
                        return e.child == child;
                    });
                if (itr == state.temp_edge_buffer.rend())
                    {
                        state.temp_edge_buffer.emplace_back(
                            typename SimplifierState::edge_t{left, right, parent,
                                                             child});
                    }
                else
                    {
                        if (itr->right == left)
                            {
                                itr->right = right;
                            }
                        else
                            {
                                state.temp_edge_buffer.emplace_back(
                                    typename SimplifierState::edge_t{left, right, parent,
                                                                     child});
                            }
                    }
            }

            template <typename SimplifierState>
            inline std::size_t
            output_buffered_edges(SimplifierState& state)
            /// Take our buffered edges and add them to the output edge table
            {
                std::stable_sort(begin(state.temp_edge_buffer),
                                 end(state.temp_edge_buffer),
                                 [](const typename SimplifierState::edge_t& a,
                                    const typename SimplifierState::edge_t& b) {
                                     return a.child < b.child;
                                 });
                state.new_edge_table.insert(end(state.new_edge_table),
                                            begin(state.temp_edge_buffer),
                                            end(state.temp_edge_buffer));
                return state.temp_edge_buffer.size();
            }

            inline void
            add_ancestry(TS_NODE_INT input_id, double left, double right,
                         TS_NODE_INT node, ancestry_list& ancestry)
            {
                if (ancestry.head[input_id] == -1)
                    {
                        ancestry.add_record(input_id, left, right, node);
                    }
                else
                    {
                        auto last_idx = ancestry.get_list_tail(input_id);
                        if (last_idx == -1)
                            {
                                throw std::runtime_error("ancestry_list data invalid");
                            }
                        auto& last = ancestry.segments[last_idx];
                        if (last.right == left && last.node == node)
                            {
                                last.right = right;
                            }
                        else
                            {
                                ancestry.add_record(input_id, left, right, node);
                            }
                    }
            }

            template <typename TableCollectionType>
            inline void
            merge_ancestors(
                double maxlen,
                const typename TableCollectionType::node_table& input_node_table,
                const TS_NODE_INT parent_input_id,
                simplifier_internal_state<TableCollectionType>& state,
                std::vector<TS_NODE_INT>& idmap)
            {
                auto output_id = idmap[parent_input_id];
                bool is_sample = (output_id != TS_NULL_NODE);
                if (is_sample == true)
                    {
                        state.ancestry.nullify_list(parent_input_id);
                    }
                double previous_right = 0.0;
                state.overlapper.init();
                TS_NODE_INT ancestry_node = TS_NULL_NODE;
                state.temp_edge_buffer.clear();
                while (state.overlapper() == true)
                    {
                        if (state.overlapper.num_overlaps() == 1)
                            {
                                ancestry_node = state.overlapper.overlap_front().node;
                                if (is_sample)
                                    {
                                        buffer_edge(state, state.overlapper.left(),
                                                    state.overlapper.right(), output_id,
                                                    ancestry_node);
                                        ancestry_node = output_id;
                                    }
                            }
                        else
                            {
                                if (output_id == TS_NULL_NODE)
                                    {
                                        state.new_node_table.emplace_back(
                                            typename TableCollectionType::node_t{
                                                input_node_table[parent_input_id].deme,
                                                input_node_table[parent_input_id].time});
                                        output_id = state.new_node_table.size() - 1;
                                        // update sample map
                                        idmap[parent_input_id] = output_id;
                                    }
                                ancestry_node = output_id;
                                for (auto& x : state.overlapper)
                                    {
                                        buffer_edge(state, state.overlapper.left(),
                                                    state.overlapper.right(), output_id,
                                                    x.node);
                                    }
                            }
                        if (is_sample && state.overlapper.left() != previous_right)
                            {
                                add_ancestry(parent_input_id, previous_right,
                                             state.overlapper.left(), output_id,
                                             state.ancestry);
                            }
                        add_ancestry(parent_input_id, state.overlapper.left(),
                                     state.overlapper.right(), ancestry_node,
                                     state.ancestry);
                        previous_right = state.overlapper.right();
                    }
                if (is_sample && previous_right != maxlen)
                    {
                        add_ancestry(parent_input_id, previous_right, maxlen, output_id,
                                     state.ancestry);
                    }
                if (output_id != TS_NULL_NODE)
                    {
                        auto n = output_buffered_edges(state);
                        if (!n && !is_sample)
                            {
                                state.new_node_table.erase(begin(state.new_node_table)
                                                               + output_id,
                                                           end(state.new_node_table));
                                idmap[parent_input_id] = TS_NULL_NODE;
                            }
                    }
            }

            template <typename Iterator, typename SimplifierState>
            inline Iterator
            find_parent_child_segment_overlap(double maxlen, Iterator edge_ptr,
                                              const Iterator edge_end, TS_NODE_INT u,
                                              SimplifierState& state)
            {
                state.overlapper.clear_queue();
                for (; edge_ptr < edge_end && edge_ptr->parent == u; ++edge_ptr)
                    {
                        // For each edge corresponding to this parent,
                        // we look at all segments from the child.
                        // If the two segments overlap, we add the
                        // minimal overlap to our queue.
                        auto idx = state.ancestry.head[edge_ptr->child];
                        while (idx != -1)
                            {
                                auto& seg = state.ancestry.segments[idx];
                                if (seg.right > edge_ptr->left
                                    && edge_ptr->right > seg.left)
                                    {
                                        state.overlapper.enqueue(
                                            std::max(seg.left, edge_ptr->left),
                                            std::min(seg.right, edge_ptr->right),
                                            seg.node);
                                    }
                                idx = state.ancestry.next[idx];
                            }
                    }
                state.overlapper.finalize_queue(maxlen);
                return edge_ptr;
            }

            template <typename TableCollectionType>
            inline void
            record_sample_nodes(const std::vector<TS_NODE_INT>& samples,
                                const TableCollectionType& tables,
                                simplifier_internal_state<TableCollectionType>& state,
                                std::vector<TS_NODE_INT>& idmap)
            /// \version 0.7.1 Throw exception if a sample is recorded twice
            {
                for (const auto& s : samples)
                    {
                        // See GitHub issue 158
                        // for background
                        if (idmap[s] != TS_NULL_NODE)
                            {
                                throw std::invalid_argument("invalid sample list");
                            }
                        state.new_node_table.emplace_back(
                            typename TableCollectionType::node_t{tables.nodes[s].deme,
                                                                 tables.nodes[s].time});
                        add_ancestry(
                            s, 0, tables.genome_length(),
                            static_cast<TS_NODE_INT>(state.new_node_table.size() - 1),
                            state.ancestry);
                        idmap[s]
                            = static_cast<TS_NODE_INT>(state.new_node_table.size() - 1);
                    }
            }

            template <typename SiteTable, typename Mutation>
            inline void
            record_site(const SiteTable& sites, SiteTable& new_site_table, Mutation& mr)
            {
                double pos = sites[mr.site].position;
                if (new_site_table.empty() || new_site_table.back().position != pos)
                    {
                        new_site_table.push_back(sites[mr.site]);
                    }
                mr.site = new_site_table.size() - 1;
            }

            template <typename TableCollectionType>
            inline void
            prep_mutation_simplification(
                const TableCollectionType& input_tables,
                std::vector<mutation_node_map_entry>& mutation_map)
            {
                mutation_map.clear();
                mutation_map.reserve(input_tables.mutations.size());
                for (std::size_t i = 0; i < input_tables.mutations.size(); ++i)
                    {
                        mutation_map.emplace_back(input_tables.mutations[i].node,
                                                  input_tables.mutations[i].site, i);
                    }

                std::sort(
                    begin(mutation_map), end(mutation_map),
                    [&input_tables](const mutation_node_map_entry& a,
                                    const mutation_node_map_entry& b) {
                        return std::tie(a.node, input_tables.sites[a.site].position)
                               < std::tie(b.node, input_tables.sites[b.site].position);
                    });
            }

            template <typename TableCollectionType, typename PreservedVariantIndexes>
            inline void
            simplify_mutations(simplifier_internal_state<TableCollectionType>& state,
                               TableCollectionType& input_tables,
                               PreservedVariantIndexes& preserved_variants)
            // Remove all mutations that do not map to nodes
            // in the simplified tree.  The key here is
            // that ancestry contains the history of
            // each node, which we use for the remapping.
            {
                static_assert(
                    std::is_integral<
                        typename PreservedVariantIndexes::value_type>::value,
                    "PreservedVariantIndexes::value_type must be an integer type");
                prep_mutation_simplification(input_tables, state.mutation_map);
                // Set all output nodes to null for now.
                for (auto& mr : input_tables.mutations)
                    {
                        mr.node = TS_NULL_NODE;
                    }

                // Map the input node id of a mutation to
                // its output node id.  If no output ID exists,
                // then the mutation will be removed by the
                // call to erase below.
                auto map_itr = begin(state.mutation_map);
                const auto map_end = end(state.mutation_map);

                while (map_itr < map_end)
                    {
                        auto n = map_itr->node;
                        auto seg_idx = state.ancestry.head[n];
                        for (; map_itr < map_end && map_itr->node == n;)
                            {
                                if (seg_idx == -1)
                                    {
                                        ++map_itr;
                                        break;
                                    }
                                while (seg_idx != -1 && map_itr < map_end
                                       && map_itr->node == n)
                                    {
                                        auto& seg = state.ancestry.segments[seg_idx];
                                        auto pos
                                            = input_tables.sites[map_itr->site].position;
                                        if (seg.left <= pos && pos < seg.right)
                                            {
                                                input_tables.mutations[map_itr->location]
                                                    .node
                                                    = seg.node;
                                                ++map_itr;
                                            }
                                        else if (pos >= seg.right)
                                            {
                                                seg_idx = state.ancestry.next[seg_idx];
                                            }
                                        else
                                            {
                                                ++map_itr;
                                            }
                                    }
                            }
                    }

                // Any mutations with null node values do not have
                // ancestry and may be removed.
                auto itr = std::remove_if(
                    begin(input_tables.mutations), end(input_tables.mutations),
                    [](const typename TableCollectionType::mutation_t& mr) {
                        return mr.node == TS_NULL_NODE;
                    });
                preserved_variants.clear();
                preserved_variants.reserve(
                    std::distance(itr, input_tables.mutations.end()));
                for (auto i = begin(input_tables.mutations); i != itr; ++i)
                    {
                        record_site(input_tables.sites, state.new_site_table, *i);
                        preserved_variants.push_back(i->key);
                    }

                input_tables.mutations.erase(itr, input_tables.mutations.end());
                input_tables.sites.swap(state.new_site_table);
                //TODO: replace assert with exception
                assert(std::is_sorted(
                    begin(input_tables.mutations), end(input_tables.mutations),
                    [&input_tables](const typename TableCollectionType::mutation_t& a,
                                    const typename TableCollectionType::mutation_t& b) {
                        return input_tables.sites[a.site].position
                               < input_tables.sites[b.site].position;
                    }));
            }

            template <typename Iterator, typename TableCollectionType>
            inline void
            transfer_new_nodes_and_edges(
                const Iterator new_edge_destination,
                simplifier_internal_state<TableCollectionType>& state,
                TableCollectionType& tables)
            /// Update the tables.  To keep memory use as sane as possible,
            /// we use resize-and-move here.  In theory, we can also do
            /// vector swaps, but that has a side-effect of keeping
            /// far too much RAM allocated compared to what we need.
            {
                tables.edges.resize(
                    std::distance(begin(tables.edges), new_edge_destination));
                tables.nodes.resize(state.new_node_table.size());
                std::move(begin(state.new_node_table), end(state.new_node_table),
                          begin(tables.nodes));
                // TODO: allow for exception instead of assert
                assert(tables.edges_are_minimally_sorted());
                // NOTE: this will be moot by removing edge sorting.
                tables.update_offset();
            }

            template <typename Iterator, typename ConstIterator,
                      typename SimplifierState>
            inline Iterator
            transfer_simplified_edges(Iterator new_edge_destination,
                                      ConstIterator input_edge_table_location,
                                      std::size_t minsize, SimplifierState& state)
            {
                if (state.new_edge_table.size() >= minsize
                    && new_edge_destination + state.new_edge_table.size()
                           < input_edge_table_location)
                    {
                        auto rv
                            = std::copy(begin(state.new_edge_table),
                                        end(state.new_edge_table), new_edge_destination);
                        state.new_edge_table.clear();
                        return rv;
                    }
                return new_edge_destination;
            }
        }

    }
}

#endif
