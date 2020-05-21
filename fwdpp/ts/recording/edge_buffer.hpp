#ifndef FWDPP_TS_RECORDING_EDGE_BUFFER_HPP
#define FWDPP_TS_RECORDING_EDGE_BUFFER_HPP

#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <limits>
#include <cstdint>
#include <vector>
#include "../definitions.hpp"

namespace fwdpp
{
    namespace ts
    {
        constexpr std::int64_t EDGE_BUFFER_NULL = -1;

        struct birth_data
        {
            double left, right;
            TS_NODE_INT child;
            std::int64_t next;

            birth_data(double l, double r, TS_NODE_INT c)
                : left{l}, right{r}, child{c}, next{EDGE_BUFFER_NULL}
            {
            }
        };

        struct edge_buffer
        {
            std::vector<std::int64_t> head;
            std::vector<birth_data> births;

            edge_buffer() : head{}, births{}
            {
            }
        };

        struct parent_location
        {
            TS_NODE_INT parent;
            std::size_t start, stop;
            parent_location(TS_NODE_INT p, std::size_t start_, std::size_t stop_)
                : parent{p}, start{start_}, stop{stop_}
            {
            }
        };

        inline std::int64_t
        get_buffer_end(const edge_buffer& new_edges, std::size_t i)
        {
            if (i >= new_edges.head.size())
                {
                    throw std::runtime_error("invalid parent index");
                }
            auto f = new_edges.head[i];
            while (f != EDGE_BUFFER_NULL && new_edges.births[f].next != EDGE_BUFFER_NULL)
                {
                    f = new_edges.births[f].next;
                    if (f != EDGE_BUFFER_NULL && f >= new_edges.births.size())
                        {
                            throw std::runtime_error("invalid next value");
                        }
                }
            return f;
        }

        inline std::int64_t
        buffer_new_edge(TS_NODE_INT parent, double left, double right, TS_NODE_INT child,
                        edge_buffer& new_edges)
        {
            if (parent == TS_NULL_NODE || child == TS_NULL_NODE)
                {
                    throw std::runtime_error("bad node IDs passed to buffer_new_edge");
                }
            if (parent >= new_edges.head.size())
                {
                    new_edges.head.resize(parent + 1, EDGE_BUFFER_NULL);
                }
            new_edges.births.emplace_back(left, right, child);
            if (new_edges.head[parent] == EDGE_BUFFER_NULL)
                {
                    new_edges.head[parent] = new_edges.births.size() - 1;
                    if (new_edges.births[new_edges.head[parent]].next
                        != EDGE_BUFFER_NULL)
                        {
                            std::ostringstream o;
                            o << "invalid next entry for head birth: " << parent << ' '
                              << new_edges.head[parent] << ' '
                              << new_edges.births[new_edges.head[parent]].next << ' '
                              << new_edges.head.size() << ' ' << new_edges.births.size();
                            throw std::runtime_error(o.str());
                        }
                }
            else
                {
                    auto l = get_buffer_end(new_edges, parent);
                    new_edges.births[l].next = new_edges.births.size() - 1;
                }
            return new_edges.births.size() - 1;
        }

        inline std::int64_t
        buffer_new_edge_at(std::int64_t loc, double left, double right,
                           TS_NODE_INT child, edge_buffer& new_edges)
        {
            if (loc >= new_edges.births.size())
                {
                    throw std::runtime_error("bad location");
                }
            new_edges.births.emplace_back(left, right, child);
            new_edges.births[loc].next = new_edges.births.size() - 1;
            return new_edges.births.size() - 1;
        }

        // Below are functions for liftover of an edge buffer
        // to a table collection

        template <typename TableCollectionType>
        inline std::vector<parent_location>
        find_pre_existing_edges(
            const TableCollectionType& tables,
            const std::vector<TS_NODE_INT>& alive_at_last_simplification,
            const fwdpp::ts::edge_buffer& new_edges)
        // FIXME: the indexing step need go no farther than the time of the most
        // recent node in alive_at_last_simplification.
        {
            std::vector<TS_NODE_INT> alive_with_new_edges;
            for (auto a : alive_at_last_simplification)
                {
                    if (new_edges.head[a] != EDGE_BUFFER_NULL)
                        {
                            alive_with_new_edges.push_back(a);
                        }
                }
            if (alive_with_new_edges.empty()) // get out early
                {
                    return {};
                }

            // index where each node already has edges.
            std::vector<std::size_t> starts(tables.num_nodes(),
                                            std::numeric_limits<std::size_t>::max()),
                stops(tables.num_nodes(), std::numeric_limits<std::size_t>::max());
            for (std::size_t i = 0; i < tables.num_edges(); ++i)
                {
                    if (starts[tables.edges[i].parent]
                        == std::numeric_limits<std::size_t>::max())
                        {
                            starts[tables.edges[i].parent] = i;
                            stops[tables.edges[i].parent]
                                = i; // FIXME: idiomatically, this should be i+1
                        }
                    else
                        {
                            stops[tables.edges[i].parent]
                                = i; // FIXME: idiomatically, this should be i+1
                        }
                }

            std::vector<parent_location> existing_edges;
            for (auto a : alive_with_new_edges)
                {
                    existing_edges.emplace_back(a, starts[a], stops[a]);
                }

            // Our only sort!!
            std::sort(begin(existing_edges), end(existing_edges),
                      [&tables](const parent_location& lhs, const parent_location& rhs) {
                          // lexical comparison of tuple elements just like in Python
                          return std::tie(tables.nodes[lhs.parent].time, lhs.start,
                                          lhs.parent)
                                 < std::tie(tables.nodes[rhs.parent].time, rhs.start,
                                            rhs.parent);
                      });

            // FIXME: this should be debug only
            for (std::size_t i = 1; i < existing_edges.size(); ++i)
                {
                    auto t0 = tables.nodes[existing_edges[i - 1].parent].time;
                    auto t1 = tables.nodes[existing_edges[i].parent].time;
                    if (t0 > t1)
                        {
                            throw std::runtime_error(
                                "existing edges not properly sorted by time");
                        }
                }

            return existing_edges;
        }

        template <typename TableCollectionType>
        std::size_t
        handle_pre_existing_edges(
            const TableCollectionType& tables, const edge_buffer& new_edges,
            const std::vector<parent_location>& existing_edges,
            typename TableCollectionType::edge_table& edge_liftover)
        {
            std::size_t offset = 0;
            for (const auto& ex : existing_edges)
                {
                    // FIXME: this while loop is repeated 2x just w/different
                    // ranges
                    while (offset < tables.num_edges()
                           && tables.nodes[tables.edges[offset].parent].time
                                  < tables.nodes[ex.parent].time)
                        {
                            edge_liftover.emplace_back(tables.edges[offset]);
                            //edge_liftover.add_edge(tables.edges.left[offset],
                            //                       tables.edges.right[offset],
                            //                       tables.edges.parent[offset],
                            //                       tables.edges.child[offset]);
                            ++offset;
                        }
                    if (ex.start != std::numeric_limits<std::size_t>::max())
                        {
                            while (offset < ex.start
                                   && tables.nodes[tables.edges[offset].parent].time
                                          <= tables.nodes[ex.parent].time)
                                {
                                    edge_liftover.emplace_back(tables.edges[offset]);
                                    //edge_liftover.add_edge(tables.edges[offset].left,
                                    //                       tables.edges[offset].right,
                                    //                       tables.edges[offset].parent,
                                    //                       tables.edges[offset].child);
                                    ++offset;
                                }
                            // FIXME: stop condition isn't idiomatic
                            for (decltype(ex.start) i = ex.start; i < ex.stop + 1; ++i)
                                {
                                    edge_liftover.emplace_back(tables.edges[i]);
                                    //edge_liftover.add_edge(
                                    //    tables.edges.left[i], tables.edges.right[i],
                                    //    tables.edges.parent[i], tables.edges.child[i]);
                                }
                            offset = ex.stop + 1;
                        }
                    auto n = new_edges.head[ex.parent];
                    while (n != EDGE_BUFFER_NULL)
                        {
                            edge_liftover.emplace_back(
                                typename TableCollectionType::edge_t{
                                    new_edges.births[n].left, new_edges.births[n].right,
                                    ex.parent, new_edges.births[n].child});
                            n = new_edges.births[n].next;
                        }
                }
            return offset;
        }

        template <typename TableCollectionType>
        inline void
        copy_births_since_last_simplification(
            const edge_buffer& new_edges, const TableCollectionType& tables,
            double max_time, typename TableCollectionType::edge_table& edge_liftover)
        {

            // TODO: should validate data in new_edges
            edge_liftover.clear(); // Re-use this buffer b/c it'll get big.

            // Go backwards through new births, and add them
            // to our temporary edge table if they are newer
            // than the last simplification time

            for (auto b = new_edges.head.rbegin(); b < new_edges.head.rend(); ++b)
                {
                    auto d = std::distance(new_edges.head.rbegin(), b);
                    auto parent = new_edges.head.size() - d - 1;
                    auto ptime = tables.nodes[parent].time;
                    if (*b != EDGE_BUFFER_NULL && ptime < max_time)
                        {
                            auto n = *b;
                            while (n != EDGE_BUFFER_NULL)
                                {
                                    edge_liftover.emplace_back(
                                        typename TableCollectionType::edge_t{
                                            new_edges.births[n].left,
                                            new_edges.births[n].right, parent,
                                            new_edges.births[n].child});
                                    //edge_liftover.add_edge(new_edges.births[n].left,
                                    //                       new_edges.births[n].right,
                                    //                       parent,
                                    //                       new_edges.births[n].child);
                                    n = new_edges.births[n].next;
                                }
                        }
                    else if (*b != EDGE_BUFFER_NULL && ptime >= max_time)
                        {
                            break;
                        }
                }
        }

        template <typename TableCollectionType>
        void
        stitch_together_edges(
            const std::vector<TS_NODE_INT>& alive_at_last_simplification,
            double max_time, edge_buffer& new_edges,
            typename TableCollectionType::edge_table& edge_liftover,
            TableCollectionType& tables)
        {
            copy_births_since_last_simplification(new_edges, tables, max_time,
                                                  edge_liftover);
            auto existing_edges = find_pre_existing_edges(
                tables, alive_at_last_simplification, new_edges);
            auto offset = handle_pre_existing_edges(tables, new_edges, existing_edges,
                                                    edge_liftover);
            for (; offset < tables.num_edges(); ++offset)
                {
                    edge_liftover.emplace_back(tables.edges[offset]);
                    //tables.edges.left[offset], tables.edges.right[offset],
                    //tables.edges.parent[offset], tables.edges.child[offset]);
                }
            tables.edges.assign(begin(edge_liftover), end(edge_liftover));
            // This resets sizes to 0, but keeps the memory allocated.
            edge_liftover.clear();
            // TODO: move this cleanup to function
            new_edges.head.resize(tables.num_nodes());
            std::fill(begin(new_edges.head), end(new_edges.head), EDGE_BUFFER_NULL);
            new_edges.births.clear();
        }
    }
}

#endif
