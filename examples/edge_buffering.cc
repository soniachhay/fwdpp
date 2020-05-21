#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <gsl/gsl_randist.h>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpp/ts/recording.hpp>

namespace po = boost::program_options;

struct command_line_options
{
    unsigned N;
    double psurvival;
    unsigned nsteps;
    unsigned simplification_interval;
    double rho;
    std::string treefile;
    bool buffer_new_edges;
    unsigned seed;

    command_line_options();
};

command_line_options::command_line_options()
    : N{1000}, psurvival{0.}, nsteps{1000}, simplification_interval{100}, rho{0.},
      treefile{"treefile.trees"}, buffer_new_edges{false}, seed{42}
{
}

void
validate_cli(const command_line_options& options)
{
    if (options.N == 0)
        {
            throw std::invalid_argument("Population size must be > 0");
        }

    if (options.psurvival < 0. || options.psurvival >= 1.
        || std::isfinite(options.psurvival) == false)
        {
            throw std::invalid_argument("psurvival must be 0.0 <= p < 1.0");
        }

    if (options.rho < 0.0 || std::isfinite(options.rho) == false)
        {
            throw std::invalid_argument("rho must be >= 0.0");
        }

    if (options.treefile.empty())
        {
            throw std::invalid_argument("treefile must not be an empty string");
        }
}

po::options_description
generate_main_options(command_line_options& o)
{
    po::options_description options("Simulation options");
    options.add_options()("help", "Display help");
    options.add_options()("N", po::value<decltype(command_line_options::N)>(&o.N),
                          "Diploid population size. Default = 1000.");

    options.add_options()(
        "psurvival", po::value<decltype(command_line_options::psurvival)>(&o.psurvival),
        "Survival probability. Default = 0.0");
    options.add_options()("nsteps",
                          po::value<decltype(command_line_options::nsteps)>(&o.nsteps),
                          "Number of time steps to evolve. Default = 1000.");
    options.add_options()(
        "simplify",
        po::value<decltype(command_line_options::simplification_interval)>(
            &o.simplification_interval),
        "Time steps between simplifications.  Default = 100.");
    options.add_options()("rho", po::value<decltype(command_line_options::rho)>(&o.rho),
                          "Scaled recombination rate, 4Nr.  Default=0.");
    options.add_options()(
        "treefile", po::value<decltype(command_line_options::treefile)>(&o.treefile),
        "Ouput file name.  Default = treefile.trees");
    options.add_options()(
        "buffer", po::bool_switch(&o.buffer_new_edges),
        "If true, use edge buffering algorithm. If not, sort and simplify. Default = "
        "false");
    options.add_options()("seed",
                          po::value<decltype(command_line_options::seed)>(&o.seed),
                          "Random number seed.  Default = 42.");

    return options;
}

// Simulation types and functions now

struct parent
{
    std::size_t index;
    fwdpp::ts::TS_NODE_INT node0, node1;
    parent(std::size_t i, fwdpp::ts::TS_NODE_INT n0, fwdpp::ts::TS_NODE_INT n1)
        : index(i), node0(n0), node1(n1)
    {
    }
};

struct birth
{
    std::size_t index;
    fwdpp::ts::TS_NODE_INT p0node0, p0node1, p1node0, p1node1;
    birth(std::size_t i, const parent& p0, const parent& p1)
        : index(i), p0node0(p0.node0), p0node1(p0.node1), p1node0(p1.node0),
          p1node1(p1.node1)
    {
    }
};

void
deaths_and_parents(const fwdpp::GSLrng_mt& rng, const std::vector<parent>& parents,
                   double psurvival, std::vector<birth>& births)
{
    births.clear();
    for (std::size_t i = 0; i < parents.size(); ++i)
        {
            if (gsl_rng_uniform(rng.get()) > psurvival)
                {
                    std::size_t parent0 = gsl_ran_flat(rng.get(), 0, parents.size());
                    std::size_t parent1 = gsl_ran_flat(rng.get(), 0, parents.size());
                    births.emplace_back(i, parents[parent0], parents[parent1]);
                }
        }
}

void
recombination_breakpoints(const fwdpp::GSLrng_mt& rng, double littler, double maxlen,
                          std::vector<double>& breakpoints)
{
    breakpoints.clear();
    auto nxovers = gsl_ran_poisson(rng.get(), littler);
    for (decltype(nxovers) i = 0; i < nxovers; ++i)
        {
            breakpoints.push_back(gsl_ran_flat(rng.get(), 0., maxlen));
        }
    std::sort(begin(breakpoints), end(breakpoints));

    // Remove all values that do not exist an odd number of times
    auto itr = std::adjacent_find(begin(breakpoints), end(breakpoints));
    if (itr != end(breakpoints))
        {
            std::vector<double> temp;
            auto start = begin(breakpoints);
            while (itr < end(breakpoints))
                {
                    auto not_equal
                        = std::find_if(itr, breakpoints.end(),
                                       [itr](const double d) { return d != *itr; });
                    int even = (std::distance(itr, not_equal) % 2 == 0.0);
                    temp.insert(temp.end(), start, itr + 1 - even);
                    start = not_equal;
                    itr = std::adjacent_find(start, std::end(breakpoints));
                }
            temp.insert(end(temp), start, breakpoints.end());
            breakpoints.swap(temp);
        }
}

void
recombine_and_record_edges(const fwdpp::GSLrng_mt& rng, double littler,
                           std::vector<double>& breakpoints,
                           fwdpp::ts::TS_NODE_INT parental_node0,
                           fwdpp::ts::TS_NODE_INT parental_node1,
                           fwdpp::ts::TS_NODE_INT child,
                           fwdpp::ts::std_table_collection& tables)
// NOTE: this is an improvement on what I do in fwdpp?
{
    recombination_breakpoints(rng, littler, tables.genome_length(), breakpoints);
    double left = 0.;
    std::size_t breakpoint = 1;
    auto pnode0 = parental_node0;
    auto pnode1 = parental_node1;
    for (; breakpoint < breakpoints.size(); ++breakpoint)
        {
            auto rv
                = tables.emplace_back_edge(left, breakpoints[breakpoint], pnode0, child);
            std::swap(pnode0, pnode1);
            left = breakpoints[breakpoint];
        }
    auto rv = tables.emplace_back_edge(left, tables.genome_length(), pnode0, child);
}

void
recombine_and_buffer_edges(const fwdpp::GSLrng_mt& rng, double littler,
                           std::vector<double>& breakpoints,
                           fwdpp::ts::TS_NODE_INT parental_node0,
                           fwdpp::ts::TS_NODE_INT parental_node1,
                           fwdpp::ts::TS_NODE_INT child, double maxlen,
                           fwdpp::ts::edge_buffer& new_edges)
{
    recombination_breakpoints(rng, littler, maxlen, breakpoints);
    double left = 0.;
    std::size_t breakpoint = 1;
    auto pnode0 = parental_node0;
    auto pnode1 = parental_node1;
    auto end = fwdpp::ts::EDGE_BUFFER_NULL, other_end = fwdpp::ts::EDGE_BUFFER_NULL;
    if (pnode0 < new_edges.head.size())
        {
            end = get_buffer_end(new_edges, pnode0);
        }
    if (pnode1 < new_edges.head.size())
        {
            other_end = get_buffer_end(new_edges, pnode1);
        }

    for (; breakpoint < breakpoints.size(); ++breakpoint)
        {
            if (end == -1)
                {
                    end = buffer_new_edge(pnode0, left, breakpoints[breakpoint], child,
                                          new_edges);
                }
            else
                {
                    end = buffer_new_edge_at(end, left, breakpoints[breakpoint], child,
                                             new_edges);
                }
            std::swap(pnode0, pnode1);
            std::swap(end, other_end);
            left = breakpoints[breakpoint];
        }
    if (end == -1)
        {
            end = buffer_new_edge(pnode0, left, maxlen, child, new_edges);
        }
    else
        {
            end = buffer_new_edge_at(end, left, maxlen, child, new_edges);
        }
}

static void
generate_births(const fwdpp::GSLrng_mt& rng, const std::vector<birth>& births,
                double littler, std::vector<double>& breakpoints, double birth_time,
                bool buffer_new_edges, fwdpp::ts::edge_buffer& new_edges,
                std::vector<parent>& parents, fwdpp::ts::std_table_collection& tables)
{
    for (auto& b : births)
        {
            auto new_node_0 = tables.emplace_back_node(0, birth_time);
            auto new_node_1 = tables.emplace_back_node(0, birth_time);
            auto p0n0 = b.p0node0;
            auto p0n1 = b.p0node1;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    std::swap(p0n0, p0n1);
                }
            auto p1n0 = b.p1node0;
            auto p1n1 = b.p1node1;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    std::swap(p1n0, p1n1);
                }
            if (buffer_new_edges == false)
                {
                    recombine_and_record_edges(rng, littler, breakpoints, p0n0, p0n1,
                                               new_node_0, tables);
                    recombine_and_record_edges(rng, littler, breakpoints, p1n0, p1n1,
                                               new_node_1, tables);
                }
            else
                {
                    double ptime = tables.nodes[p0n0].time;
                    double ctime = tables.nodes[new_node_0].time;
                    if (ctime >= ptime)
                        {
                            throw std::runtime_error("bad parent/child time");
                        }
                    recombine_and_buffer_edges(rng, littler, breakpoints, p0n0, p0n1,
                                               new_node_0, tables.genome_length(),
                                               new_edges);
                    ptime = tables.nodes[p1n0].time;
                    ctime = tables.nodes[new_node_1].time;
                    if (ctime >= ptime)
                        {
                            throw std::runtime_error("bad parent/child time");
                        }
                    recombine_and_buffer_edges(rng, littler, breakpoints, p1n0, p1n1,
                                               new_node_1, tables.genome_length(),
                                               new_edges);
                }
            parents[b.index] = parent(b.index, new_node_0, new_node_1);
        }
}
//
// NOTE: seems like samples could/should be const?
static void
sort_n_simplify(double last_time_simplified,
                std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                std::vector<fwdpp::ts::TS_NODE_INT>& node_map,
                fwdpp::ts::std_table_collection& tables)
{
}

static void
flush_buffer_n_simplify(
    std::vector<fwdpp::ts::TS_NODE_INT>& alive_at_last_simplification,
    std::vector<fwdpp::ts::TS_NODE_INT>& samples,
    std::vector<fwdpp::ts::TS_NODE_INT>& node_map, fwdpp::ts::edge_buffer& new_edges,
    fwdpp::ts::std_table_collection::edge_table& edge_liftover,
    fwdpp::ts::std_table_collection& tables)
{
    double max_time = std::numeric_limits<double>::max();
    for (auto a : alive_at_last_simplification)
        {
            max_time = std::min(max_time, tables.nodes[a].time);
        }

    stitch_together_edges(alive_at_last_simplification, max_time, new_edges,
                          edge_liftover, tables);
    //FIXME: simplify
}

void
simulate(const command_line_options& options)
{
    fwdpp::GSLrng_mt rng(options.seed);
    fwdpp::ts::std_table_collection tables(1.0);
    fwdpp::ts::edge_buffer buffer;
    fwdpp::ts::std_table_collection::edge_table edge_liftover;
    std::vector<parent> parents;
    for (unsigned i = 0; i < options.N; ++i)
        {
            auto id0 = tables.emplace_back_node(0, 0.);
            auto id1 = tables.emplace_back_node(0, 0.);
            parents.emplace_back(i, id0, id1);
        }

    // The next bits are all for buffering
    std::vector<fwdpp::ts::TS_NODE_INT> alive_at_last_simplification;

    if (options.buffer_new_edges)
        {
            if (buffer.head.size() != 2 * options.N)
                {
                    throw std::runtime_error("bad setup of edge_buffer_ptr");
                }
        }

    std::vector<birth> births;
    std::vector<fwdpp::ts::TS_NODE_INT> samples, node_map;
    bool simplified = false;
    double last_time_simplified = options.nsteps;
    double littler = options.rho / (2. * static_cast<double>(options.N));
    std::vector<double> breakpoints;
    for (unsigned step = 1; step <= options.nsteps; ++step)
        {
            deaths_and_parents(rng, parents, options.psurvival, births);
            generate_births(rng, births, littler, breakpoints, options.nsteps - step,
                            options.buffer_new_edges, buffer, parents, tables);
            if (step % options.simplification_interval == 0.)
                {
                    samples.clear();
                    for (auto& p : parents)
                        {
                            samples.push_back(p.node0);
                            samples.push_back(p.node1);
                        }
                    node_map.resize(tables.num_nodes());

                    if (options.buffer_new_edges == false)
                        {
                            sort_n_simplify(last_time_simplified, samples, node_map,
                                            tables);
                        }
                    else
                        {
                            flush_buffer_n_simplify(alive_at_last_simplification,
                                                    samples, node_map, buffer,
                                                    edge_liftover, tables);
                        }
                    simplified = true;
                    last_time_simplified = options.nsteps - step;
                    //remap parent nodes
                    for (auto& p : parents)
                        {
                            p.node0 = node_map[p.node0];
                            p.node1 = node_map[p.node1];
                        }
                    if (options.buffer_new_edges == true)
                        {
                            alive_at_last_simplification.clear();
                            for (auto& p : parents)
                                {
                                    alive_at_last_simplification.push_back(p.node0);
                                    alive_at_last_simplification.push_back(p.node1);
                                }
                        }
                }
            else
                {
                    simplified = false;
                }
        }
    if (simplified == false)
        {
            samples.clear();
            for (auto& p : parents)
                {
                    samples.push_back(p.node0);
                    samples.push_back(p.node1);
                }
            node_map.resize(tables.num_nodes());
            if (options.buffer_new_edges == false)
                {
                    sort_n_simplify(last_time_simplified, samples, node_map, tables);
                }
            else
                {
                    flush_buffer_n_simplify(alive_at_last_simplification, samples,
                                            node_map, buffer, edge_liftover, tables);
                }
        }
}

int
main(int argc, char** argv)
{
    command_line_options options;
    auto cli = generate_main_options(options);
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cli), vm);
    po::notify(vm);
    validate_cli(options);

    if (vm.count("help"))
        {
            std::cout << cli << '\n';
            std::exit(1);
        }

    simulate(options);
}
