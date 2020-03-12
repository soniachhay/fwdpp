// Example of reading in a table_collection that has
// been dumped to an uncompressed file.  For example,
// the output from wfts_integration_test

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <fwdpp/ts/serialization.hpp>

int
main(int argc, char** argv)
{
    if (argc < 2)
        {
            std::cerr << "Usage: " << argv[0] << " infile\n";
            std::exit(1);
        }

    std::ifstream in(argv[1]);
    if (!in)
        {
            throw std::invalid_argument("cannot open file for reading");
        }

    auto tables = fwdpp::ts::io::deserialize_tables(in);

    std::cout << "Nodes:\n";
    for (auto&& n : tables.node_table)
        {
            std::cout << n.time << ' ' << n.deme << '\n';
        }

    std::cout << "Edges:\n";
    for (auto&& e : tables.edge_table)
        {
            std::cout << e.left << ' ' << e.right << ' ' << e.parent << ' '
                      << e.child << '\n';
        }

    std::cout << "Sites:\n";
    for (auto&& s : tables.site_table)
        {
            std::cout << s.position << ' ' << s.ancestral_state << '\n';
        }

    std::cout << "Mutations:\n";
    for (auto&& m : tables.mutation_table)
        {
            std::cout << m.key << ' ' << m.node << ' ' << m.site << ' '
                      << m.neutral << ' ' << m.derived_state << '\n';
        }
}
