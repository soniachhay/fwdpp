#ifndef FWDPP_TESTSUITE_EMPTY_TABLE_COLLECTION_HPP
#define FWDPP_TESTSUITE_EMPTY_TABLE_COLLECTION_HPP

#include <fwdpp/ts/table_collection.hpp>

struct empty_table_collection
{
    fwdpp::ts::table_collection tables;

    empty_table_collection() : tables(1.) {}
};

#endif
