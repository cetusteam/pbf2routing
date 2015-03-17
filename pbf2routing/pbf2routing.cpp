/*
 This file is part of the osmpbf library.
 
 Copyright(c) 2014 Daniel Bahrdt.
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <osmpbf/osmfile.h>
#include <osmpbf/inode.h>
#include <osmpbf/iway.h>
#include <osmpbf/irelation.h>
#include <osmpbf/filter.h>
#include <osmpbf/primitiveblockinputadaptor.h>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <limits>
#include <sqlite3.h>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry.hpp>

namespace bg = boost::geometry;
typedef bg::cs::spherical_equatorial<bg::degree> coord_system_t;
typedef bg::model::point<double, 2, coord_system_t> degree_point;
using namespace std;

const double m = 0.000000001;
const double granularity = 1000;
const double m_by_g = m * granularity;

// Converts the cost from a double to an integer.
static uint32_t encode_cost(double cost) {
    return cost * 10.;
}

class StatusUpdater {
    uint32_t count, last_count;
    uint32_t total;
    clock_t start_time;
    clock_t last_time;
    string msg;
    string eol;
public:
    StatusUpdater() {
        if (isatty(fileno(stdout))) {
            // the output is a terminal. use interactive mode
            eol = "           \r";
        } else {
            // output is a file or pile. use simple mode.
            eol = "\n";
        }
    }
    void start(const string& msg, uint32_t total) {
        this->msg = msg;
        this->total = total;
        count = last_count = 0;
        start_time = last_time = clock();
        cerr << msg << ": starting" << eol << flush;
    }

    void absolute(uint32_t value) {
        count = value;
        print_message();
    }

    void progress(uint32_t increment=1) {
        count += increment;
        print_message();
    }

    void print_message() {
        clock_t end = clock();
        if ((end - last_time) < CLOCKS_PER_SEC) {
            return;
        }
        double tot_time = double(end - last_time) / CLOCKS_PER_SEC;
        double ops_per_sec = double(count - last_count) / tot_time;
        uint32_t eta = (total - count) / ops_per_sec;
        cerr << msg << ": " << count << "/" << total
        << "  ops/sec:" << ops_per_sec
        << " ETA:" << eta << eol << flush;
        last_time = end;
        last_count = count;
    }

    void finish() {
        clock_t end = clock();
        double tot_time = double(end - start_time) / CLOCKS_PER_SEC;
        double ops_per_sec = double(count) / tot_time;
        cerr << msg << ": finished. ops/sec:" << ops_per_sec
        << " time(s):" << tot_time << eol << flush;
        if (isatty(fileno(stdout))) {
            // in a terminal we need a read \n
            cerr << endl;
        }
    }
};

class MinMax {
    double _min = numeric_limits<double>::infinity();
    double _max = -numeric_limits<double>::infinity();
    int32_t _offset;
public:
    void feed(double v) {
        if (v > _max) {
            _max = v;
        }
        if (v < _min) {
            _min = v;
        }
    }
    double min() const { return _min; }
    double max() const { return _max; }
    int encode(double value) const {
        int encoded  = round(value / m_by_g - _offset);
        return encoded;
    }
    double decode(int value) const {
        double decoded =  m * (granularity * (_offset + static_cast<double>(value)));
        return decoded;
    }
    void calc_offset() {
        _offset = round((min() + max()) / (2.0 * m_by_g));
        cerr << "Calculating offset. min:" << _min
             << " max:" << _max << " offset:" << _offset << endl;
    }
    int32_t offset() const {
        return _offset;
    }
};

class RoutingNode : public degree_point {
    const double EARTH_RADIUS = 6372798.2;
    int mapped_id_;
public:
    RoutingNode (int mapped_id, double lon, double lat)
    : degree_point(lon, lat)
    , mapped_id_(mapped_id)
    {}
    double lon() const {
        return this->get<0>();
    }
    double lat() const {
        return this->get<1>();
    }
    int32_t mapped_id() const {
        return mapped_id_;
    }
    double distance(const RoutingNode& other) const {
        return bg::distance(static_cast<const degree_point&>(*this),
                            static_cast<const degree_point&>(other)) * EARTH_RADIUS;
    }
};

struct RoutingWay {
    bool one_way;
    vector<int64_t> nodes;
};

class RoutingGraphGenerator {
public:
    RoutingGraphGenerator(std::ofstream& _os)
    : os(_os) {
        // place holder for inserting the header
        for (int i=0; i < 30; ++i) {
            os << " ";
        }
        os << "\n";

    }
    void insert(const RoutingNode& n1, const RoutingNode& n2, double cost, bool oneWay) {
        ++edges_;
        string oneWayStr = oneWay ? "1" : "0";
        os << n1.mapped_id() << " "
           << n2.mapped_id() << " "
           << encode_cost(cost) << " "
           << oneWayStr << "\n";
    }
    void finish(uint32_t nodes) {
        // creates the header for the file
        os.seekp(0);
        os << "d\n";
        os << nodes << " " << edges_ << "\n";
        os.close();
    }
    uint32_t edges() {
        return edges_;
    }
private:
    std::ofstream& os;
    uint32_t edges_ = 0;
};

// Maps: pbf-node-id -> sequential-id
class NodeMapper {
    typedef std::unordered_map<int64_t, RoutingNode> map_type_t;
    typedef map_type_t::const_iterator node_iterator_t;
public:

    const RoutingNode& get(int64_t nodeId) const {
        return nodeIdMapper.find(nodeId)->second;
    }

    void insert(int64_t nodeId, double lon, double lat) {
        int32_t mappedId = get_mapped_id(nodeId);
        nodeIdMapper.insert({nodeId, RoutingNode(mappedId, lon, lat)});
    }

    bool exists(int64_t nodeId) const {
        return nodeIdMapper.find(nodeId) != nodeIdMapper.end();
    }

    uint32_t nodes_count() const {
        return nodeIdMapper.size();
    }

    node_iterator_t begin() const {
        return nodeIdMapper.begin();
    }

    node_iterator_t end() const {
        return nodeIdMapper.end();
    }

private:
    int32_t get_mapped_id(int64_t nodeId) {
        node_iterator_t it = nodeIdMapper.find(nodeId);
        if (it != nodeIdMapper.end()) {
            return it->second.mapped_id();
        }
        return nodeIdMapper.size() + 1;
    }
    map_type_t nodeIdMapper;
};

static bool is_highway(const osmpbf::IWayStream& way) {
    for (int i=0; i < way.tagsSize(); ++i) {
        if ("highway" == way.key(i)){
            return true;
        }
    }
    return false;
}

template<typename Func>
static void iterate_file(osmpbf::OSMFileIn* in, Func fun) {
    in->reset();
}

void read_file(osmpbf::OSMFileIn* in,
               unordered_map<int64_t,RoutingNode>* all_nodes,
               unordered_map<int64_t,RoutingWay>* all_ways)
{
    using namespace osmpbf;

    PrimitiveBlockInputAdaptor pbi;
    StatusUpdater upd;
    upd.start("read nodes and ways from file (kb)", in->totalSize() / 1024);
    while (in->parseNextBlock(pbi)) {
        upd.absolute(in->dataPosition() / 1024);
        for (osmpbf::INodeStream node = pbi.getNodeStream(); !node.isNull(); node.next()) {
            all_nodes->insert({node.id(), RoutingNode(-1, node.lond(), node.latd())});
        }
        for (osmpbf::IWayStream way = pbi.getWayStream(); !way.isNull(); way.next()) {
            if (!is_highway(way)) {
                // cerr << "Skipping way because it's not a highway " << way.id() << endl;
                continue;
            }
            vector<int64_t> v;
            for(RefIterator refIt = way.refBegin(); refIt != way.refEnd(); ++refIt) {
                v.push_back(*refIt);
            }
            bool one_way = "yes" == way.valueByKey("oneway");
            all_ways->insert({way.id(), {one_way, v}});
        }
    }
    upd.finish();

    cerr << "all_nodes:" << all_nodes->size()
         << " all_ways:" << all_ways->size() << endl;
}

void enumerate_valid_ways(const unordered_map<int64_t,RoutingNode>& all_nodes,
                          const unordered_map<int64_t,RoutingWay>& all_ways,
                          unordered_set<int64_t>* valid_ways,
                          unordered_set<int64_t>* nodes_used)
{
    using namespace osmpbf;
    cerr << "SECOND STEP started: enumerate valid ways" << endl;
    for (const auto& waykv : all_ways) {
        int64_t way_id = waykv.first;
        bool valid_way = true;
        const vector<int64_t>& v = waykv.second.nodes;
        for (const auto& node_id : v) {
            if (all_nodes.count(node_id) == 0) {
                cerr << "WARNING: Invalid way: " << way_id << endl;
                valid_way = false;
                break;
            }
        }

        if (!valid_way) {
            break;
        }

        valid_ways->insert(way_id);
        for (const auto& node_id : v) {
            nodes_used->insert(node_id);
        }
    }
    cerr << "SECOND STEP finished: valid_ways: " << valid_ways->size()
         << " nodes_used:" << nodes_used->size()
         << endl;
}

void feed_node_mapper(const unordered_map<int64_t,RoutingNode>& all_nodes,
                      const unordered_set<int64_t>& nodes_used,
                      NodeMapper* nodeMapper,
                      MinMax* lonMinMax,
                      MinMax* latMinMax)
{
    using namespace osmpbf;
    cerr << "THIRD STEP started: retrieve detailed information about important nodes" << endl;
    for (const auto& node_id : nodes_used) {
        const RoutingNode& node = all_nodes.at(node_id);
        lonMinMax->feed(node.lon());
        latMinMax->feed(node.lat());
        nodeMapper->insert(node_id, node.lon(), node.lat());
    }
}

static void
generate_ddsg(const string& ddsg_fn,
              const unordered_map<int64_t,RoutingWay>& all_ways,
              const unordered_set<int64_t>& valid_ways,
              const NodeMapper& node_mapper)
{
    using namespace osmpbf;
    cerr << "FOURTH STEP: create the ddsg graph" << endl;

    ofstream ddsg(ddsg_fn, std::ios::out | std::ios::trunc);
    if (!ddsg) {
        cerr << "Failed to open: " << ddsg_fn << endl;
        exit(1);
    }

    StatusUpdater upd;
    upd.start("Generating the ddsg file", all_ways.size());
    RoutingGraphGenerator graphGenerator(ddsg);
    for (const auto& waykv : all_ways) {
        upd.progress();
        int64_t way_id = waykv.first;
        if (valid_ways.count(way_id) == 0) {
            continue;
        }
        const RoutingWay& routing_way = waykv.second;
        const vector<int64_t>& way_nodes = routing_way.nodes;
        for (int i=1; i < way_nodes.size(); ++i) {
            const RoutingNode& p1 = node_mapper.get(way_nodes[i-1]);
            const RoutingNode& p2 = node_mapper.get(way_nodes[i]);
            double cost = p1.distance(p2);

            graphGenerator.insert(p1, p2, cost, routing_way.one_way);
        }
    }
    graphGenerator.finish(node_mapper.nodes_count());
    upd.finish();

    cerr << "FOURTH STEP finished: Total edges:" << graphGenerator.edges() << endl;
}

static void
exec_sql(sqlite3* db, const char* sql) {
    char *errmsg;
    sqlite3_exec(db, sql, NULL, NULL, &errmsg);
    if (errmsg) {
        cerr << "Error executing '" << sql << ": " << errmsg << endl;
        sqlite3_free(errmsg);
    }
}

static void
insert_coords(sqlite3* db,
              const NodeMapper& node_mapper,
              const MinMax& lon_min_max,
              const MinMax& lat_min_max)
{
    sqlite3_stmt *stmt;
    int errc = sqlite3_prepare_v2(db, "INSERT INTO coords VALUES(?1,?2,?2,?3,?3)", -1, &stmt, NULL);
    if (stmt == NULL) {
        cerr << "Error preparing statement: insert_coords: " << errc << endl;
        return;
    }

    StatusUpdater upd;
    upd.start("inserting coords into the r-tree table", node_mapper.nodes_count());
    for (const auto& it : node_mapper) {
        const RoutingNode& routingNode = it.second;

        sqlite3_bind_int(stmt, 1, routingNode.mapped_id());
        sqlite3_bind_int(stmt, 2, lon_min_max.encode(routingNode.lon()));
        sqlite3_bind_int(stmt, 3, lat_min_max.encode(routingNode.lat()));
        sqlite3_step(stmt);
        sqlite3_reset(stmt);

        upd.progress();
    }
    sqlite3_finalize(stmt);
    upd.finish();
}

static void
insert_nodes_map(sqlite3* db, const NodeMapper& node_mapper) {
    sqlite3_stmt *stmt;
    int errc = sqlite3_prepare_v2(db, "INSERT INTO nodes_map VALUES(?1,?2)", -1, &stmt, NULL);
    if (stmt == NULL) {
        cerr << "Error preparing statement insert_nodes_map: " << errc << endl;
        return;
    }
    StatusUpdater upd;
    upd.start("inserting data into the nodes-map table", node_mapper.nodes_count());
    for (const auto& it : node_mapper) {
        const RoutingNode& routingNode = it.second;

        sqlite3_bind_int64(stmt, 1, it.first);
        sqlite3_bind_int(stmt, 2, routingNode.mapped_id());
        sqlite3_step(stmt);
        sqlite3_reset(stmt);

        upd.progress();
    }
    sqlite3_finalize(stmt);
    upd.finish();

}

static void
insert_param(sqlite3_stmt *stmt, const char* key, double val) {
    // in sqlite we can store anything we want. it's very flexible
    sqlite3_bind_text(stmt, 1, key, -1, SQLITE_STATIC);
    sqlite3_bind_double(stmt, 2, val);
    sqlite3_step(stmt);
    sqlite3_reset(stmt);
}


static void
insert_params(sqlite3* db,
              const MinMax& lon_min_max,
              const MinMax& lat_min_max) {

    sqlite3_stmt *stmt;
    int errc = sqlite3_prepare_v2(db, "INSERT INTO params VALUES(?1,?2)", -1, &stmt, NULL);
    if (stmt == NULL) {
        cerr << "Error preparing statement insert_nodes_map: " << errc << endl;
        return;
    }

    insert_param(stmt, "m", m);
    insert_param(stmt, "granularity", granularity);
    insert_param(stmt, "lon_min", lon_min_max.min());
    insert_param(stmt, "lon_max", lon_min_max.max());
    insert_param(stmt, "lat_min", lat_min_max.min());
    insert_param(stmt, "lat_max", lat_min_max.max());

    sqlite3_finalize(stmt);
}

static void
save_nodes_db(const string& fn,
              const NodeMapper& node_mapper,
              const MinMax& lon_min_max,
              const MinMax& lat_min_max)
{
    // delete any old db
    unlink(fn.c_str());
    sqlite3 * db;

    sqlite3_open(fn.c_str(), &db);
    exec_sql(db, "PRAGMA synchronous = OFF");
    exec_sql(db, "PRAGMA journal_mode = OFF");

    exec_sql(db, "BEGIN TRANSACTION");
    exec_sql(db, "CREATE VIRTUAL TABLE coords USING rtree_i32(id,ilon,alon,ilat,alat)");
    exec_sql(db, "CREATE TABLE nodes_map(old_id,new_id)");
    exec_sql(db, "CREATE TABLE params(key, value)");

    insert_coords(db, node_mapper, lon_min_max, lat_min_max);
    insert_nodes_map(db, node_mapper);
    insert_params(db, lon_min_max, lat_min_max);

    exec_sql(db, "END TRANSACTION");
    sqlite3_close(db);

}

int main(int argc, char ** argv) {
    using namespace osmpbf;
    if (argc < 4) {
        cerr << "usage: " << argv[0] << " <in> <out-nodes> <out-ddsg>" << endl;
        return 1;
    }

    string input_fn(argv[1]);
    string nodes_fn(argv[2]);
    string ddsg_fn(argv[3]);

    cerr << "input file" << input_fn << endl;

    osmpbf::OSMFileIn in(input_fn, false);
    if (!in.open()) {
        cerr << "Failed to open " <<  input_fn << endl;
        return -1;
    }

    unordered_map<int64_t,RoutingNode> all_nodes;
    unordered_map<int64_t,RoutingWay> all_ways;
    unordered_set<int64_t> valid_ways, nodes_used;

    read_file(&in, &all_nodes, &all_ways);

    enumerate_valid_ways(all_nodes, all_ways, &valid_ways, &nodes_used);

    NodeMapper node_mapper;
    MinMax latMinMax, lonMinMax;
    feed_node_mapper(all_nodes, nodes_used, &node_mapper, &lonMinMax, &latMinMax);

    generate_ddsg(ddsg_fn, all_ways, valid_ways, node_mapper);

    // Generate the nodes structure
    lonMinMax.calc_offset();
    latMinMax.calc_offset();

    cerr << "Creating nodes file: " << nodes_fn << endl;
    save_nodes_db(nodes_fn, node_mapper, lonMinMax, latMinMax);

    cerr << "Conversion process finished successfully" << endl;
    
    return 0;
}
