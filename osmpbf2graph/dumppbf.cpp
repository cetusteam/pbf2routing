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
#if defined(__OPENMP)
#include <omp.h>
#endif
#include <osmpbf/osmfile.h>
#include <osmpbf/inode.h>
#include <osmpbf/iway.h>
#include <osmpbf/irelation.h>
#include <osmpbf/filter.h>
#include <osmpbf/primitiveblockinputadaptor.h>
#include <map>
#include <fstream>
#include <limits>

#include "../osmpbf/graph-node.pb.h"
#include "../osmpbf/routing-graph.pb.h"


#include <boost/geometry/core/cs.hpp>
#include <boost/geometry.hpp>

namespace bg = boost::geometry;

static std::map<int64_t,std::pair<double,double>> mapLatLng;
static std::pair<double,double>  maxMinLat;
static std::pair<double,double>  maxMinLng;

// To encode the lat,lon
double m = 0.000000001;
int granularity = 5000;

inline std::string primitiveTypeToString(osmpbf::PrimitiveType t) {
    switch (t) {
        case osmpbf::PrimitiveType::NodePrimitive:
            return "node";
        case osmpbf::PrimitiveType::WayPrimitive:
            return "way";
        case osmpbf::PrimitiveType::RelationPrimitive:
            return "relation";
        default:
            return "invalid";
    };
}

double distance(double lat1, double lon1, double lat2,double lon2){
    typedef bg::cs::spherical_equatorial<bg::degree> coord_system_t;
    typedef bg::model::point<double, 2, coord_system_t> degree_point;
    // use lon, lat
    degree_point p1(lon1,lat1);
    degree_point p2(lon2,lat2);
    boost::geometry::strategy::distance::haversine<double> const haversine(6372.7982);
    return bg::distance(p1, p2, haversine);
}

void updateMaxMin(std::pair<double,double> *maxMinPair,double value){
    if(value > maxMinPair->first){
        maxMinPair->first = value;
    } else if(value < maxMinPair->second ){
        maxMinPair->second = value;
    }
}

int encode(std::pair<double,double> &maxMinPair,double value){
    double offset = (maxMinPair.first - maxMinPair.second) / (2.0 * m);
    int encoded  = round( ((value / m) - offset  ) / granularity);
    return encoded;
}


double decode(std::pair<double,double> &maxMinPair,int value){
    double offset = (maxMinPair.first - maxMinPair.second) / (2.0 * m);
    double decoded =  m * (offset + (granularity * static_cast<double>(value)));
    return decoded;
}

void parseBlock(osmpbf::PrimitiveBlockInputAdaptor & pbi,RoutingGraph *routingGraph) {
    //using filters is straight forward
    osmpbf::MultiStringTagFilter * hwFilter = new osmpbf::MultiStringTagFilter("highway");
    hwFilter->setValues(std::set<std::string>({"path", "residential","service","unclassified","footway"}));
    //andFilter takes ownership of hwFilter
    osmpbf::AndTagFilter andFilter;
    andFilter.addChild(new osmpbf::KeyOnlyTagFilter("name"));
    andFilter.addChild(hwFilter);
    //build the id cache for faster queries (this is not neccessary)
    if (!andFilter.buildIdCache()) {
        // 		std::cout << "No matching elements in this block" << std::endl;
    }
    auto maxMinLat = std::make_pair(std::numeric_limits<double>::min(),std::numeric_limits<double>::max());
    auto maxMinLng = std::make_pair(std::numeric_limits<double>::min(),std::numeric_limits<double>::max());
    
    for (osmpbf::INodeStream node = pbi.getNodeStream(); !node.isNull(); node.next()) {
        
        //to check if a primitive matches the filter use
        if (andFilter.matches(node)) {
            continue;//do something
        }
        
        std::map<int64_t,std::pair<double,double>>::iterator itLatLon = mapLatLng.find(node.id());
        if(itLatLon != mapLatLng.end()){
            // Node was already processed
            continue;
        }
        //
        updateMaxMin(&maxMinLat,node.latd());
        updateMaxMin(&maxMinLng,node.lond());
        //
        auto latLon = std::make_pair(node.latd(),node.lond());
        mapLatLng.insert(std::pair<int64_t,std::pair<double,double>>(node.id(),latLon));
    }
    
    for (osmpbf::IWayStream way = pbi.getWayStream(); !way.isNull(); way.next()) {
        int64_t currentRef=-1;
        int64_t prevRef=-1;
        bool highway = false;
        std::string highwayValue = "";
        std::string oneWay = way.valueByKey("oneway");
        for(int i=0;i<way.tagsSize();i++){
            if("highway" == way.key(i)){
                highway = true;
                highwayValue = way.value(i);
                break;
            }
        }
        if(!highway ){
            continue;
        }
        if (andFilter.matches(way)) {
            //std::cout<<"FILTERED WAY "<<highwayValue<<std::endl;
            continue;//do something
        }
        
        for(osmpbf::RefIterator refIt(way.refBegin()), refEnd(way.refEnd()); refIt != refEnd; ++refIt) {
            prevRef = currentRef;
            currentRef = *refIt;
            std::map<int64_t,std::pair<double,double>>::iterator itLatLonCurrent = mapLatLng.find(currentRef);
            std::map<int64_t,std::pair<double,double>>::iterator itLatLonPrev = mapLatLng.find(prevRef);
            if(itLatLonCurrent != mapLatLng.end() && itLatLonPrev != mapLatLng.end()){
                double cost;
                cost = distance(
                                itLatLonPrev->second.first,itLatLonPrev->second.second,
                                itLatLonCurrent->second.first, itLatLonCurrent->second.second);
                //
                routingGraph->add_sources(prevRef);
                routingGraph->add_targets(currentRef);
                if(oneWay == "yes"){
                    //routingGraph->add_direction(RoutingGraph_DirectionType_FORWARD);
                }else{
                    //routingGraph->add_direction(RoutingGraph_DirectionType_BOTH);
                }
                routingGraph->add_costs(cost);
                routingGraph->add_added_by_ch(false);
                
            }
        }
        
        if (andFilter.matches(way)) {
            //std::cout<<"Matches "<<highway<<" "<<oneWay<<std::endl;
        }
    }
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        std::cout << "Need in file" << std::endl;
    }
    
    std::cout << "in file"<<argv[1] << std::endl;
    std::string inputFileName(argv[1]);
    
    osmpbf::OSMFileIn inFile(inputFileName, false);
    
    if (!inFile.open()) {
        std::cout << "Failed to open " <<  inputFileName << std::endl;
        return -1;
    }
    
    /*
     With some additionaly effort it's possible to parse the pbf in parallel which of course uses more memory
     */
    
    RoutingGraph routingGraph;
    Nodes nodesOut;
    
#if defined(__OPENMP)
    uint32_t readBlobCount = omp_get_num_procs();
    bool processedFile = false;
    while (!processedFile) {
        std::vector<osmpbf::BlobDataBuffer> pbiBuffers = inFile.getNextBlocks(readBlobCount);
        uint32_t pbiCount = pbiBuffers.size();
        processedFile = (pbiCount < readBlobCount);
#pragma omp parallel for
        for(uint32_t i = 0; i < pbiCount; ++i) {
            osmpbf::PrimitiveBlockInputAdaptor pbi(pbiBuffers[i].data, pbiBuffers[i].availableBytes);
            pbiBuffers[i].clear();
            if (pbi.isNull()) {
                continue;
            }
            parseBlock(pbi,&routingGraph);
        }
    }
#else
    osmpbf::PrimitiveBlockInputAdaptor pbi;
    while (inFile.parseNextBlock(pbi)) {
        if (pbi.isNull()){
            continue;
        }
        parseBlock(pbi,&routingGraph);
    }
#endif
    
    // Update the lat,lon
    int nodesCount=0;
    std::map<int64_t,std::pair<double,double>>::iterator itLatLon = mapLatLng.begin();
    while(itLatLon != mapLatLng.end()){
        Node *nodeGraph = nodesOut.add_nodes();
        nodeGraph->set_id(itLatLon->first);
        nodeGraph->set_new_id(nodesCount);
        if(itLatLon != mapLatLng.end()){
            nodeGraph->set_lng(encode(maxMinLng,itLatLon->second.second));
            nodeGraph->set_lat(encode(maxMinLat,itLatLon->second.first));
            
            //
            //std::cout<<""<<nodeGraph->lng()<<" "<<decode(maxMinLng, nodeGraph->lng())<<" vs "<<itLatLon->second.second<<std::endl;
            //std::cout<<""<<nodeGraph->lat()<<" "<<decode(maxMinLat, nodeGraph->lat())<<" vs "<<itLatLon->second.first<<std::endl;
        }
        nodesCount++;
        itLatLon++;
    }
    
    /**
     Writes the graph data
     */
    
    
    //Write nodes
    std::fstream outputNodes("nodes", std::ios::out | std::ios::trunc | std::ios::binary);
    if (!nodesOut.SerializeToOstream(&outputNodes)) {
        std::cerr << "Failed to write graph." << std::endl;
        return -1;
    }
    
    //Write graph
    std::fstream output("graph", std::ios::out | std::ios::trunc | std::ios::binary);
    if (!routingGraph.SerializeToOstream(&output)) {
        std::cerr << "Failed to write graph." << std::endl;
        return -1;
    }
    
    std::cout<<"O Nodes Size:: "<<nodesOut.nodes_size()<<std::endl;
    
    
    std::cout<<"Sources Size:: "<<routingGraph.sources_size()<<std::endl;
    std::cout<<"Targets Size:: "<<routingGraph.targets_size()<<std::endl;
    return 0;
}
