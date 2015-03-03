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
#include <osmpbf/parsehelpers.h>
#include <osmpbf/inode.h>
#include <osmpbf/iway.h>
#include <osmpbf/irelation.h>
#include <osmpbf/filter.h>


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

void parseBlock(osmpbf::PrimitiveBlockInputAdaptor & pbi) {
	//using filters is straight forward
	osmpbf::MultiStringTagFilter * hwFilter = new osmpbf::MultiStringTagFilter("highway");
	hwFilter->setValues(std::set<std::string>({"path", "residential"}));
	//andFilter takes ownership of hwFilter
	osmpbf::AndTagFilter andFilter;
	andFilter.addChild(new osmpbf::KeyOnlyTagFilter("name"));
	andFilter.addChild(hwFilter);

	//build the id cache for faster queries (this is not neccessary)
	if (!andFilter.buildIdCache()) {
// 		std::cout << "No matching elements in this block" << std::endl;
	}
	
	if (pbi.nodesSize()) {
		for (osmpbf::INodeStream node = pbi.getNodeStream(); !node.isNull(); node.next()) {
			std::cout << "<node id=" << node.id() << " lat=" << node.latd() << " lon=" << node.lond()  << ">"<< std::endl;
			for(uint32_t i = 0, s = node.tagsSize();  i < s; ++i) {
				std::cout << "\t<tag k=" << node.key(i) << " v=" << node.value(i) << ">" << std::endl;
			}
			std::cout << "</node>" << std::endl;
			//to check if a primitive matches the filter use
			if (andFilter.matches(node)) {
				;//do something
			}
		}
	}
		
	if (pbi.waysSize()) {
		for (osmpbf::IWayStream way = pbi.getWayStream(); !way.isNull(); way.next()) {
			std::cout << "<way id=" << way.id() << ">" << std::endl;
			for(osmpbf::IWay::RefIterator refIt(way.refBegin()), refEnd(way.refEnd()); refIt != refEnd; ++refIt) {
				std::cout << "\t<nd ref=" << *refIt << "/>" << std::endl;
			}
			for(uint32_t i = 0, s = way.tagsSize();  i < s; ++i) {
				std::cout << "\t<tag k=" << way.key(i) << " v=" << way.value(i) << ">" << std::endl;
			}
			std::cout << "</way>" << std::endl;
		}
	}
	
	if (pbi.relationsSize()) {
		for(osmpbf::IRelationStream relation = pbi.getRelationStream(); !relation.isNull(); relation.next()) {
			std::cout << "<relation id=" << relation.id() << ">" << std::endl;
			for(osmpbf::IMemberStream mem(relation.getMemberStream()); !mem.isNull(); mem.next()) {
				std::cout << "\t<member type=" << primitiveTypeToString( mem.type() ) << " ref=" << mem.id() << " role=" << mem.role() << "/>" <<  std::endl;
			}
			for(uint32_t i = 0, s = relation.tagsSize();  i < s; ++i) {
				std::cout << "\t<tag k=" << relation.key(i) << " v=" << relation.value(i) << ">" << std::endl;
			}
			std::cout << "</relation>" << std::endl;
		}
	}
}

int main(int argc, char ** argv) {
	if (argc < 3) {
		std::cout << "Need parse type and in file" << std::endl;
		std::cout << "Parse type may be any of s=single threaded,o=OpenMP,c=C++11 Threads" << std::endl;
	}
	
	std::string parseType(argv[1]);
	std::string inputFileName(argv[2]);

	osmpbf::OSMFileIn inFile(inputFileName, false);

	if (!inFile.open()) {
		std::cout << "Failed to open " <<  inputFileName << std::endl;
		return -1;
	}
	auto parseFunc = [](osmpbf::PrimitiveBlockInputAdaptor & pbi){parseBlock(pbi);};
	
	if (parseType == "s") {
		osmpbf::parseFile(inFile, parseFunc);
	}
	else if (parseType == "o") {
		osmpbf::parseFileOmp(inFile, parseFunc);
	}
	else if (parseType == "c") {
		osmpbf::parseFileCPPThreads(inFile, parseFunc);
	}

	return 0;
}
