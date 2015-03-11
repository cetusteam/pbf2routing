# MAP_FILE=maps/delaware-latest.osm.pbf
MAP_FILE=maps/monaco-latest.osm.pbf

./build/osmpbf2graph/Debug/pbf2route \
    $MAP_FILE \
    out-nodes.pb.bin \
    out-graph.txt
