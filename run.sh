INPUT=maps/monaco-latest.osm.pbf

MODE=debug
if [ ${#@} -eq 0 ]; then
    :
elif [ ${#@} -eq 1 ]; then
    INPUT=$1
elif [ ${#@} -eq 2 ]; then
    INPUT=$1
    MODE=$2
else
    echo "usage $0 [map.osm.pbf [debug|release]]"
    exit 1
fi

if [[ $MODE != "debug" && $MODE != "release" ]]; then
    echo "error: mode $MODE must be debug or release. invalid: $MODE"
    exit 1
fi

echo "Generating graph: input $INPUT in $MODE"


CONVERTERS_PATH=./$MODE/converters
SOLVER_PATH=./$MODE/src

./$MODE/pbf2routing/pbf2routing $INPUT $INPUT.sqlite.db $INPUT.ddsg

echo "Execution finished of:$INPUT outputs:$TMP"
