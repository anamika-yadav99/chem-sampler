INPUT_FILE=$1
OUTPUT_FILE=$2
TRACE_FILE=$3
CHECKPOINT_FILE=$4

if [ $# -eq 3 ] ; then
    echo "Filename not supplied."
else

python molecule_generation/cli/preprocess.py $INPUT_FILE $OUTPUT_FILE $TRACE_FILE --pretrained-model-path $CHECKPOINT_FILE

fi
