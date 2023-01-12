INPUT_FILE=$1
OUTPUT_FILE=$2
TRACE_FILE=$3
CHECKPOINT_FILE=$4

if [ $# -eq 3 ] ; then
    echo "Filename not supplied."
else

bash run_processing.sh $INPUT_FILE $OUTPUT_FILE $TRACE_FILE $CHECKPOINT_FILE &&
bash run_training.sh $INPUT_FILE $OUTPUT_FILE $TRACE_FILE $CHECKPOINT_FILE


fi