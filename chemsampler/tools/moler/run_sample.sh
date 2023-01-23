CHECKPOINT_DIR=$1
SAMPLES=$2
DATA_FILE=$3

if [ $# -eq 2 ] ; then
    echo "Filename not supplied."
else

python molecule_generation/cli/sample.py $CHECKPOINT_DIR $SAMPLES $DATA_FILE

fi
