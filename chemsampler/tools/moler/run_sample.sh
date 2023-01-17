CHECKPOINT_DIR=$1
SAMPLES=$2

if [ $# -eq 1 ] ; then
    echo "Filename not supplied."
else

python molecule_generation/cli/sample.py $CHECKPOINT_DIR $SAMPLES

fi
