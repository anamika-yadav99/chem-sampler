CHECKPOINT_DIR=$1
SAMPLES=$2

if [ $# -eq  ] ; then
    echo "Filename not supplied."
else

python molecule_generation.sample.py $CHECKPOINT_DIR $SAMPLES

fi
