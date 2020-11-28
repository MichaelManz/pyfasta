#!/bin/bash -e

cd /data
FILE=$(ls ./*.fasta)
OUTPUT_DIR_SPLIT=$FILE.split
mkdir "$OUTPUT_DIR_SPLIT"
cd "$OUTPUT_DIR_SPLIT"

SIZE_LIMIT=${SIZE_LIMIT:-10000}

echo "Running pyfasta split on $FILE into $OUTPUT_DIR_SPLIT"
pyfasta split --header '%(seqid)s' -r -s $SIZE_LIMIT "../$FILE"

NUM_PROC=$(grep -c ^processor /proc/cpuinfo)

OUTPUT_DIR_FOLD=$FILE.fold
mkdir "../$OUTPUT_DIR_FOLD"
cd "../$OUTPUT_DIR_FOLD"

echo "Running RNAfold on $OUTPUT_DIR_SPLIT/* into $OUTPUT_DIR_FOLD"
find ../"$OUTPUT_DIR_SPLIT/" -type f -print0 | xargs -0 -P $NUM_PROC -I {} /bin/bash -c 'if [ -f "$(basename {}).fold" ]; then printf "x"; else RNAfold -o {} -i {};printf "o"; fi'
echo "..done"

OUTPUT_DIR_CT=$FILE.ct
mkdir "../$OUTPUT_DIR_CT"
cd "../$OUTPUT_DIR_CT"
echo "Running b2ct on $OUTPUT_DIR_FOLD/* into $OUTPUT_DIR_CT"
find ../"$OUTPUT_DIR_FOLD/" -type f -name "*fold" -print0 | xargs -0 -P $NUM_PROC -I {} /bin/bash -c 'if [ -f "$(basename {}).ct" ]; then printf "x"; else b2ct < "{}" > "$(basename {}).ct";printf "o"; fi'
echo "..done"

OUTPUT_DIR_FRAMES=$FILE.frames
mkdir "../$OUTPUT_DIR_FRAMES"
cd "../$OUTPUT_DIR_FRAMES"
echo "Running ct2bFrames on $OUTPUT_DIR_CT/* into $OUTPUT_DIR_FRAMES"
find ../"$OUTPUT_DIR_CT/" -type f -print0 | xargs -0 -P $NUM_PROC -I {} /bin/bash -c 'if [ -f "$(basename {}).ss.b" ]; then printf "x"; else perl /pyfasta/scripts/ct2bFrames.pl < "{}" > "$(basename {}).ss.b";printf "o"; fi'
echo "..done"

OUTPUT_DIR_STATS=$FILE.stats
mkdir "../$OUTPUT_DIR_STATS"
cd "../$OUTPUT_DIR_STATS"
echo "Running ct2bFrameStats on $OUTPUT_DIR_CT/* into $OUTPUT_DIR_STATS"
find ../"$OUTPUT_DIR_CT/" -type f -print0 | xargs -0 -P $NUM_PROC -I {} /bin/bash -c 'if [ -f "$(basename {}).txt" ]; then printf "x"; else perl /pyfasta/scripts/ct2bFrameStats.pl < "{}" > "$(basename {}).txt";printf "o"; fi'
echo "..done"