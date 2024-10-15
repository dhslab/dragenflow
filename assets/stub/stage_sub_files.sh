#!/usr/bin/env bash

PWD=$(realpath .)

for i in demux_fastq/*.orig *.orig;
do
    sed -e 's!REPLACEPATH!'$PWD'!g' $i > ${i%.orig}
done