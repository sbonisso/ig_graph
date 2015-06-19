#!/bin/bash

./iggraph -v data/igh_refs_simple/human_IGHV.fa \
    -d data/igh_refs_simple/human_IGHD.fa \
    -j data/igh_refs_simple/human_IGHJ.fa \
    -V 21 \
    -D 10 \
    -J 21 \
    -s 0 \
    -f \
    $@ 
