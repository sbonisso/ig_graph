#!/bin/bash

./iggraph -v data/igh_refs_simple/mouse_IGHV.fa \
    -d data/igh_refs_simple/mouse_IGHD.fa \
    -j data/igh_refs_simple/mouse_IGHJ.fa \
    -V 21 \
    -D 10 \
    -J 21 \
    -s 0 \
    -f \
    $@
