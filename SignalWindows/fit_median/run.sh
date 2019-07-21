#!/bin/bash

rm *.txt *.png

root -l -b << EOF
.L Make2Dplots.C
Make2Dplots a
a.Loop(1)
EOF
