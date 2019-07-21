#!/bin/bash

rm *.txt *.png

root -l -b << EOF
.L Make2Dplots.C
Make2Dplots a
a.Loop(1)
EOF

root -l -b << EOF
.L Make2Dplots_eta.C
Make2Dplots a
a.Loop(1)
EOF
