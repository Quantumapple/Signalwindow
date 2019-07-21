#!/bin/bash

rm *.txt *.png

root -l -b << EOF
.L test.C
test a
a.Loop(6)
EOF
