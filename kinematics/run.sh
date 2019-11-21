#!/bin/bash

root -l -b << EOF
.L test.C
test a
a.Loop()
EOF
