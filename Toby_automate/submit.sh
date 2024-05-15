#!/bin/bash

# load mctdh
# module load mctdh/84.16

echo "Running MCTDH-86!"

# call mctdh
# mctdh84 -mnd -w ${1}
/home/bjb2chen/LOCAL/mctdh/mctdh86.4/bin/binary/x86_64/mctdh86 -mnd -w ${1}

echo "Finished execution!"

