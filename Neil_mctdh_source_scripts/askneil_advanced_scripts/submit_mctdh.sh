#!/bin/bash

# load mctdh
module load mctdh/84.16

# call mctdh
mctdh84 -mnd -w ${1}

echo "Finished execution!"

