#!/bin/bash
hg clone https://bitbucket.org/eigen/eigen/
mv eigen ~/
cd ~/eigen
echo $PWD
export EIGEN3_INCLUDE_DIR="$PWD"
cd -

