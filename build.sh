#!/bin/bash


make -f Makefile_gfortran fms_build/libfms.a
make -f Makefile_gfortran MOM6_ALE/build_ale/libale.a
python setup.py config_fc --f90flags="-i4 -r8 -DPY_SOLO" --fcompiler=gfortran \
          --f90flags="-fcray-pointer -fdefault-real-8 \
          -ffixed-line-length-132 -ffree-line-length-0 -DPY_SOLO" build
python setup.py install
