#!/bin/bash

rosetta_scripts.default.macosclangrelease \
  -database $ROSETTA/database/ \
  -s 6XKI.pdb \
  -parser:protocol galigand_minimal.xml \
  -nstruct 10 \
  -gen_potential \
  -ex1 -ex2 \
  -use_input_sc \
  -extra_res_fa V6D.params \
  -out:suffix _galigand \
  -out:file:scorefile score_503A.sc \
  -ignore_unrecognized_res



