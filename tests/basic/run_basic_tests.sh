#!/bin/bash

set -e

VAMP_BIN="../../vampire-cuda"

function run_test() {
  TESTNAME=$1
  INPUT=$2
  MAT=$3
  echo "Running $TESTNAME..."
  $VAMP_BIN < $INPUT > ${TESTNAME}_output.txt 2>&1
  if [ $? -eq 0 ]; then
    echo "$TESTNAME: PASS"
  else
    echo "$TESTNAME: FAIL"
    exit 1
  fi
}

run_test "basic_spin" basic_spin_input basic_spin.mat
run_test "basic_field" basic_field_input basic_field.mat
run_test "basic_temp" basic_temp_input basic_temp.mat

echo "All basic tests passed." 