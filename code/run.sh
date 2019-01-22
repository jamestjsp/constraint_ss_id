#!/bin/bash

matlab -nodisplay -nosoftwareopengl -r \
  "disp('Running test_one.m'); \
  test_one; \
  disp('Running test_multiple.m'); \
  test_multiple; \
  disp('Running test_ineq.m'); \
  test_ineq"