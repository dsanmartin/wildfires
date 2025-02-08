#!/bin/bash
# Run on background no hangup
# nohup ./bin/wildfires data/input/parameters.txt > wildfire.out 2>&1 &
# nohup ./bin/wildfires_cu data/input/f19_gaussian.txt &
nohup ./bin/wildfires_cu data/input/debug.txt > variable.out 2>&1 &