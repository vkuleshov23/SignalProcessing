#!/bin/bash

gnuplot -persist <<-EOFMarker
 set title "$2"
 plot "$1" using 1:2 with lines