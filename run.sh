#!/bin/bash

rm -r build/
rm output/*
clear

time fpm --flag "-O3 -Wextra -finit-real=zero -finit-integer=0 -fcheck=all -ffpe-trap=invalid " run

