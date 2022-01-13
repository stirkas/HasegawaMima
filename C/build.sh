#!/bin/sh
rm -r bin
rm -r build
cmake . -DCMAKE_BUILD_TYPE=Release -B build
cmake --build build --target install
