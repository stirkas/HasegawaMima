if [$1=="clean"]
then
   rm -r bin
   rm -r build
fi

cmake . -DCMAKE_BUILD_TYPE=Release -B build
cmake --build build --target install
