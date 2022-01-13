# HasegawaMima
2D Hasegawa-Mima solver in Python and C++

## Setup
```
sudo apt update
sudo apt install build-essential
sudo apt install libboost-all-dev
sudo apt install libfftw3-dev
sudo apt install snapd
sudo snap install cmake --classic
sudo apt update
sudo apt upgrade
```

## Configure and build
```
git clone https://github.com/stirkas/HasegawaMima.git
cd HasegawaMima/C
cmake . -DCMAKE_BUILD_TYPE=Release -B build
cmake --build build --target install
```
