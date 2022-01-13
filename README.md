# HasegawaMima
2D Hasegawa-Mima solver in Python and C++

## Setup
```
sudo apt update
sudo apt install build-essential
sudo apt install libboost-all-dev
sudo apt install fftw3
sudo apt install snapd
sudo snap install cmake --classic
sudo apt update
sudo apt upgrade
```

## Configure and build
```
git clone https://github.com/stirkas/HasegawaMima.git
cd HasegawaMima/C
cmake . -B build
cmake --build build --target install
```
