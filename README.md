# Introduction
The CZ project, the TSMFVC project and the TSNIVC project include benchmarking of the CZ scheme, the TSMFVC scheme and the TSNIVC scheme.
The CZ project, the TSMFVC project, the TSNIVC project, the ZW1 project, and the ZW2 project include the performance of these schemes in evaluating multivariate polynomial functions.  
The DecisionTree_BGV_CKKS project, the DecisionTree_CZ project and the DecisionTree_TSMFVC project show the efficiency of BGV, CKKS, CZ and TSMFVC for evaluating decision trees on private data.

# Dependencies  
## Install m4
    wget http://mirrors.kernel.org/gnu/m4/m4-1.4.13.tar.gz
    tar -xzvf m4-1.4.13.tar.gz
    cd m4-1.4.13
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install autoconf
    wget http://mirrors.kernel.org/gnu/autoconf/autoconf-2.65.tar.gz
    tar -xzvf autoconf-2.65.tar.gz
    cd autoconf-2.65
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install automake
    wget http://mirrors.kernel.org/gnu/automake/automake-1.11.tar.gz
    tar xzvf automake-1.11.tar.gz
    cd automake-1.11
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install libtool
    wget http://mirrors.kernel.org/gnu/libtool/libtool-2.2.6b.tar.gz
    tar xzvf libtool-2.2.6b.tar.gz
    cd libtool-2.2.6b
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install GMP
    tar -jxvf  gmp-6.2.1.tar.bz2
    cd gmp-6.2.1
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
    make check
## Install CMake
    tar -xzvf cmake-3.10.2.tar.gz
    cd cmake-3.10.2
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install NTL
    tar -xzvf ntl-11.5.1.tar.gz
    cd ntl-11.5.1/src
    ./configure 
    make
    sudo make install
##  Install flint
    sudo apt install -y libflint-dev
## Install pbc
    wget https://crypto.stanford.edu/pbc/files/pbc-0.5.14.tar.gz
    tar -xzf pbc-0.5.14.tar.gz
    cd pbc-0.5.14
    ./configure --prefix=/usr/local
    make -j$(nproc)
    sudo make install
## Install relic
    git clone https://github.com/relic-toolkit/relic.git
    cd relic
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
    make -j$(nproc)
    sudo make install
## Install mpfr
    wget https://www.mpfr.org/mpfr-current/mpfr-4.2.0.tar.gz
    tar -xzf mpfr-4.2.0.tar.gz
    cd mpfr-4.2.0
    ./configure --prefix=/usr/local --with-gmp=/usr
    make -j$(nproc)
    sudo make install
## Install libsodium
    sudo apt install -y libsodium-dev
## Install HElib
https://github.com/homenc/HElib
"# paper_code" 
