FROM opensuse/leap:15.1

RUN zypper ref && zypper up -l -y

RUN zypper in -l -y gcc8 gcc8-c++ eigen3-devel tbb-devel libboost_program_options*-devel libboost_filesystem*-devel git cmake hdf5-devel

WORKDIR /usr/src

RUN git clone --recurse-submodules -b master https://github.com/simonraschke/vesicle2.git

WORKDIR /usr/src/vesicle2

RUN mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_CXX_COMPILER=g++-8 && make -j

ENTRYPOINT ["/root/bin/vesicle2"]