pbf2routing
===========

This is a simple C++ library to parse OpenStreetMap's PBF files and generate
a ddsg file.

How to compile and run this version:

    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake ..   # alternatively you can use cmake -G Xcode ..
    make
    cd ..
    sh run.sh

If you don't use Xcode, you might need to change the run.sh script.

To clean the project just delete the `build` directory.
