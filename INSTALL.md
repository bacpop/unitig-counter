# Installation instructions
unitig-counter depends on boost. You may have this on your system
already, but it is probably to install a specific minimal version (which
will be statically linked, and therefore not conflict with existing
system installs) using the following instructions.

## boost
First, install boost:
```
git clone --depth 1 --single-branch --branch boost-1.69.0 https://github.com/boostorg/boost.git boost
cd boost
rmdir libs/regex libs/filesystem libs/system
git clone --depth 50 https://github.com/boostorg/filesystem.git libs/filesystem
git clone --depth 50 https://github.com/boostorg/system.git libs/system
git clone --depth 50 https://github.com/boostorg/regex.git libs/regex
git submodule update -q --init libs/algorithm
git submodule update -q --init tools/boostdep
git submodule update -q --init tools/build
git submodule update -q --init tools/inspect
python2 tools/boostdep/depinst/depinst.py regex --include example
python2 tools/boostdep/depinst/depinst.py filesystem --include example
python2 tools/boostdep/depinst/depinst.py system --include example
mkdir ../boost_built
./bootstrap.sh --prefix=../boost_built --with-libraries=regex,filesystem,system
./b2 variant=release architecture=x86 debug-symbols=off threading=multi runtime-link=shared link=static,shared cxxflags=-fPIC install
```

## unitig-counter
Now run cmake to install unitig-counter. First set your install location
e.g. `export PREFIX=$HOME`. Then run:
```
mkdir build && pushd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DBoost_INCLUDE_DIR=$PWD/boost_built/include -DBoost_LIBRARY_DIR=$PWD/boost_built/lib ..
make VERBOSE=1
make install
popd build
```

If that works, and `$PREFIX/bin` is on your path, you should be able to
run:
```
unitig-counter -version
unitig-counter -help
```
