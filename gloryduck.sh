source ../../build/bin/thisroot.sh

export GLIKESYS="/home/tjark/darkmatter/gloryduck/gLike/"

export DYLD_LIBRARY_PATH="${GLIKESYS}:${DYLD_LIBRARY_PATH}"

make

root


