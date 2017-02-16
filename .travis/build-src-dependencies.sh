#!/bin/bash

REMOTE_SOURCE="http://software.ligo.org/lscsoft/source"

# set prefix and pkgconfig path
PREFIX=`python -c "import sys; print(sys.prefix)"`
export PKG_CONFIG_PATH=${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}

# build a package using autotools configure/make/make install
build-package() {
    # get tarball from remote
    local remote=$1
    shift
    local tarball=`basename $remote`
    local tag="${tarball%.*}"
    [[ "${tag}" == *".tar" ]] && tag="${tag%.*}"
    wget $remote -O $tarball --quiet
    # unpack
    mkdir -p $tag
    tar -xf $tarball --strip-components=1 -C $tag
    cd $tag
    # build and install
    ./configure --enable-silent-rules --quiet --prefix=$PREFIX $@
    make -j2
    make install
    # clean up
    cd -
    rm -rf $tag $tarball
}

# build NDS2
build-package ${REMOTE_SOURCE}/nds2-client-${NDS2_CLIENT_VERSION}.tar.gz --disable-swig-java --disable-mex-matlab

# build LAL
build-package ${REMOTE_SOURCE}/lalsuite/lal-${LAL_VERSION}.tar.xz --enable-swig-python

# build LALFrame
build-package ${REMOTE_SOURCE}/lalsuite/lalframe-${LALFRAME_VERSION}.tar.xz --enable-swig-python
