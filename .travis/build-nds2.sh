#!/bin/bash

NDS2_CLIENT_VERSION="0.10.4"

# get
wget http://software.ligo.org/lscsoft/source/nds2-client-${NDS2_CLIENT_VERSION}.tar.gz -O nds2-client-${NDS2_CLIENT_VERSION}.tar.gz --quiet
# unpack
mkdir -p nds2-client-${NDS2_CLIENT_VERSION}
tar -zxf nds2-client-${NDS2_CLIENT_VERSION}.tar.gz --strip-components=1 -C nds2-client-${NDS2_CLIENT_VERSION}
# build
cd nds2-client-${NDS2_CLIENT_VERSION}
./configure --enable-silent-rules --quiet --prefix=`python -c "import sys; print(sys.prefix)"` --disable-swig-java --disable-mex-matlab
make
# install
make install
# clean up
cd -
rm -rf nds2-client-${NDS2_CLIENT_VERSION}/
rm -f nds2-client-${NDS2_CLIENT_VERSION}.tar.gz
