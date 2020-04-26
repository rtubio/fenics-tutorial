#!/bin/bash

source "conf/project.conf"
CWD="$(pwd)"
USER="$(whoami)"
echo ">>> Working directory: $CWD"

echo "### (0) Configuring operating system"
[[ -f "$CONF_DEBPKS" ]] && {
    sudo apt-get install $(grep -vE "^\s*#" $CONF_DEBPKS | tr "\n" " ")
} || {
    echo "No <$CONF_DEBPKS> available, skipping OS packages installation"
}

echo "### (1) Configuring Python environment"
virtualenv --python "$PYENV_VERSION" "$PYENV_DIR"
source "$PYENV_ACTIVATE"

[[ -f "$CONF_PYPKS" ]] && {
    pip install -r "$CONF_PYPKS"
} || {
    echo "No <$CONF_PYPKS> available, skipping PYTHON packages installation"
}

[[ ! -d "$SRCBUILD_DIR" ]] && {
    sudo mkdir -p "$SRCBUILD_DIR"
    sudo chown -R "$USER:$USER" "$SRCBUILD_DIR"
} || {
    echo "<$SRCBUILD_DIR> exists, skipping creation"
}

### Move to the source build dir to build the required libraries (dolfin, pybind, mshr)
cd "$SRCBUILD_DIR"

# echo "### (2) Pybind 11"
# PYBIND11_VERSION=2.2.3
# PYBIND11_TARGZ="$SRCBUILD_DIR/pybind11-$PYBIND11_VERSION.tar.gz"
# wget -nc --quiet\
#     "https://github.com/pybind/pybind11/archive/v${PYBIND11_VERSION}.tar.gz"\
#     -O "$PYBIND11_TARGZ"
# 
# tar -xf "$PYBIND11_TARGZ" && cd pybind11-${PYBIND11_VERSION}
# mkdir build
# cd build && cmake -DPYBIND11_TEST=off .. && sudo make install && cd ../..

echo "### (3) Configuring FENICS"
FENICS_VERSION=$(python -c"import ffc; print(ffc.__version__)") && \
    echo ">> FENICS_VERSION = $FENICS_VERSION"

echo "### (4) Cloning Dolphin"
git clone --branch=$FENICS_VERSION https://bitbucket.org/fenics-project/dolfin
echo "### (5) Building Dolphin"
mkdir dolfin/build
cd dolfin/build && cmake .. && sudo make install && cd ../..
cd dolfin/python && pip install . && cd ../..

# latest available branch is too old, getting master
# --branch=$FENICS_VERSION
echo "### (5) Cloning MSHR"
git clone https://bitbucket.org/fenics-project/mshr
git checkout "tags/$FENICS_VERSION"
echo "### (6) Building MSHR"
mkdir mshr/build
cd mshr/build && cmake .. && sudo make install && cd ../..
cd mshr/python   && pip install . && cd ../..

cd "$CWD"

deactivate
