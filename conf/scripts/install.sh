#!/bin/bash

source "conf/project.conf"
CWD="$(pwd)"
USER="$(whoami)"
echo ">>> Working directory: $CWD"

echo "### (1) Configuring operating system"
[[ -f "$CONF_DEBPKS" ]] && {
    sudo apt-get install $(grep -vE "^\s*#" $CONF_DEBPKS | tr "\n" " ")
} || {
    echo "No <$CONF_DEBPKS> available, skipping OS packages installation"
}

echo "### (2) Configuring Python environment"
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

echo "### (3) Configuring FENICS"
FENICS_VERSION=$(python -c"import ffc; print(ffc.__version__)") && \
    echo ">> FENICS_VERSION = $FENICS_VERSION"
cd "$SRCBUILD_DIR"

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
echo "### (6) Building MSHR"
mkdir mshr/build
cd mshr/build && cmake .. && sudo make install && cd ../..
cd mshr/python   && pip install . && cd ../..

cd "$CWD"

deactivate
