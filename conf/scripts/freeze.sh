#!/bin/bash
#
# This script freezes the current packages installed for the Python virtual environment,
# removing the following ones:
#
#     pkg-resources, dolfin, mshr
#
# The first one is removed for compatibility reasons, while the second and third ones are removed
# since they will have to be installed from within the manually compiled pacakges, and not from
# the Pipi repository.
#
# @author: rtpardavila@gmail.com

source "conf/project.conf"

source "$PYENV_ACTIVATE"
pip freeze | grep -v 'pkg-resources' | grep -v 'fenics-dolfin' | grep -v 'mshr' | tee "$CONF_PYPKS"
deactivate
