#!/usr/bin/env bash
#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

# Configure libMesh with the default MOOSE configuration options
#
# Separated so that it can be used across all libMesh build scripts:
# - scripts/update_and_rebuild_libmesh.sh
# - conda/libmesh/build.sh
function configure_libmesh()
{
  ../configure --enable-silent-rules \
               --enable-unique-id \
               --disable-warnings \
               --enable-glibcxx-debugging \
               --with-thread-model=openmp \
               --disable-maintainer-mode \
               --enable-petsc-hypre-required \
               --enable-metaphysicl-required \
               --with-cxx-std-min=2014

  return $?
}
