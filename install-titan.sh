#! /bin/bash

## check for MPI installation
MPCC=
MPCC=`which mpicc`  
MPCPP=
MPCPP=`which mpiCC` || MPCPP=`which mpicxx` 
mpi_home=
if [ "x$MPCPP" == x ]; then
    echo '-----------------------------------------------------------'
    echo ' MPI installtion not found. MPI is required to build titan'
    echo '-----------------------------------------------------------'
    echo 'Enter absolute path to MPI installtion directory if installed in non-standard place.'
    echo '[ctrl-c to quit]'
    read mpi_home
fi

## check for HDF5
header='include/hdf5.h'
hdf5_dir=

#look for hdf5 in some usual places
if [ -f /usr/$header ]; then
    hdf5_dir=/usr
elif [ -f /usr/local/$header ]; then
    hdf5_dir=/usr/local
# ${HOME}
elif [ -f ${HOME}/$header ]; then
    hdf5_dir=${HOME}
# /opt
elif [ -f /opt/$header ]; then
    hdf5_dir=/opt
fi

if [[ "x$hdf5_dir" == x ]]; then
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+  Unable to locate HDF5 installtion. It might not be installed  +"
    echo "+  or installed in unusual path. PARAVIEW support needs hdf5     +"
    echo "+  Please provide absolute path to hdf5 installtion.             +"
    echo "+  HDF5 is freely downloadle from http://hdf.ncsa.uiuc.edu/HDF5/ +"
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "Continue without HDF5? (yes/no) [no]"
    read response
    if [[ "x$response" != xyes ]]; then
        echo "If HDF5 is installed in non-standard path, "
        echo "  provide /full/path/to/hdf5 ...(ctrl-c to QUIT)"
        read hdf5_dir
        if [ ! -f $hdf5_dir/$header ]; then
            echo "ERROR: Could not find hdf5.h in $hdf5_dir"
            exit 1;
        fi
    else
        hdf5_dir=no
    fi
fi

# run configure
MPI=
if test -n "$MPCPP"; then
   MPI="CC=$MPCC CXX=$MPCPP"
else
   MPI="--with-mpi=$mpi_home"
fi
   
./configure $MPI --with-hdf5=$hdf5_dir
make && make install
