
def get_src(location):
    if location == 'home':
        mpirun = '/Users/krishna/Programs/build/openmpi/bin/mpirun'
        twopoint = '/Users/krishna/Research/github/twofast/TWOPOINT'
        twopoint_mpi = '/Users/krishna/Research/github/twofast/TWOPOINT_MPI'
    elif location == 'splinter':
        mpirun = '/opt/openmpi/bin/mpirun'
        twopoint = '/share/data1/knaidoo/splinter_libraries/src/twofast/TWOPOINT'
        twopoint_mpi = '/share/data1/knaidoo/splinter_libraries/src/twofast/TWOPOINT_MPI'
    return mpirun, twopoint, twopoint_mpi
