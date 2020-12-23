# FastTrip
(Fast) MPI-accelerated (Tr)iplication Waveform (I)nversion (P)ackage

## Requirements
A Linux environment and a FORTRAN compiler (e.g., gfortran or ifort) is required.
This package uses the CPU parallel strategy to effciently search for the model space. Therefore a MPI implementation is needed:
```
sudo apt install mpich
```

To ensure we can use Makefile:
```
sudo apt install make
```

We can also create a new environment for FastTrip:
```
conda create -n FastTrip python=3.7
conda activate FastTrip
```

Some Python modules are needed:
```
conda install numpy scipy mpi4py
```

## Running FastTrip
First, get a clone of the package:
```
git clone https://github.com/lijiaqi0315/FastTrip.git
cd FastTrip
```
To make sure we are using the correct python environment with all the moduled we have downloaded, we can find the absolute path:
```
conda activate FastTrip
which python

/mnt/d/anaconda/anaconda/envs/FastTrip/bin/python
```
Then we specify this path to the package:
```
vi ./go_fort

#!/bin/sh
Python_Path='/mnt/d/anaconda/anaconda/envs/FastTrip/bin/python'
$Python_Path ./Model.py

:wq
```



