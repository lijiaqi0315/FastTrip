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

Line 2: Python_Path='/mnt/d/anaconda/anaconda/envs/FastTrip/bin/python'

:wq
```

According to the status of your computer, specify the number of cores:
```
vi ./Model.py

...
Line 37:CoreNumber=5
...

:wq
```

We can simply run FastTrip for the example in the paper:
```
./Run_FastTrip
```

Batch deletion for the output files is also supported:
```
make rm
```

## Acknowledgments

* We thank Yanbin Wang (Peking University), and Shaohua Li (CEA, LanZhou) for testing the FastTrip package and Zhigang Peng (Georgia Tech) for valuable advice. We appreciate Rongjiang Wang (GFZ) for providing the source code of QSEIS.
* We thank the Institute for Cyber-Enabled Research (ICER) at Michigan State University, the Extreme Science and Engineering Discovery Environment (XSEDE supported by NSF grant ACI-1053575), and the High-performance Computing Platform of Peking University for providing the high-performance computing resources.
* This research was supported by NSF grant 1802247 and the startup fund of Min Chen at Michigan State University.
