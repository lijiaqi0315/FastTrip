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
conda install numpy scipy
```

Install mpi4py through pip install:
```
pip install mpi4py
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
./Optimization_NGA
```
(The 'Optimization_NGA' is a binary file, which is pre-compiled on a 64-bit machine. If you fail to run it, you can also type 'make' in the source code folder './src_Optimization_NGA' to compile it on your own machine.)


Batch deletion for the output files is also supported:
```
make rm
```


# Tips
## 1. Focal Mechanism
Note that the definitions of the coordinates in Qseis (Mxyz) and GCMT (Mtpr) are different. Therefore the moment tensor downloaded from GCMT needs to be converted:
```
Mxx=+Mtt
Myy=+Mpp
Mzz=+Mrr
Mxy=-Mtp
Myz=-Mrp
Mzx=+Mrt
```

## 2. Normalization
If you prefer use the array normalization (default in the code), this means that you trust the absolute amplitudes of your stations. If, somtimes, the absolute amplitudes maybe not that reliable, trace normalization works better. You can swift to trace normalization in this way:

Original codes:
```
##Array Normalize All the Traces According to the Reference Trace
for Station_index in np.arange(1,Station_Num+1,1):  
    data_filter[Station_index,:] = data_filter[Station_index,:].copy()/Syn_max/4
```

Trace-normalization codes:
```
##Trace Normalize All the Traces According to the Reference Trace
for Station_index in np.arange(1,Station_Num+1,1):  
    Trace_max=max(map(abs,data_filter[Station_index,:].copy()))                                        
    data_filter[Station_index,:] = data_filter[Station_index,:].copy()/Trace_max/4
```

Note that here I choose a normalization factor of 4 for the synthetic waveforms (e.g, /Trace_max/4). In this way, you need to also divide your data by 4.







## Acknowledgments

* We thank Yanbin Wang (Peking University) for testing the FastTrip package and Zhigang Peng (Georgia Tech) for valuable advice.
* We appreciate Rongjiang Wang (GFZ) for providing the source code of QSEIS (Wang, 1999).
* We thank the Institute for Cyber-Enabled Research (ICER) at Michigan State University, the Extreme Science and Engineering Discovery Environment (XSEDE supported by NSF grant ACI-1053575), and the High-performance Computing Platform of Peking University for providing the high-performance computing resources.
* This research was supported by NSF grant 1802247 and the startup fund of Min Chen at Michigan State University.
