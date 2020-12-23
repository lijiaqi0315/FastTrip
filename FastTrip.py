import os
import numpy as np
from scipy import interpolate
from mpi4py import MPI
from scipy import signal
import sys
import argparse




def cal_phi(Station_Num, Length,lowf,highf,delta, folder_index):


    Tstar_1s   =[
      0.0615,   0.3089,    0.5303,    0.5571,    0.4776,    0.3750,    0.2862,    0.2188,    0.1686,    0.1323,    0.1052,    0.0852,
      0.0699,   0.0583,    0.0491,    0.0418,    0.0360,    0.0313,    0.0274,    0.0242,    0.0214,    0.0192,    0.0172,    0.0155,
      0.0141,   0.0128,    0.0117,    0.0108,    0.0099,    0.0092,    0.0085,    0.0079,    0.0074,    0.0069,    0.0064,    0.0060,
      0.0057,   0.0053,    0.0050,    0.0048,    0.0045,    0.0043,    0.0041,    0.0039,    0.0037,    0.0035,    0.0033,    0.0032,
      0.0030,   0.0029,    0.0028,    0.0027,    0.0026,    0.0025,    0.0024,    0.0023,    0.0022,    0.0021,    0.0021,    0.0020,
      0.0019,   0.0019,    0.0018,    0.0017,    0.0017,    0.0016,    0.0016,    0.0015,    0.0015,    0.0014,    0.0014,    0.0014,
      0.0013,   0.0013,    0.0013,    0.0012,    0.0012,    0.0012,    0.0011,    0.0011,    0.0011,    0.0011,    0.0010,    0.0010,
      0.0010,   0.0010,    0.0009,    0.0009,    0.0009,    0.0009,    0.0009,    0.0008,    0.0008,    0.0008,    0.0008,    0.0008,
      0.0008,   0.0008,    0.0007,    0.0007,    0.0007,    0.0007,    0.0007,    0.0007,    0.0007,    0.0006,    0.0006,    0.0006,
      0.0006,   0.0006,    0.0006,    0.0006,    0.0006,    0.0006,    0.0006,    0.0006,    0.0005,    0.0005,    0.0005,    0.0005,
      0.0005,   0.0005,    0.0005,    0.0005,    0.0005,    0.0005,    0.0005,    0.0005,    0.0005,    0.0005,    0.0004,    0.0004,
      0.0004,   0.0004,    0.0004,    0.0004,    0.0004,    0.0004,    0.0004,    0.0004,    0.0004,    0.0004,    0.0004,    0.0004,
      0.0004,   0.0004,    0.0004,    0.0004,    0.0004,    0.0004,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,
      0.0003,   0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,
      0.0003,   0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,    0.0003,
      0.0003,   0.0003,    0.0003,    0.0003,    0.0003,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,
      0.0002,   0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,    0.0002,
      0.0036
    ] 
 
 
 
    ##File Names 
    absolute_path = os.getcwd()
    folder_name = folder_index + '_folder'
    syn_Z_name = absolute_path + '/' + folder_name + '/syn_wave.tz'
    phi_name   = absolute_path + '/' + folder_name + '/phi.dat'
    infor_name = absolute_path + '/' + folder_name + '/infor.dat'


    
    ##Read Windows and Weights
    Station_file = open('Window.txt', 'r')
    Station_name    = [None for i in range(Station_Num)]
    Station_begin1  = [0    for i in range(Station_Num)]
    Station_end1    = [0    for i in range(Station_Num)]
    Station_weight1 = [0.0  for i in range(Station_Num)]
    Station_begin2  = [0    for i in range(Station_Num)]
    Station_end2    = [0    for i in range(Station_Num)]
    Station_weight2 = [0.0  for i in range(Station_Num)]
    Station_begin3  = [0    for i in range(Station_Num)]
    Station_end3    = [0    for i in range(Station_Num)]
    Station_weight3 = [0.0  for i in range(Station_Num)]
    for i in range(Station_Num):
        Station_name[i]    = Station_file.readline()[:-1]
        line1=Station_file.readline()[:-1].split()
        Station_begin1[i]  = int(line1[0])
        Station_end1[i]    = int(line1[1])
        Station_weight1[i] = float(line1[2])
        line2=Station_file.readline()[:-1].split()
        Station_begin2[i]  = int(line2[0])
        Station_end2[i]    = int(line2[1])
        Station_weight2[i] = float(line2[2])
        line3=Station_file.readline()[:-1].split()
        Station_begin3[i]  = int(line3[0])
        Station_end3[i]    = int(line3[1])
        Station_weight3[i] = float(line3[2])
    Station_file.close()


    #Read Obs Data
    Obs_Data = np.zeros([Station_Num,Length])
    Syn_Data = np.zeros([Station_Num,Length])
    Obs_name = [None for i in range(Station_Num)]
    for i in range(Station_Num):
        Obs_name[i]=absolute_path + '/data/' + Station_name[i]
        Obs_txt = np.loadtxt(Obs_name[i])
        Obs_Data[i, :] = Obs_txt[:].copy()

    


    ##Read Syn Data
    data = np.loadtxt(absolute_path + '/' + folder_name + '/syn_seis.tz')
    data_reverse=data.T.copy()

    ##Filter and Attenuation
    Station=np.zeros([Station_Num+1,Length])
    data_filter=np.zeros([Station_Num+1,Length])
    for Station_index in np.arange(0,Station_Num+1,1):
        Station[Station_index,:]=data_reverse[Station_index,:].copy()

    for Station_index in np.arange(1,Station_Num+1,1):
      b, a = signal.butter(1, [2*delta*lowf, 2*delta*highf], 'bandpass')
      tempt_in=Station[Station_index,:].copy()
      tempt_out=signal.filtfilt(b, a, tempt_in)
      Tempt_convolve=np.convolve(tempt_out.copy(),Tstar_1s[:],'full')*delta
      data_filter[Station_index,:]=Tempt_convolve[0:Length].copy()
        


    ##Align the Reference Trace and Find the Max Value for it
    norm_index=int(15-1)   
    res = np.zeros([Station_Num])
    index_align = [0    for i in range(Station_Num)]
    norm_2 = np.zeros([Station_Num])
    norm_3 = np.zeros([Station_Num])
    value_align = 9999
    index_align[norm_index] = 0
    for syn_shift in range(-60,20,1):
        Obs_Data_1 = Obs_Data[norm_index][Station_begin1[norm_index]:Station_end1[norm_index]+1].copy()/max(abs(Obs_Data[norm_index][100:160].copy()))
        Syn_Data_1 = data_filter[norm_index+1,Station_begin1[norm_index]+syn_shift:Station_end1[norm_index]+1+syn_shift].copy()/max(abs(data_filter[norm_index+1,100:160].copy()))
        norm_shift = np.linalg.norm((Obs_Data_1-Syn_Data_1), ord=2)/ np.linalg.norm((Obs_Data_1+Syn_Data_1)/2, ord=2)
        if norm_shift < value_align:
            value_align = norm_shift
            index_align[norm_index] = syn_shift    
    Syn_max=max(map(abs,data_filter[norm_index+1,100+index_align[norm_index]:160+index_align[norm_index]].copy()))
   


    ##Array Normalize All the Traces According to the Reference Trace
    for Station_index in np.arange(1,Station_Num+1,1):
    		data_filter[Station_index,:] = data_filter[Station_index,:].copy()/Syn_max/4
    
    ##Write the Processed Syn Data
    data_filter[0,:]=Station[0,:].copy()
    data_output=data_filter.T.copy()
    np.savetxt(absolute_path + '/' + folder_name + '/syn_wave.tz', data_output, fmt='%.4e')


    ##Read the Processed Syn Data
    Syn_txt = np.loadtxt(syn_Z_name)
    for i in range(Station_Num):
        Syn_Data[i, :] = Syn_txt[:, i + 1].copy()

    ##Misfit Calculation
    for i in range(Station_Num):
        value_align = 9999
        res[i]=0
        index_align[i] = 0
        ##Find the Time Shift for Alignment Using Window1
        for syn_shift in range(-60,20,1):
            Obs_Data_1 = Obs_Data[i][Station_begin1[i]:Station_end1[i]+1].copy()
            Syn_Data_1 = Syn_Data[i][Station_begin1[i]+syn_shift:Station_end1[i]+1+syn_shift].copy()
            norm_shift = np.linalg.norm((Obs_Data_1-Syn_Data_1), ord=2)/ np.linalg.norm((Obs_Data_1+Syn_Data_1)/2, ord=2)
            if norm_shift < value_align:
                value_align = norm_shift
                index_align[i] = syn_shift

        ##After Shifting, Calculate the Misfit for Window2 and Window3
        Obs_Data_2 = Obs_Data[i][Station_begin2[i]:Station_end2[i]+1].copy()
        Syn_Data_2 = Syn_Data[i][Station_begin2[i]+index_align[i]:Station_end2[i]+1+index_align[i]].copy()
        Obs_Data_3 = Obs_Data[i][Station_begin3[i]:Station_end3[i]+1].copy()
        Syn_Data_3 = Syn_Data[i][Station_begin3[i]+index_align[i]:Station_end3[i]+1+index_align[i]].copy()

        ##Weight the Windows
        norm_2[i] = np.linalg.norm((Obs_Data_2-Syn_Data_2), ord=2) / np.linalg.norm((Obs_Data_2+Syn_Data_2)/2, ord=2) / Station_weight2[i]
        norm_3[i] = np.linalg.norm((Obs_Data_3-Syn_Data_3), ord=2) / np.linalg.norm((Obs_Data_3+Syn_Data_3)/2, ord=2) / Station_weight3[i]
 
        ##Store the Misfit
        res[i] = res[i] + norm_2[i]
        res[i] = res[i] + norm_3[i]


    ##Write the Misfit
    residual=0
    for i in range(Station_Num):
        residual=residual+res[i]
    residual=residual/Station_Num

    phi_file=open(phi_name,'w')
    phi_file.write('%f\n' %(residual))
    phi_file.close()


    ##Write Information
    infor_file=open(infor_name,'w')
    for i in range(Station_Num):
        infor_file.write('Station Num: %2d   shift: %3d\n' %(i,index_align[i]))
        infor_file.write('res : %f  weight: %4.2f\n' % (norm_2[i], Station_weight2[i]))
        infor_file.write('res : %f  weight: %4.2f\n' % (norm_3[i], Station_weight3[i]))
        infor_file.write('%f\n' % (res[i]))
    infor_file.write('\n')
    infor_file.write('Final Residual :\n')
    infor_file.write('%f\n' % (residual))
    infor_file.close()

def get_args(args=None):
    parser = argparse.ArgumentParser(
        description='pre-proc obs and syn data (depth pert and/or cmt3d pert)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--loop_core',
                        help='loop core',
                        required=True)

    results = parser.parse_args(args)
    return results.loop_core


def main():
    
    loop_core= get_args(sys.argv[1:])
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    arg = int(rank) + int(loop_core)*int(size)
    absolute_path=os.getcwd()

    folder_index = format("%03d" % (arg + 000))
    folder_name= folder_index + '_folder'

    ##Run Simulation
    os.system(absolute_path + '/' + folder_name + '/qseis ' + folder_index)  

    ##Parameters
    Station_Num=17
    Length=401
    delta=0.25
    lowf=0.05
    highf=1

    ##Calculate Misfit
    cal_phi(Station_Num,Length,lowf,highf,delta,folder_index)

main()
