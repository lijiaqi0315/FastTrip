import os
import numpy as np
from scipy import interpolate
from scipy import signal



def main():

    ##Interp Rho
    depth_iasp=[
          0,         3,         3,         3,         3,         10,        10,        18,        18,        43,         80,       80,       120,       165,       210,
        210,       260,       310,       360,       410,        410,       460,       510,       560,       610,        660,      660,       710,       760,       809,
        859,       908,       958,      1007,      1057,       1106,      1156,      1205,      1255,      1304,       1354,     1403,      1453,      1502,      1552,
       1601,      1651,      1700,      1750,      1799,       1849,      1898,      1948,      1997,      2047,       2096,     2146,      2195,      2245,      2294,
       2344,      2393,      2443,      2492,      2542,       2591,      2640,      2690,      2740,      2740,       2789,     2839,      2891
    ]

    ro_iasp   =[
      2.6000,    2.6000,    2.6000,    2.6000,    2.6000,    2.6000,    2.6000,    2.6000,    2.6000,    3.3787,    3.3747,    3.3747,    3.3704,    3.3655,    3.3606,
      3.4298,    3.4597,    3.4895,    3.5194,    3.5492,    3.7364,    3.7994,    3.8624,    3.9254,    3.9781,    3.9898,    4.3745,    4.4056,    4.4365,    4.4667,
      4.4966,    4.5262,    4.5555,    4.5845,    4.6132,    4.6417,    4.6700,    4.6980,    4.7257,    4.7532,    4.7805,    4.8076,    4.8345,    4.8612,    4.8877,
      4.9141,    4.9402,    4.9662,    4.9921,    5.0178,    5.0434,    5.0689,    5.0942,    5.1194,    5.1446,    5.1696,    5.1946,    5.2195,    5.2443,    5.2691,
      5.2938,    5.3185,    5.3431,    5.3678,    5.3924,    5.4170,    5.4412,    5.4661,    5.4910,    5.4910,    5.5158,    5.5406,    9.9042
    ]




    absolute_path=os.getcwd()
    output_file_name='.fortran_output_file'
    input_file_name='.fortran_input_file'

    ##Basic settings
    layers=22
    ka=1.732
    CoreNumber=5
    ##


    ##Python Path
    go_fort = open ( absolute_path + '/go_fort', 'r')
    readpath=[]
    for line in go_fort:
    	    readpath.append(line)   
    go_fort.close()
    Python_Path=readpath[1][12:-1]
    ##


    ##Compile QSEIS
    if not os.path.isfile('qseis'): 
    	##Compile Qseis##
    	qseis_main_backup = open ( absolute_path + '/' + 'qseis_code' + '/qsmain.f_backup', 'r')
    	qseis_main = open(absolute_path + '/' + 'qseis_code' + '/qsmain.f', 'w')
    	new = []
    	for line in qseis_main_backup:
    	    new.append(line)
    	new[39] = '     &\'%s/\'\n'%absolute_path
    	new[48] = '     &\'%s/Qseis.input\',status=\'old\')\n'%absolute_path
    	new[14] = '      character(%d) absolute_path\n' % (len(absolute_path) + 1)
    	new[16] = '      character(%d) absolute_name\n' % (len(absolute_path) + 12)
    	for n in new:
    	    qseis_main.write(n)
    	qseis_main_backup.close()
    	qseis_main.close()
    	
    	qsfftinv_backup = open ( absolute_path + '/' + 'qseis_code' + '/qsfftinv.f_backup', 'r')
    	qsfftinv = open(absolute_path + '/' + 'qseis_code' + '/qsfftinv.f', 'w')
    	new_2 = []
    	for line in qsfftinv_backup:
    	    new_2.append(line)
    	new_2[15] = '      character(%d) absolute_name\n' % (len(absolute_path) + 12)
    	for n in new_2:
    	    qsfftinv.write(n)
    	qsfftinv_backup.close()
    	qsfftinv.close()
    	
    	os.chdir('qseis_code')
    	os.system('make clean')
    	os.system('make')
    	os.system('cp qseis ..')
    	os.chdir('..')
    	os.system('chmod u+x qseis')
    ##

    
    ##Generate Model Folders
    input_file = open(input_file_name, 'r')
    first_line_of_input=input_file.readline()[:-1].split()
    npar=int(first_line_of_input[0])
    nmod=int(first_line_of_input[1])
    ndeme=int(first_line_of_input[2])
    folder_index=[None for i in range(nmod*ndeme)]
    folder_name=[None for i in range(nmod*ndeme)]
    phi_name=[None for i in range(nmod*ndeme)]
    misfit_name=[None for i in range(nmod*ndeme)]

    for j in range(nmod*ndeme):
      folder_index[j] = format("%03d" %(j+000))
      folder_name[j] = folder_index[j] + '_folder'
      phi_name[j] = absolute_path + '/' + folder_name[j] + '/phi.dat'
      misfit_name[j] = absolute_path + '/' + folder_name[j] + '/misfit.dat'
      if not os.path.exists(folder_name[j]):
          os.system('mkdir '+ folder_name[j])
          os.system('cp qseis '+ folder_name[j])

      ##Set Model parameters
      para = [0 for i in range(npar)]
      for k in range(npar):
        para[k] = float(input_file.readline()[:-1])

      depth = np.zeros([layers])
      vp    = np.zeros([layers])
      vs    = np.zeros([layers])
      ro    = np.zeros([layers])
      qp    = np.zeros([layers])
      qs    = np.zeros([layers])
      velocity_file = open('velocity.dat', 'w')
      velocity_file.write('%d\n'%layers)

      depth[0]=0;          depth[1]= 20;         depth[2]= 20;         depth[3]= 35;         depth[4]= 35;
      vp[0]=5.8;           vp[1]=5.8;            vp[2]=6.5;            vp[3]=6.5;            vp[4]=8.04;
      vs[0]=vp[0]/ka;      vs[1]=vp[1]/ka;       vs[2]=vp[2]/ka;       vs[3]=vp[3]/ka;       vs[4]=vp[4]/ka;

      depth[5]=120;        depth[6]=165;         depth[7]=210;         depth[8]=250;         depth[9]=290;
      vp[5]=8.05;          vp[6]=8.175;          vp[7]=8.3;            vp[8]=8.45+para[0];   vp[9]=8.60+para[1];
      vs[5]=vp[5]/ka;      vs[6]=vp[6]/ka;       vs[7]=vp[7]/ka;       vs[8]=vp[8]/ka;       vs[9]=vp[9]/ka;

      depth[10]=330;       depth[11]=370;        depth[12]=410+para[9];depth[13]=411+para[9];depth[14]=450;
      vp[10]=8.74+para[2]; vp[11]=8.88+para[3];  vp[12]=9.03+para[4];  vp[13]=9.36+para[5];  vp[14]=9.50+para[6];
      vs[10]=vp[10]/ka;    vs[11]=vp[11]/ka;     vs[12]=vp[12]/ka;     vs[13]=vp[13]/ka;     vs[14]=vp[14]/ka;

      depth[15]=490;       depth[16]=530;        depth[17]=570;        depth[18]=600;        depth[19]=660;
      vp[15]=9.63+para[7]; vp[16]=9.76+para[8];  vp[17]=9.901;         vp[18]=10.01;         vp[19]=10.24;
      vs[15]=vp[15]/ka;    vs[16]=vp[16]/ka;     vs[17]=vp[17]/ka;     vs[18]=vp[18]/ka;     vs[19]=vp[19]/ka;
      
      depth[20]=710;       depth[21]=800;        
      vp[20]=10.42;        vp[21]=10.53;         
      vs[20]=vp[20]/ka;    vs[21]=vp[21]/ka;         
      

      for kk in range(layers):
        ro_interp = interpolate.interp1d(depth_iasp,ro_iasp,kind="slinear")
        ro[kk] = ro_interp(depth[kk])
        velocity_file.write('%2d %9.3f %9.4f %9.4f %9.4f 9999.999 9999.999\n' %(kk+1,depth[kk],vp[kk],vs[kk],ro[kk]))

      velocity_file.close()
      os.system('cp velocity.dat '+ folder_name[j])

    input_file.close()



    ##Run Simulation and Misfit
    for loop_core in range(0,int(100/CoreNumber)):
        os.system('mpirun -n '+str(CoreNumber)+' '+Python_Path+' FastTrip.py --loop_core '+ str(loop_core))

    ##


    ##Write Misfit
    out_file = open(output_file_name, 'w')
    for j in range(nmod*ndeme):
        phi_file = open(phi_name[j], 'r')
        phi_file_content = phi_file.readline()[:-1]
        out_file.write('%f\n' %(float(phi_file_content)))
        os.system('cat '+phi_name[j]+' >> '+misfit_name[j])
        phi_file.close()
    out_file.close()



main()
