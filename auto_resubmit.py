"""
Auto resubmission to avoid Wheeler hanging.

To run the process:
nohup python auto_resubmit.py &

To kill the process:
kill -9 $(ps aux | grep python | grep auto_resubmit.py | awk '{print $2}')
"""


import subprocess
import re
import os
import time
import numpy as np
import glob
import gizmo_analysis as ga


MRs = [
[1e4, 5, 128],
[1e5, 5, 128],
[1e6, 5, 128],
[1e6, 50, 128],
[1e7, 50, 128],
[1e8, 50, 128],
[1e8, 500, 128],
[1e9, 500, 128],
[1e10, 500, 128],
#[1e6, 20, 64],
[1e6, 20, 128],
[1e10, 500, 256],

[1e4, 5, 64],
[1e6, 50, 64],
[1e8, 500, 64],

[1e5, 5, 64],
[1e7, 50, 64],
[1e9, 500, 64],

[1e6, 5, 64],
[1e8, 50, 64],
[1e10, 500, 64],

[1e8, 50, 256],
[3e7, 50, 128],
]

timeout = 2400

while 1:
    finished = True
    print(time.ctime(), '\n')
    MRs = np.loadtxt('auto_resubmit_params.txt')
    for MR in MRs:
        M = MR[0]
        R = MR[1]
        Res = MR[2]
        turb_seed = MR[3]
        folder = ga.set_folder_name(M, R, Res, turb_seed=turb_seed)
        t_gizmo = os.path.getctime(folder+'gizmo.out')

        slurms = glob.glob1(folder, 'slurm*')
        t_slurm = 0.0
        for slurm in slurms:
            if os.path.getctime(folder+slurm) > t_slurm:
                t_slurm = os.path.getctime(folder+slurm)
                job_id = slurm[-10:-4]
        #print(job_id)
        
        t_latest = 0
        for i in range(len(glob.glob1(folder+'output/', 'snap*'))):
            temp = os.path.getctime(folder+'output/snapshot_%03d.hdf5'%i)
            if temp > t_latest:
                n_snapshot = i
                t_latest = temp
        #print(n_snapshot)
        #if n_snapshot < 100:
        #    finished = False

        if R==20:
            n_snap_tot = 500
        else:
            n_snap_tot = 100

        if n_snapshot < n_snap_tot:
            finished = False
        print(folder)

        job_name = ga.set_job_name(folder)
        cmd = "squeue -u yanlong |grep %s"%job_name
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout.decode('UTF-8')
        try:
            res_ = res.split()
            job_id_submitted = res_[0]
            status = res_[4]
        except:
            job_id_submitted = ''
            status = ''
            pass
        # print(job_name, job_id, job_id_submitted, status)


        if 'PD' in status:
            print('...pending')
            continue
        if 'ReqNodeNotAvail' in res:
            subprocess.run(["scancel", job_id])
            print('...broken; cancelled submission:', job_id)

        if time.time() - t_gizmo > timeout and n_snapshot < n_snap_tot:
            res = subprocess.run(['squeue', '-j', job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            squeue = (res.stdout).decode('UTF-8')
            
            # If hanging or stopped, resubmit
            if job_id in squeue:
                print('...hanging')
                subprocess.run(["scancel", job_id])
            else:
                print('...stopped')
            time.sleep(5)
            os.chdir(folder)
            res = subprocess.run(["sbatch", "resubmit.sh"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            newjobid = (res.stdout).decode('UTF-8')[-7:-1]

            # If broken, cancel
            res = subprocess.run(['squeue', '-j', newjobid], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            squeue = (res.stdout).decode('UTF-8')
            if 'ReqNodeNotAvail' in squeue:
                subprocess.run(["scancel", newjobid])
                print('...cancelled submission:', newjobid)
            else:
                print('...submitted:', newjobid)
            os.chdir('../')
        elif n_snapshot >= n_snap_tot:
            print('...completed')
        else:
            print('...running')
    
    time.sleep(600)
    if finished == True:
        break
    print('\n')
