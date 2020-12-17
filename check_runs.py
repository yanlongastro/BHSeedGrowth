import glob
import time
import subprocess
import os

folders = glob.glob1('./', '*_42')
for f in folders:
    #print(f)
    gizmo = f+'/gizmo.out'
    #result = subprocess.run(['ls', '-l', '-h', gizmo], stdout=subprocess.PIPE)
    #s = result.stdout
    #print(s[29:-1].decode('utf-8'))
    print('%-50s'%f, end='  ')
    snaps = glob.glob1(f+'/output/', 'snap*')
    ns = len(snaps)-1
    print('%4d'%(ns), end='  ')

    t_final = os.path.getmtime(f+'/output/'+'snapshot_%03d.hdf5'%ns)
    t_start = os.path.getmtime(f+'/output/'+'snapshot_000.hdf5')
    print(time.ctime(t_final), end='  ')
    t = t_final-t_start
    print('%6.2f d' % (t/86400))
    
