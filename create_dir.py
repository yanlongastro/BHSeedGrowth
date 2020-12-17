import os
import glob
import shutil
import sys
from distutils.dir_util import copy_tree
import gizmo_analysis as ga
import re
import numpy as np


def replace_lines(lines, key, value, comment='\t %for BH Seed Growth Run\t'):
    key_ = key+'\\b'
    for line in lines:
        if re.search(key_, line):
            print(line, end='')
            i = lines.index(line)
            line = key+'\t'+value+comment+'\n'
            lines[i] = line
            print(line, end='\n')
    return lines


if '.hdf5' not in sys.argv[1]:
    print('Not a hdf5 file!')
    exit()

ics = glob.glob1(os.getcwd(), sys.argv[1])
#print(ics)
path = os.getcwd()
for ic in ics:
    print(ic, '\n')
    newpath = os.path.join(path, ic[:-5])
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)
    if not os.path.exists('output'):
        os.makedirs('output')
    if not os.path.exists('spcool_tables'):
        os.makedirs('spcool_tables')
    os.chdir(path)
    #shutil.copy2('submit.sh', newpath)
    
    
    with open('submit.sh', 'r') as f:
        lines = f.readlines()
    
    job_name = ga.set_job_name(ic)
    lines = replace_lines(lines, r'#SBATCH -J', job_name, comment='')

    with open(os.path.join(path, ic[:-5], 'submit.sh'), 'w+') as f:
        f.writelines("%s" % l for l in lines)
    
    lines = replace_lines(lines, r'mpiexec', '-np 72 ../gizmo/GIZMO ./params.txt 1 1>gizmo.out 2>gizmo.err', comment='') 
    
    with open(os.path.join(path, ic[:-5], 'resubmit.sh'), 'w+') as f:
        f.writelines("%s" % l for l in lines)

    #shutil.copy2(os.path.join(path, ic[:-5], 'submit.sh'), os.path.join(path, ic[:-5], 'resubmit.sh'))
    shutil.copy2('TREECOOL', newpath)
    copy_tree(os.path.join(path, 'spcool_tables'), os.path.join(newpath, 'spcool_tables'))
    
    
    # params.txt file
    par = ga.par_path(ic, skip=0)
    M = par[0]/1e10
    R = par[1]/1e3
    Res = int(par[5])
    t_ff = ga.t_ff(M, R)
    r_sink_star = ga.set_star_softening(M*1e10, Res)
    r_sink_bh = ga.set_bh_softening(M*1e10, R*1e3, Res)
    
    with open('params_template.txt', 'r') as f:
        lines = f.readlines()
        
    lines = replace_lines(lines, r'InitCondFile', '../'+ic[:-5] )
    lines = replace_lines(lines, r'TimeMax', '%.12e'%(t_ff*5) )
    lines = replace_lines(lines, r'TimeBetSnapshot', '%.12e'%(t_ff*5/100) )
    #lines = replace_lines(lines, r'TimeBetStatistics', '%.12e'%(t_ff*5/100) )
    lines = replace_lines(lines, r'Softening_Type4', '%.12e'%(r_sink_star/1e3) )
    lines = replace_lines(lines, r'Softening_Type5', '%.12e'%(r_sink_bh/1e3) )

    with open(os.path.join(path, ic[:-5], 'params.txt'), 'w+') as f:
        #f.write(''.join(lines))
        f.writelines("%s" % l for l in lines)
