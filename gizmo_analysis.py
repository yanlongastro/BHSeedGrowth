"""
Analysis code for gizmo, used in BH accretion project

yanlong@caltech.edu
"""


import numpy as np
import re
import h5py
from scipy import spatial
import pandas as pd


tunit = 206265*1000*1.5e8/(86400*365)
G = 4*np.pi**2/(206265000**3/1e10/tunit**2)
pc_in_cm = 3.086e+18
Msun_in_g = 2e33
Na = 6.02e23


def t_ff(M, R):
    """
    Free fall time: everything in code unit
    """
    return np.pi/2 *np.sqrt(R**3/G/M/2)


def show_info(fin):
    for k in list(fin.keys()):
        print(k)
        for atrs in list(fin[k].attrs.keys()):
            print(' ', atrs, '=', fin[k].attrs[atrs])
        for sk in list(fin[k].keys()):
            print(' ', end=' ')
            print(sk)
            if k == 'PartType5':
                print('  ', end=' ')
                print(fin[k][sk][()])
            

            
def par_path(path, skip=1):
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    num_list = [float(x) for x in re.findall(match_number, path)]
    return num_list[skip:]

# Unit: Msun, pc, yr
def set_star_softening(M, Res):
    dM = M/Res**3
    drho = 2227.359223512722/Na
    drho /= Msun_in_g /pc_in_cm**3
    dr = (dM/drho)**(1/3)
    # print(dM, drho)
    return dr


# Unit: Msun, pc, yr
def set_bh_softening(M, R, Res, debug=False):
    dr_v = 100/M*R
    if debug:
        print('dr_v =', dr_v)
    dr_star = set_star_softening(M, Res)
    if debug:
        print('dr_star =', dr_star)
    cs2 = 10000*8.31/0.001*5/3
    
    # Cooling: this is arbitrary
    cs2 /= 100
    
    # to km/s
    cs2 /= 1000**2
    dr_cs = G*100/1e10/cs2
    # to pc
    dr_cs *= 1000
    if debug:
        print('dr_cs =', dr_cs)
    dr = min(dr_star, dr_cs)
    if dr == dr_star:
        print('Used star softenning length.')
    if dr == dr_cs:
        print('Used sound speed limit.')
    if dr == dr_v:
        print('Used velocity limit.')
    #print(dr)
    return dr


def check_vel(m_bh, xyz_bh, vxyz_bh, m_gas, pos_gas, vel_gas, den):
    r = np.linalg.norm(xyz_bh-pos_gas)
    v = np.linalg.norm(vxyz_bh-vel_gas)
    m_add = 4*np.pi/3*r**3*den
    return (v**2<2*G*(m_bh+m_gas+m_add)/r)


def check_ang_mom(m_bh, xyz_bh, vxyz_bh, m_gas, pos_gas, vel_gas, sink_radius):
    dr = xyz_bh-pos_gas
    dv = vxyz_bh-vel_gas
    drdv = np.sum(dr*dv)
    r = np.linalg.norm(dr)
    v = np.linalg.norm(dv)
    spec_mom = (r*v)**2 - drdv**2
    return (spec_mom < G*(m_bh+m_gas)*sink_radius)


def check_boundedness(xyz_bh, pos_gas, sink_radius):
    r = np.linalg.norm(xyz_bh-pos_gas)
    return (r<sink_radius)

def set_job_name(ic, skip=0):
    if 'BH' in ic:
        job_name = 'B'
    else:
        job_name = 'M'

    par = par_path(ic, skip=skip)
    M = par[0]
    R = par[1]
    Res = int(par[5])

    ind = int(np.log10(M))
    coe = int(M/10**ind)
    job_name += str(coe)+'%d'%ind

    ind = int(np.log10(R))
    coe = int(R/10**ind)
    job_name += str(coe)+str(ind)

    job_name += '%d'%(np.log(Res)/np.log(2))
    return job_name
 
def set_folder_name(M, R, Res):
    power = int(np.log10(M))
    coeff = int(M/10**power)
    folder = 'BH_M%de%d_R%d_S0_T1_B0.01_Res%d_n2_sol0.5_42/'%(coeff, power, R, Res)
    return folder

class snapshot:
    def __init__(self, file, showinfo=False):
        self.f =  h5py.File(file, 'r')
        self.gas_number = self.f['Header'].attrs['NumPart_ThisFile'][0]
        self.star_number = self.f['Header'].attrs['NumPart_ThisFile'][4]
        self.bh_number = self.f['Header'].attrs['NumPart_ThisFile'][5]
        self.time = self.f['Header'].attrs['Time']
        
        if showinfo == True:
            show_info(self.f)
    
    def close(self):
        self.f.close()
        
    def gas(self, attr, partial=[]):
        res = self.f['PartType0'][attr][()]
        if partial == []:
            return res
        else:
            return res[partial]
    
    def star(self, attr):
        if 'PartType4' in list(self.f.keys()):
            return self.f['PartType4'][attr][()]
        else:
            return None
        
    def bh(self, attr):
        return self.f['PartType5'][attr][()]
    
    def single_bh(self, bhid, attr):
        bhpid_base = min(self.f['PartType5']['ParticleIDs'][()])-1
        bhpid = bhpid_base + bhid
        bhpid = np.where(self.f['PartType5']['ParticleIDs'][()]==bhpid)[0][0]
        return self.f['PartType5'][attr][()][bhpid]
            
    def find_gas_near_bh(self, bhid, kneighbor=96, drmax=10086):
        pos_bh = self.single_bh(bhid, 'Coordinates')
        pos_gas = self.gas('Coordinates')
        kdtree = spatial.cKDTree(pos_gas)
        dist, inds = kdtree.query(pos_bh, k=kneighbor, eps=0, distance_upper_bound=drmax)
        if kneighbor==1:
            dist = np.array([dist])
            inds = np.array([inds])
        return inds[dist!=np.inf]
        

def pass_row_header(fname):
    """
    delete 'BH=' in every line
    """
    with open(fname, 'r') as fin:
        if 'BH=' in fin.readline():
            for line in fin:
                try:
                    yield line[3:]
                except IndexError:
                    continue
        else:
            for line in fin:
                try:
                    yield line
                except IndexError:
                    continue

        
        
class blackhole_details:
    def __init__(self, outputdir, filename='blackhole_details', tasks=72, io_reduced_mode=False):
        self.io_reduced_mode = io_reduced_mode
        df = None
        for task in range(tasks):
            bhdetail = outputdir+'blackhole_details/'+filename+'_%d.txt'%task
            if io_reduced_mode:
                data = np.loadtxt(pass_row_header(bhdetail))
                sort_key = 1
            else:
                data = np.loadtxt(bhdetail)
                sort_key = 0
            pdf = pd.DataFrame(data=data)
            df = pd.concat([df, pdf])
            print(task, end=' ')
        df = df.sort_values(by=[sort_key])
        self.df = df
    def get_detail(self, bhpid, column):
        if self.io_reduced_mode:
            mask_key = 0
        else:
            mask_key = 1
        mask = self.df[mask_key] == bhpid
        res = self.df[mask][column]
        return res.values
