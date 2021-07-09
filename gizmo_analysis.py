"""
Analysis code for gizmo, used in BH accretion project

yanlong@caltech.edu
"""


import numpy as np
import re
import h5py
from scipy import spatial
import pandas as pd
import glob
import os


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

def v_esc(M, R):
    return np.sqrt(G*M/R)


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
    if debug:
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
 
def set_folder_name(M, R, Res, pre='BH_', turb_seed=42):
    power = int(np.log10(M))
    coeff = int(M/10**power)
    folder = '%sM%de%d_R%d_S0_T1_B0.01_Res%d_n2_sol0.5_%d/'%(pre, coeff, power, R, Res, turb_seed)
    return folder

class snapshot:
    def __init__(self, file, showinfo=False):
        self.f =  h5py.File(file, 'r')
        self.gas_number = self.f['Header'].attrs['NumPart_ThisFile'][0]
        self.star_number = self.f['Header'].attrs['NumPart_ThisFile'][4]
        self.bh_number = self.f['Header'].attrs['NumPart_ThisFile'][5]
        self.time = self.f['Header'].attrs['Time']
        self.bh_sink_radius = self.f['Header'].attrs['Fixed_ForceSoftening_Keplerian_Kernel_Extent'][5]/2.8
        
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
    
    def star(self, attr, partial=[]):
        if 'PartType4' in list(self.f.keys()):
            if partial == []:
                return self.f['PartType4'][attr][()]
            else:
                return self.f['PartType4'][attr][()][partial]
        else:
            return np.array([])
        
    def bh(self, attr):
        return self.f['PartType5'][attr][()]
    
    def bh_sorted(self, attr):
        ids = self.f['PartType5']['ParticleIDs'][()]
        target = self.f['PartType5'][attr][()]
        a = target[ids.argsort()]
        return a
    
    def single_bh(self, bhid, attr):
        bhpid_base = min(self.f['PartType5']['ParticleIDs'][()])-1
        bhpid = bhpid_base + bhid
        bhpid = np.where(self.f['PartType5']['ParticleIDs'][()]==bhpid)[0][0]
        return self.f['PartType5'][attr][()][bhpid]
            
    def find_gas_near_bh(self, bhid=1, kneighbor=96, drmax=10086, p_norm=2, center=None):
        if center is None:
            pos_bh = self.single_bh(bhid, 'Coordinates')
        else:
            pos_bh = center
        pos_gas = self.gas('Coordinates')
        kdtree = spatial.cKDTree(pos_gas)
        dist, inds = kdtree.query(pos_bh, k=kneighbor, eps=0, distance_upper_bound=drmax, p=p_norm)
        if kneighbor==1:
            dist = np.array([dist])
            inds = np.array([inds])
        return inds[dist!=np.inf]
    
    def find_star_near_bh(self, bhid, kneighbor=96, drmax=10086, p_norm=2):
        pos_bh = self.single_bh(bhid, 'Coordinates')
        pos_gas = self.star('Coordinates')
        kdtree = spatial.cKDTree(pos_gas)
        dist, inds = kdtree.query(pos_bh, k=kneighbor, eps=0, distance_upper_bound=drmax, p=p_norm)
        if kneighbor==1:
            dist = np.array([dist])
            inds = np.array([inds])
        return inds[dist!=np.inf]
    
    
    
def get_num_snaps(path, snap='snapshot_*.hdf5'):
    fns = glob.glob1(path, snap)
    imax = 0
    tmax = 0.
    for i in range(1, len(fns)):
        t = os.path.getmtime(os.path.join(path, 'snapshot_%03d.hdf5'%i))
        if t>tmax:
            imax += 1
            tmax = t
        else:
            break
            
    return imax
        

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

    
    
    
    
class simulation:
    """
    A series of snapshots in a single simulation.
    """
    def __init__(self, folder):
        self.folder = folder
        self.last = get_num_snaps(folder)
        
    def snapshot(self, i):
        return snapshot(self.folder+'snapshot_%03d.hdf5'%i)
        
    def find_interesting_BHs(self, num=5, sort_by_ratio=True):
        sp0 = snapshot(self.folder+'snapshot_%03d.hdf5'%0)
        sp1 = snapshot(self.folder+'snapshot_%03d.hdf5'%(self.last))
        if sort_by_ratio:
            dm = sp1.bh_sorted('Masses')/sp0.bh_sorted('Masses')
        else:
            dm = sp1.bh_sorted('Masses')-sp0.bh_sorted('Masses')
        return (-dm).argsort()[:num]+1
    
    def find_fastest_growth_snapshot(self, num=5, sort_by_ratio=False):
        ms = []
        for i in range(0, self.last+1):
            sp = snapshot(self.folder+'snapshot_%03d.hdf5'%i)
            ms.append(np.sum(sp.bh('Masses')))
        ms = np.array(ms)
        dm = []
        for i in range(self.last):
            if sort_by_ratio:
                dm.append(ms[i+1]/ms[i])
            else:
                dm.append(ms[i+1]-ms[i])
        dm = np.array(dm)
        return (-dm).argsort()[:num]+1
    
    def find_fastest_growth_bh(self, spid, bhs):
        sp0 = snapshot(self.folder+'snapshot_%03d.hdf5'%(spid-1))
        sp = snapshot(self.folder+'snapshot_%03d.hdf5'%spid)
        dbh = []
        for bh in bhs:
            dbh.append(sp.single_bh(bh, 'Masses') - sp0.single_bh(bh, 'Masses'))
        return bhs[np.argmax(dbh)]
    
    def get_bh_history(self, bhid=None, attr='Masses', method=''):
        age = []
        history = []
        for i in range(self.last+1):
            sp = snapshot(self.folder+'snapshot_%03d.hdf5'%i)
            age.append(sp.time)
            if bhid is not None:
                temp = sp.single_bh(bhid, attr)
            else:
                temp = sp.bh(attr)
            
            if 'norm' in method:
                temp = np.linalg.norm(temp, axis=1)
            if 'log' in method:
                temp = np.log10(temp)
            
            if 'sum' in method:
                temp = np.sum(temp)
            if 'inverse_sum' in method:
                temp = np.sum(1/temp)
            if 'average' in method:
                temp = np.mean(temp)
            if 'stdev' in method:
                temp = np.std(temp)    
            
            history.append(temp)
        return np.array(age), np.array(history)
    
    def get_sf_history(self, attr='Masses'):
        age = []
        history = []
        for i in range(self.last+1):
            sp = snapshot(self.folder+'snapshot_%03d.hdf5'%i)
            age.append(sp.time)
            history.append(np.sum(sp.star(attr)))
        return np.array(age), np.array(history)
    
    def get_gas_history(self, attr='Masses', method='sum', crit_density=None):
        age = []
        history = []
        for i in range(self.last+1):
            sp = snapshot(self.folder+'snapshot_%03d.hdf5'%i)
            age.append(sp.time)
            if 'Volume' in attr:
                temp = sp.gas('Masses')/sp.gas('Density')
            elif 'SoundSpeed' in attr:
                temp = np.sqrt(10./9.*sp.gas('InternalEnergy'))
            else:
                temp = sp.gas(attr)
                
            if crit_density is not None:
                temp = temp[sp.gas('Density')>crit_density]
                
            if 'norm' in method:
                temp = np.linalg.norm(temp, axis=1)
            if 'log' in method:
                temp = np.log10(temp)
            
            if 'sum' in method:
                temp = np.sum(temp)
            if 'inverse_sum' in method:
                temp = np.sum(1/temp)
            if 'average' in method:
                temp = np.mean(temp)
            if 'stdev' in method:
                temp = np.std(temp)
            history.append(temp)
        return np.array(age), np.array(history)