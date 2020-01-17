# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 14:58:01 2017

@author: pyh
"""
import os
import shutil
import sys
sys.path.append('/public/home/users/ruc001/bin/aBest')
import vasp

NP=sys.argv[1]
PBS_NODEFILE=sys.argv[2]
num1=int(sys.argv[3])
nvtstep=1999

run='/public/software/mpi/openmpi/1.6.5/intel/bin/mpirun -np '+NP+' -machinefile '+PBS_NODEFILE+' --mca btl self,sm,openib --bind-to-core  ~/bin/vasp5c  | tee sout'
try:
    os.mkdir('NVEmin')
    os.chdir('NVEmin')
    os.symlink('../INCAR_NVE','./INCAR')
    os.symlink('../KPOINTS','./KPOINTS')
    os.symlink('../POTCAR','./POTCAR')
    shutil.copy('../pos_vmax','./pos_vmax')
    shutil.copy('../record','./record')
    os.chdir('..')
except:
    pass

os.chdir('NVEmin')
ap,avel=vasp.readpos('pos_vmax')
avel[num1-1,2]=avel[num1-1,2]-0.1
vasp.writepos('pos_vmax','POSCAR',ap,avel)

fid = open('record','rt')
line=fid.readlines()
fid.close()

while 'end\n' not in line:
    vasp.filemv('INCAR')
    ind=line[-1].replace('\n','')
    os.chdir(ind)
    shutil.copy('../record','./record')
    os.system(run)
    
    potim,nsw,natom,lc,posa=vasp.readvasprun()
    vela1=vasp.velcal(lc,posa,potim,num1)
    check1=vasp.checkaway(lc,posa,vela1,num1)

    ns=vasp.bisec(check1)

    if ns:
        ap,avel=vasp.readpos('POSCAR')
        avel[num1-1,2]=avel[num1-1,2]-ns
        vasp.writepos('POSCAR','pos_next',ap,avel)
    try:    
        shutil.copy('./pos_next','../POSCAR')
    except:
        pass
    shutil.copy('./record','../record')
    
    fid = open('record','rt')
    line=fid.readlines()
    fid.close()
    
    os.chdir('../')
    
os.chdir('../')

try:
    os.mkdir('NVEmax')
    os.chdir('NVEmax')
    os.symlink('../INCAR_NVE','./INCAR')
    os.symlink('../KPOINTS','./KPOINTS')
    os.symlink('../POTCAR','./POTCAR')
    shutil.copy('../pos_vmin','./pos_vmin')
    shutil.copy('../record','./record')
    os.chdir('..')
except:
    pass

os.chdir('NVEmax')
ap,avel=vasp.readpos('pos_vmin')
avel[num1-1,2]=avel[num1-1,2]-0.1
vasp.writepos('pos_vmin','POSCAR',ap,avel)

fid = open('record','rt')
line=fid.readlines()
fid.close()

while 'end\n' not in line:
    vasp.filemv('INCAR')
    ind=line[-1].replace('\n','')
    os.chdir(ind)
    shutil.copy('../record','./record')
    os.system(run)
    
    potim,nsw,natom,lc,posa=vasp.readvasprun()
    vela1=vasp.velcal(lc,posa,potim,num1)
    check1=vasp.checkaway(lc,posa,vela1,num1)

    ns=vasp.bisec(check1)

    if ns:
        ap,avel=vasp.readpos('POSCAR')
        avel[num1-1,2]=avel[num1-1,2]-ns
        vasp.writepos('POSCAR','pos_next',ap,avel)
    
    try:
        shutil.copy('./pos_next','../POSCAR')
    except:
        pass

    shutil.copy('./record','../record')
    
    fid = open('record','rt')
    line=fid.readlines()
    fid.close()
    
    os.chdir('../')
