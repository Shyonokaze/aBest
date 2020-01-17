# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:45:43 2017

@author: pyh
"""


def readvasprun():
#????vasprun.xml,???????????????????????????
    import xml.etree.ElementTree as ET 
    import numpy as np
    tree = ET.parse('vasprun.xml')
    root = tree.getroot()   
    basis = True

    lc = np.mat(np.zeros((3,3)))
    
    n=0
    child = root[n]    
    while child in root and basis:
        child = root[n]
        n += 1
        for grandson in child:
            if grandson.tag == 'atoms':
                natom = int(grandson.text)
            for ggrandson in grandson:
                if 'NSW' in ggrandson.attrib.values():
                    NSW = int(ggrandson.text)
                if 'POTIM' in ggrandson.attrib.values():
                    POTIM = float(ggrandson.text)
                if 'basis' in ggrandson.attrib.values():
                    for l in range(3):
                        lc[l,:] = np.mat(ggrandson[l].text) 
                    basis = False
    atom=0
    step=0
    
    posa=[None]*NSW
    for i in range(NSW):
        posa[i]=np.mat(np.zeros((natom,3)))
    
    for child in root:
        for grandson in child:
            for ggrandson in grandson:
                if 'positions' in ggrandson.attrib.values():
                    for atom in range(natom):
                        posa[step][atom,:]=np.mat(ggrandson[atom].text)
                    step += 1
    return POTIM,NSW,natom,lc,posa

def readpos(file): #??????selective dynamiscs???????
#??POSCAR/CONTCAR???????????
    import numpy as np
    lc=np.mat(np.zeros((3,3)));
    fid = open(file,'rt')
    line = fid.readline()
    line = fid.readline()
    for i in range(0,3):
        lc[i,:]=np.mat(fid.readline())
    line = fid.readline()
    natom=np.mat(fid.readline())
    number=natom.sum()
    ap=np.mat(np.zeros((number,3)))
    avel=np.mat(np.zeros((number,3)))
    while not 'Direct' in line:
        line = str(fid.readline())
    for i in range(0,number):  
        line = fid.readline()
        if 'T' in line:
            line=line.replace('T','');
        if 'F' in line:
            line=line.replace('F','')
        ap[i,:]=np.mat(line)
    line = fid.readline()
    for i in range(0,number):  
        line = fid.readline()
        avel[i,:]=np.mat(line)
    fid.close()
    return ap,avel

def writepos(file_input,file_output,pos,vel):
#??????????(???)??POSCAR?????,?????POSCAR????
    fido = open(file_input,'rt')
    fidn = open(file_output,'wt')
    while True:
        line=fido.readline()
        print(line[0:(len(line)-1)],file=fidn)
        if 'Direct' in line:
            break
    fido.close()
    for i in range(len(pos)):
        for j in range(3):
            print('%.16f' % pos[i,j],file=fidn,end=' ') 
        print('',file=fidn)
    print('',file=fidn)
 
    for i in range(len(pos)):
        for j in range(3):
            print('%.12e' % vel[i,j],file=fidn,end=' ') 
        print('',file=fidn)
            

def filemv(incar):
#????????vasp??
    import os
    import shutil
    fid = open('record','rt')
    line=fid.readlines()
    ind=line[-1].replace('\n','')
#    ind=line[-1]
    try:
        os.mkdir(ind[:])
        shutil.copy('./POSCAR','./'+ind[:]+'/POSCAR')
        os.symlink('../'+incar,'./'+ind[:]+'/INCAR')
        os.symlink('../KPOINTS','./'+ind[:]+'/KPOINTS')
        os.symlink('../POTCAR','./'+ind[:]+'/POTCAR')
        if os.path.exists('vdw_kernel.bindat'):
            os.symlink('../vdw_kernel.bindat','./'+ind[:]+'/vdw_kernel.bindat')
    except:
        pass
    fid.close()
    
def nearby(lc,pos1,pos2): 
#????????????????????????????
    import numpy as np
    new_pos=pos2
    dp1=pos1*lc
    dp2=pos2*lc
    dis=(dp1-dp2)*np.transpose(dp1-dp2)
    old_dis=0
    n_pos=new_pos.copy()
    while np.all(old_dis != dis):
        old_dis=dis
        for i in range(0,3):
            for j in range(0,3):
                for k in range(0,3):
                    dp2=(new_pos+np.mat([i-1,j-1,k-1]))*lc
                    n_dis=(dp1-dp2)*np.transpose(dp1-dp2)
                    if np.all(n_dis < dis):
                        n_pos=new_pos+np.mat([i-1,j-1,k-1])
                        dis=n_dis
        new_pos=n_pos.copy()
    return new_pos
          

def bond_length(ap,lc,arg1,arg2):
#????????
    import numpy as np
    import vasp
    import math
    arg1 -=1
    arg2 -=1
    ap[arg2,:]=vasp.nearby(lc,ap[arg1,:],ap[arg2,:])
    dp1=ap[arg1,:]*lc
    dp2=ap[arg2,:]*lc
    dis=(dp1-dp2)*np.transpose(dp1-dp2)
    dis=math.sqrt(dis)
    return dis   

def bond_angle(ap,lc,arg1,arg2,arg3):
#????????
    import numpy as np
    import vasp
    import math
    arg1 -=1
    arg2 -=1
    arg3 -=1
    ap[arg1,:]=vasp.nearby(lc,ap[arg2,:],ap[arg1,:])
    ap[arg3,:]=vasp.nearby(lc,ap[arg2,:],ap[arg3,:])
    dp1=ap[arg1,:]*lc
    dp2=ap[arg2,:]*lc
    dp3=ap[arg3,:]*lc    
    ll1=dp2-dp1
    ll2=dp2-dp3
    cita=math.acos(ll1*np.transpose(ll2)/math.sqrt((ll1*np.transpose(ll1))*(ll2*np.transpose(ll2))))*180/math.pi
    return cita

def velcal(lc,posa,potim,num=None,bound=0.1): 
#????????????????,?????????(??????2~n-1??)
    import numpy as np
    import vasp
    if num == None:
        timestep = 2*potim
        vela = [0]*(len(posa)-2)
        check = False
        for i in range(len(vela)):
            vela[i]=np.mat(np.zeros((len(posa[1]),3)))
        for i in range(len(posa)-2):
#            print(i)  #??
            for j in range(len(posa[1])):
                for k in range(3):
#                    if abs(posa[i+2][j,k]) < 0.05 or abs(posa[i+2][j,k]-1) < 0.05 or abs(posa[i+2][j,k]) or abs(posa[i][j,k])> 1: #????
                    if abs(posa[i+2][j,k]) < bound or abs(posa[i+2][j,k]-1) < bound: #????
                        check = True
                        break
                if check:
                    posa[i+2][j,:]=vasp.nearby(lc,posa[i][j,:],posa[i+2][j,:])
                    check = False
                vela[i][j,:]=(posa[i+2][j,:]-posa[i][j,:])*lc/timestep
        return vela
    else:
        num -=1
        timestep = 2*potim
        vela = [0]*(len(posa)-2)
        check = False
        for i in range(len(vela)):
            vela[i]=np.mat(np.zeros((1,3)))
        for i in range(len(vela)):
#            print(i)  #??
            for k in range(3):
#                if abs(posa[i+2][num,k]) < 0.05 or abs(posa[i+2][num,k]-1) < 0.05 or abs(posa[i+2][num,k]) > 1 or abs(posa[i][num,k]) > 1: #????
                if abs(posa[i+2][num,k]) < bound or abs(posa[i+2][num,k]-1) < bound: #????
                    check = True
                    break
            if check:
                posa[i+2][num,:]=vasp.nearby(lc,posa[i][num,:],posa[i+2][num,:])
                check = False
            vela[i][:]=(posa[i+2][num,:]-posa[i][num,:])*lc/timestep
        return vela

def findlimit(vela,num=None,direct=3):
#??????,????,??????????????????????
    if num == None:
        if direct == 'x' or direct == 'X':
            direct = 1
        elif direct == 'y' or direct == 'Y':
            direct = 2
        elif  direct == 'z' or direct == 'Z':
            direct = 3
        direct -=1
        fvmax=0
        fvmin=0
        vmax=float(vela[0][0,direct])
        vmin=float(vela[0][0,direct])
        for i in range(len(vela)):
            if vmax < float(vela[i][0,direct]):
                vmax = float(vela[i][0,direct])
#                print('vmax:',vmax,float(vela[i][num,direct]))  #??
                fvmax = i
            if vmin > float(vela[i][0,direct]):
                vmin = float(vela[i][0,direct])
#                print('vmin:',vmin,float(vela[i][num,direct]))  #??
                fvmin = i
        return fvmax,fvmin

    else:
        num -=1
        if direct == 'x' or direct == 'X':
            direct = 1
        elif direct == 'y' or direct == 'Y':
            direct = 2
        elif  direct == 'z' or direct == 'Z':
            direct = 3
        direct -=1
        fvmax=0
        fvmin=0
        vmax=float(vela[0][num,direct])
        vmin=float(vela[0][num,direct])
        for i in range(len(vela)):
            if vmax < float(vela[i][num,direct]):
                vmax = float(vela[i][num,direct])
#                print('vmax:',vmax,float(vela[i][num,direct]))  #??
                fvmax = i
            if vmin > float(vela[i][num,direct]):
                vmin = float(vela[i][num,direct])
#                print('vmin:',vmin,float(vela[i][num,direct]))  #??
                fvmin = i
        return fvmax,fvmin
    
def checkaway(lc,posa,vela,num,oneatom=True,dis_lim=5,v_d_lim=0.10):#????
#???????????
    import numpy as np
    import math
    import vasp
    num -= 1
    nstep = len(posa)
    check = False
    if oneatom:
        for n in range(7):
            step=int((n+1)*nstep/8)
            posa[step][num,:]=vasp.nearby(lc,posa[0][num,:],posa[step][num,:]) 
            dp=(posa[step][num,:]-posa[0][num,:])*lc
            dis = math.sqrt(dp*np.transpose(dp))
            disv=dp*np.transpose(vela[step])
            if dis > dis_lim:
                check = True
                break
            elif disv > v_d_lim:
                check = True
                break       
        return check        
    else:
        for n in range(7):
            step=int((n+1)*nstep/8)
            posa[step][num,:]=vasp.nearby(lc,posa[0][num,:],posa[step][num,:])
            dp=(posa[step][num,:]-posa[0][num,:])*lc
            dis = math.sqrt(dp*np.transpose(dp))
            disv=dp*np.transpose(vela[step][num,:])
            if dis > dis_lim:
                check = True
                break
            elif disv > v_d_lim:
                check = True
                break       
        return check


def bisec(check):
    fid = open('record','rt')
    line=fid.readlines()
    ran=[float(line[-2].split()[i]) for i in range(len(line[-2].split()))]
    ind=float(line[-1])
    fid.close()
    fidr = open('record','a+')
    if '\n' not in line[-1]:
        print('',file=fidr)
    if abs(abs(ran[1]-ind)-0.001)<=1e-6:
        print(check,file=fidr)
        print('end',file=fidr)
        return 0
    elif check:
        print('%.3f %.3f'% (ran[0],ind),file=fidr)
        print('%.3f'% (float(int((ind+ran[0])*500))/1000),file=fidr)
        return float(int((ind+ran[0])*500))/1000-ind
    elif abs(ind-ran[1])<=1e-6:
        if abs(ind-0.5)<=1e-6:
            print('end',file=fidr)
            return 0
        else:
            print('%.3f %.3f'% (ran[0]+0.1,ran[1]+0.1),file=fidr)
            print('%.3f'% (ind+0.1),file=fidr)
            return 0.1
    else:
        print('%.3f %.3f'% (ind,ran[1]),file=fidr)
        print('%.3f'% (float(int((ind+ran[1])*500))/1000),file=fidr)
        return float(int((ind+ran[1])*500))/1000-ind
    fidr.close()   

if __name__=='__main__':
#    vasprun=readvasprun()
#    potim=vasprun[0]
#    nsw=vasprun[1]
#    natom=vasprun[2]
#    lc=vasprun[3]
#    posa=vasprun[4]
#    vela=velcal(lc,posa,potim)
#    filemv()
    a=bisec(False)   
  
