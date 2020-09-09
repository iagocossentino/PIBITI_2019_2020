# Simulação de oximetria baseada nos trabalhos de Gerhard Magnus e Joel Keizer:
# 1) Minimal model of b-cell mitochondrial Ca2+ handling (DOI: 10.1152/ajpcell.1997.273.2.C717)
# 2) Model of b-cell mitochondrial calcium handling and electrical activity. I. Cytoplasmic variables (DOI: 10.1152/ajpcell.1998.274.4.C1158)
# Código referente ao projeto de iniciação científica dos estudantes: Amanda dos Santos Pereira, Elysa Beatriz de Oliveira Damas e Iago Cossentino de Andrade
# Orientador: Jair Trapé Goulart 
# Instituição: Universidade de Brasília (UnB)

# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style

import numpy as np
from scipy.integrate import odeint
import tkinter as tk
import tkinter.scrolledtext as st
import os.path 

# Parâmetros tabelados --------------------------------------------------------------------------------------------------------------------------------------------------
outros = {'F':96485,'R':8314,'NADtot':8,'Atot':12,'Atoti':1.5}    
tabela1 = {'deltapH':-0.4,'T':310,'gH':0.2}
tabela2 = {'Kres':1.35e18,'pres':0.4,'r1':2.077e-18,'r2':1.728e-9,'r3':1.059e-26,'ra':6.394e-10,'rb':1.762e-13,'rc1':2.656e-19,'rc2':8.632e-27,'dpB':50,'g':0.85}
tabela3 = {'KF1':1.71e6,'Pim':20,'pF1':0.7,'p1':1.346e-8,'p2':7.739e-7,'p3':6.65e-15,'pa':1.656e-5,'pb':3.373e-7,'pc1':9.651e-14,'pc2':4.845e-19,'dpB':50}
tabela4 = {'JmaxANT':1000,'f':0.5}
tabela5 = {'Kcat':6.0,'Kact':0.38,'Jmaxuni':400,'dpuni':91,'nact':2.8,'Lunimax':50,'KNa':9.4,'KCa':0.003,'Nai':30,'dpNaCa':91,'b':0.5,'JmaxNaCa':5.5,'n':3}
betas = {'betamax':126,'beta1':1.66,'beta2':0.0249,'beta3':4,'beta4':2.83,'beta5':1.3,'beta6':2.66,'beta7':0.16}
table3 = {'u1':15,'u2':1.1,'KCa2':0.05,'Jredbasal':20,'CCa2m':3000,'Cmito':1.45e-3,'khyd':41,'tauhyd':50,'dJhydmax':30.1,'kGlc':8.7,'nhyd':2.7}
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

def atualiza_par(t):
  global u, user_u, old_u, old_t, inclinacao_u, min

  dy = 2.5/min

  for i in range(0,3):
    if(user_u[i] == 0):
      user_u[i] = 0.0001

    if(u[i] != user_u[i]):
      if(u[i] == old_u[i]):
        inclinacao_u[i] = (user_u[i] - u[i])/dt
        old_t[i] = t
        u[i] = old_u[i] + (inclinacao_u[i])
      else:
        if(abs(1-(t - old_t[i])/dt) < dy):
          u[i] = user_u[i]
          old_u[i] = user_u[i]
          inclinacao_u[i] = 0
        else:
          u[i] = old_u[i] + (inclinacao_u[i]*(t - old_t[i]))
          if (u[i] == user_u[i]):
              old_u[i] = u[i]
              
  #print(str(u[1]) + "  " + str(user_u[1]) + "  " + str(old_u[1]) + "  " + str(1-(t - old_t[1])/dt))

def animate(j):
    
  global i, n, z0, y0, u, t, figure, root, figureplot;
  global NADHm, ADPm, Cam, dp, ADPi, dJhyd, Jo, JHres, JHleak, JpF1, JHF1, Juni, JNaCa, Jant, Jred, JpTCA, Jpgly, Jhyd, dJhydss

  batch_size = 10
  script = False

  while ((i in range(1,N_ADP)) & enable == 1):
    # span for next time step
    tspan = [t[i-1],t[i]]
    # solve for next step
    z = odeint(modelo1,z0,tspan,args=(outros,tabela1,tabela2,tabela3,tabela4,tabela5,betas,table3,inputs,u,))
    # store solution for plotting
    NADHm[i] = z[1][0]
    ADPm[i] = z[1][1]
    Cam[i] = z[1][2]
    dp[i] = z[1][3]
    ADPi[i] = z[1][4]
    dJhyd[i] = z[1][5]
    Jo[i] = find_Jo(outros,tabela1,tabela2,NADHm[i],dp[i])
    JHres[i] = find_JHres(outros,tabela1,tabela2,NADHm[i],dp[i])
    JHleak[i] = find_JHleak(outros,tabela1,dp[i],u[inputs.get("fccp")])
    JpF1[i] = find_JpF1(outros,tabela1,tabela3,ADPm[i],dp[i],u[inputs.get("oligo")])
    JHF1[i] = find_JHF1(outros,tabela1,tabela3,ADPm[i],dp[i],u[inputs.get("oligo")])
    Juni[i] = find_Juni(outros,tabela1,tabela5,Ca_c,dp[i])
    JNaCa[i] = find_JNaCa(outros,tabela1,tabela5,Cam[i],dp[i])
    Jant[i] = find_Jant(outros,tabela1,tabela4,ADPm[i],ADPi[i],dp[i])
    Jred[i] = find_Jred(outros,betas,table3,NADHm[i],((outros['Atoti'])-ADPi[i]),Cam[i],u[inputs.get("Glc")])
    JpTCA[i] = find_JpTCA(outros,betas,table3,NADHm[i],((outros['Atoti'])-ADPi[i]),Cam[i],u[inputs.get("Glc")])
    Jpgly[i] = find_Jpgly(betas,((outros['Atoti'])-ADPi[i]),u[inputs.get("Glc")])
    Jhyd[i] = find_Jhyd(table3,((outros['Atoti'])-ADPi[i]),dJhyd[i])
    dJhydss[i] = find_dJhydss(table3,u[inputs.get("Glc")])
    # next initial condition
    z0 = z[1]
    if script == True:
        u = tscript(t[i], u)
    atualiza_par(t[i])
    i += 1;
        
    if i%batch_size == 0:
      break;

  if i == N_ADP:  
    y0 = [NADHm[i-1],ADPm[i-1],Cam[i-1],dp[i-1]];

  while ((i in range(N_ADP,N)) & enable == 1):
  # span for next time step
    tspan = [t[i-1],t[i]]
    # solve for next step
    y = odeint(modelo2,y0,tspan,args=(outros,tabela1,tabela2,tabela3,tabela4,tabela5,betas,table3,inputs,u,))
    # store solution for plotting
    NADHm[i] = y[1][0]
    ADPm[i] = y[1][1]
    Cam[i] = y[1][2]
    dp[i] = y[1][3]
    ADPi[i] = ADP_c
    dJhyd[i] = d_Jhyd
    Jo[i] = find_Jo(outros,tabela1,tabela2,NADHm[i],dp[i])
    JHres[i] = find_JHres(outros,tabela1,tabela2,NADHm[i],dp[i])
    JHleak[i] = find_JHleak(outros,tabela1,dp[i],u[inputs.get("fccp")])
    JpF1[i] = find_JpF1(outros,tabela1,tabela3,ADPm[i],dp[i],u[inputs.get("oligo")])
    JHF1[i] = find_JHF1(outros,tabela1,tabela3,ADPm[i],dp[i],u[inputs.get("oligo")])
    Juni[i] = find_Juni(outros,tabela1,tabela5,Ca_c,dp[i])
    JNaCa[i] = find_JNaCa(outros,tabela1,tabela5,Cam[i],dp[i])
    Jant[i] = find_Jant2(outros,tabela1,tabela4,ADPm[i],ADPi[i],tot,dp[i])
    Jred[i] = find_Jred(outros,betas,table3,NADHm[i],ATP_c,Cam[i],u[inputs.get("Glc")])
    JpTCA[i] = find_JpTCA(outros,betas,table3,NADHm[i],ATP_c,Cam[i],u[inputs.get("Glc")])
    # next initial condition
    y0 = y[1]
    if script == True:
        u = tscript(t[i], u)
    atualiza_par(t[i])
    i += 1;
        
    if i%batch_size == 0:
      break;

  i -= 1;


  max_of_y = 22.6;
  min_of_y = 20.8;
  delta = 0.01

  if i > 10:
    amax = np.amax(Jo[10:i])
    if amax > max_of_y:
      max_of_y = amax + abs(amax)*delta;
      
    amin = np.amin(Jo[10:i])
    if amin < min_of_y:
      min_of_y = amin - abs(amin)*delta;

  figureplot.clear()
  figureplot.plot(t_n[10:i],Jo[10:i],'r',label='Jo', linewidth=2.5)
  figureplot.set_xlabel('time (min)',fontsize=16)
  figureplot.set_ylabel('Jo (nmol/min/mg-protein)',fontsize=16)
  figureplot.set_ylim(min_of_y, max_of_y)
  figureplot.legend(loc='best')
  figureplot.set_title('Simulação de Oximetria', fontdict=None, loc='center')
  
  i += 1;

#==============================================================================================================================================================

# MODELOS ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# Antes da adição de ADP:
def modelo1(ci,t,outros,tabela1,tabela2,tabela3,tabela4,tabela5,betas,table3,inputs,u):
  Atoti = outros['Atoti']
  CCa2m = table3['CCa2m']
  Cmito = table3['Cmito']
  tauhyd = table3['tauhyd']
  Cai = 0.2         
  NADHm = ci[0]
  ADPm = ci[1]
  Cam = ci[2]
  dp = ci[3]
  ADPi = ci[4]
  dJhyd = ci[5]
  ATPi = Atoti-ADPi
  nadhm = (find_Jred(outros,betas,table3,NADHm,ATPi,Cam,u[inputs.get("Glc")]))-(find_Jo(outros,tabela1,tabela2,NADHm,dp))
  adpm = (find_Jant(outros,tabela1,tabela4,ADPm,ADPi,dp))-(find_JpTCA(outros,betas,table3,NADHm,ATPi,Cam,u[inputs.get("Glc")]))-(find_JpF1(outros,tabela1,tabela3,ADPm,dp,u[inputs.get("oligo")]))
  cam = (1/CCa2m)*((find_Juni(outros,tabela1,tabela5,Cai,dp))-(find_JNaCa(outros,tabela1,tabela5,Cam,dp)))
  dpm = (-1/Cmito)*((find_JHF1(outros,tabela1,tabela3,ADPm,dp,u[inputs.get("oligo")]))+(find_JHleak(outros,tabela1,dp,u[inputs.get("fccp")]))+(find_Jant(outros,tabela1,tabela4,ADPm,ADPi,dp))+(2*(find_Juni(outros,tabela1,tabela5,Cai,dp)))+(find_JNaCa(outros,tabela1,tabela5,Cam,dp))-(find_JHres(outros,tabela1,tabela2,NADHm,dp)))
  adpi = (0.09/60000)*((find_Jhyd(table3,ATPi,dJhyd))-(find_Jpgly(betas,ATPi,u[inputs.get("Glc")]))-(find_Jant(outros,tabela1,tabela4,ADPm,ADPi,dp)))
  djhyd = (1/tauhyd)*((find_dJhydss(table3,u[inputs.get("Glc")]))-dJhyd)
  dzdt = [nadhm,adpm,cam,dpm,adpi,djhyd]
  return dzdt
# Depois da adição de ADP:
def modelo2(ci,t,outros,tabela1,tabela2,tabela3,tabela4,tabela5,betas,table3,inputs,u):
    CCa2m = table3['CCa2m']
    Cmito = table3['Cmito']
    Cai = 0.2         
    NADHm = ci[0]
    ADPm = ci[1]
    Cam = ci[2]
    dp = ci[3]
    # Ajuste da razão ATPi/ADPi após a adição de ADP -----------------------------------------------------
    ADPi = 2.0
    ATPi = 8.0
    total = ADPi+ATPi
    # ----------------------------------------------------------------------------------------------------
    nadhm = (find_Jred(outros,betas,table3,NADHm,ATPi,Cam,u[inputs.get("Glc")]))-(find_Jo(outros,tabela1,tabela2,NADHm,dp))
    adpm = (find_Jant2(outros,tabela1,tabela4,ADPm,ADPi,total,dp))-(find_JpTCA(outros,betas,table3,NADHm,ATPi,Cam,u[inputs.get("Glc")]))-(find_JpF1(outros,tabela1,tabela3,ADPm,dp,u[inputs.get("oligo")]))
    cam = (1/CCa2m)*((find_Juni(outros,tabela1,tabela5,Cai,dp))-(find_JNaCa(outros,tabela1,tabela5,Cam,dp)))
    dpm = (-1/Cmito)*((find_JHF1(outros,tabela1,tabela3,ADPm,dp,u[inputs.get("oligo")]))+(find_JHleak(outros,tabela1,dp,u[inputs.get("fccp")]))+(find_Jant2(outros,tabela1,tabela4,ADPm,ADPi,total,dp))+(2*(find_Juni(outros,tabela1,tabela5,Cai,dp)))+(find_JNaCa(outros,tabela1,tabela5,Cam,dp))-(find_JHres(outros,tabela1,tabela2,NADHm,dp)))
    dzdt = [nadhm,adpm,cam,dpm]
    return dzdt
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FLUXOS ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def find_Jo(outros,tabela1,tabela2,NADHm,dp):
    F = outros['F']
    R = outros['R']
    NADtot = outros['NADtot']
    deltapH = tabela1['deltapH']
    T = tabela1['T']
    Kres = tabela2['Kres']
    pres = tabela2['pres']
    r1 = tabela2['r1']
    r2 = tabela2['r2']
    r3 = tabela2['r3']
    ra = tabela2['ra']
    rc1 = tabela2['rc1']
    rc2 = tabela2['rc2']
    dpB = tabela2['dpB']
    g = tabela2['g']
    FRT = F/(R*T)
    NADm = NADtot-NADHm
    Ares = ((R*T)/F)*(np.log((Kres*(np.sqrt(NADHm)) + 0.0001)/(np.sqrt(NADm) + 0.0001)))
    n = (((ra*(np.power(10,(6*deltapH))))+(rc1*(np.exp(6*FRT*dpB))))*(np.exp(FRT*Ares)))-(ra*(np.exp(g*6*FRT*dp)))+(rc2*(np.exp(FRT*Ares))*(np.exp(g*6*FRT*dp)))
    d = ((1+(r1*(np.exp(FRT*Ares))))*(np.exp(6*FRT*dpB)))+((r2+(r3*(np.exp(FRT*Ares))))*(np.exp(g*6*FRT*dp)))
    Jo = 30*pres*(n/d)
    return Jo

def find_JHres(outros,tabela1,tabela2,NADHm,dp):
    F = outros['F']
    R = outros['R']
    NADtot = outros['NADtot']
    deltapH = tabela1['deltapH']
    T = tabela1['T']
    Kres = tabela2['Kres']
    pres = tabela2['pres']
    r1 = tabela2['r1']
    r2 = tabela2['r2']
    r3 = tabela2['r3']
    ra = tabela2['ra']
    rb = tabela2['rb']
    dpB = tabela2['dpB']
    g = tabela2['g']
    FRT = F/(R*T)
    NADm = NADtot-NADHm
    Ares = ((R*T)/F)*(np.log((Kres*(np.sqrt(NADHm)))/(np.sqrt(NADm))))
    n = (ra*(np.power(10,(6*deltapH)))*(np.exp(FRT*Ares)))-((ra+rb)*(np.exp(g*6*FRT*dp)))
    d = ((1+(r1*(np.exp(FRT*Ares))))*(np.exp(6*FRT*dpB)))+((r2+(r3*(np.exp(FRT*Ares))))*(np.exp(g*6*FRT*dp)))
    JHres = 360*pres*(n/d)
    return JHres

def find_JHleak(outros,tabela1,dp,fccp):
    F = outros['F']
    R = outros['R']
    deltapH = tabela1['deltapH']
    T = tabela1['T']
    gH = tabela1['gH']
    Z = 2.303*((R*T)/F)
    deltap = dp-(Z*deltapH)
    JHleak = gH*deltap*fccp
    return JHleak
    
def find_JpF1(outros,tabela1,tabela3,ADPm,dp,oligo):
    F = outros['F']
    R = outros['R']
    Atot = outros['Atot']
    deltapH = tabela1['deltapH']
    T = tabela1['T']
    KF1 = tabela3['KF1']
    Pim = tabela3['Pim']
    pF1 = tabela3['pF1']
    p1 = tabela3['p1']
    p2 = tabela3['p2']
    p3 = tabela3['p3']
    pa = tabela3['pa']
    pc1 = tabela3['pc1']
    pc2 = tabela3['pc2']
    dpB = tabela3['dpB']
    FRT = F/(R*T)
    ATPm = Atot-ADPm
    ADPmfree = 0.8*ADPm 
    AF1 = ((R*T)/F)*(np.log((KF1*ATPm)/(ADPmfree*Pim)))
    n = (((pa*(np.power(10,(3*deltapH))))+(pc1*(np.exp(3*FRT*dpB))))*(np.exp(FRT*AF1)))-(pa*(np.exp(3*FRT*dp)))+(pc2*(np.exp(FRT*AF1))*(np.exp(3*FRT*dp)))
    d = ((1+(p1*(np.exp(FRT*AF1))))*(np.exp(3*FRT*dpB)))+((p2+(p3*(np.exp(FRT*AF1))))*(np.exp(3*FRT*dp)))
    JpF1 = -60*pF1*(n/d)*oligo
    return JpF1
    
def find_JHF1(outros,tabela1,tabela3,ADPm,dp,oligo):
    F = outros['F']
    R = outros['R']
    Atot = outros['Atot']
    deltapH = tabela1['deltapH']
    T = tabela1['T']
    KF1 = tabela3['KF1']
    Pim = tabela3['Pim']
    pF1 = tabela3['pF1']
    p1 = tabela3['p1']
    p2 = tabela3['p2']
    p3 = tabela3['p3']
    pa = tabela3['pa']
    pb = tabela3['pb']
    dpB = tabela3['dpB']
    FRT = F/(R*T)
    ATPm = Atot-ADPm
    ADPmfree = 0.8*ADPm 
    AF1 = ((R*T)/F)*(np.log((KF1*ATPm + 0.0001)/(ADPmfree*Pim)))
    n = (pa*(np.power(10,(3*deltapH)))*(np.exp(FRT*AF1)))-((pa+pb)*(np.exp(3*FRT*dp)))
    d = ((1+(p1*(np.exp(FRT*AF1))))*(np.exp(3*FRT*dpB)))+((p2+(p3*(np.exp(FRT*AF1))))*(np.exp(3*FRT*dp)))
    JHF1 = -180*pF1*(n/d)*oligo
    return JHF1

def find_Jant(outros,tabela1,tabela4,ADPm,ADPi,dp):
    F = outros['F']
    R = outros['R']
    Atot = outros['Atot']
    Atoti = outros['Atoti']
    T = tabela1['T']
    JmaxANT = tabela4['JmaxANT']
    f = tabela4['f']
    FRT = F/(R*T)
    ATPm = Atot-ADPm
    ATPi = Atoti-ADPi
    ATP4m = 0.05*ATPm
    ATP4i = 0.05*ATPi
    ADP3m = 0.36*ADPm
    ADP3i = 0.135*ADPi
    n = 1-((ATP4i/ADP3i)*(ADP3m/ATP4m)*(np.exp(-1*FRT*dp)))
    d = (1+((ATP4i/ADP3i)*(np.exp(-1*f*FRT*dp))))*(1+(ADP3m/ATP4m))
    Jant = JmaxANT*(n/d)
    return Jant   

def find_Jant2(outros,tabela1,tabela4,ADPm,ADPi,total,dp):              
    F = outros['F']
    R = outros['R']
    Atot = outros['Atot']
    T = tabela1['T']
    JmaxANT = tabela4['JmaxANT']
    f = tabela4['f']
    FRT = F/(R*T)
    ATPm = Atot-ADPm
    ATPi = total-ADPi
    ATP4m = 0.05*ATPm
    ATP4i = 0.05*ATPi
    ADP3m = 0.36*ADPm
    ADP3i = 0.135*ADPi
    n = 1-((ATP4i/ADP3i)*(ADP3m/ATP4m)*(np.exp(-1*FRT*dp)))
    d = (1+((ATP4i/ADP3i)*(np.exp(-1*f*FRT*dp))))*(1+(ADP3m/ATP4m))
    Jant = JmaxANT*(n/d)
    return Jant
    
def find_Juni(outros,tabela1,tabela5,Cai,dp):
    F = outros['F']
    R = outros['R']
    T = tabela1['T']
    Kcat = tabela5['Kcat']
    Kact = tabela5['Kact']
    Jmaxuni = tabela5['Jmaxuni']
    dpuni = tabela5['dpuni']
    nact = tabela5['nact']
    Lunimax = tabela5['Lunimax']
    FRT = F/(R*T)
    t1 = ((Cai/Kcat)*(np.power((1+(Cai/Kcat)),3)))/((np.power((1+(Cai/Kcat)),4))+(Lunimax/(np.power((1+(Cai/Kact)),nact))))
    t2 = (2*FRT*(dp-dpuni))/(1-(np.exp(-2*FRT*(dp-dpuni))))
    Juni = t1*Jmaxuni*t2
    return Juni
    
def find_JNaCa(outros,tabela1,tabela5,Cam,dp):
    F = outros['F']
    R = outros['R']
    T = tabela1['T']
    KNa = tabela5['KNa']
    KCa = tabela5['KCa']
    Nai = tabela5['Nai']
    dpNaCa = tabela5['dpNaCa']
    b = tabela5['b']
    JmaxNaCa = tabela5['JmaxNaCa']
    n = tabela5['n']
    FRT = F/(R*T)
    JNaCa = JmaxNaCa*(np.exp(b*FRT*(dp-dpNaCa)))/((np.power((1+(KNa/Nai)),n))*(1+(KCa/Cam)))
    return JNaCa

def find_dJglytotal(betas,ATPi,Glc):
    betamax = betas['betamax']
    beta1 = betas['beta1']
    beta2 = betas['beta2']
    beta3 = betas['beta3']
    beta4 = betas['beta4']
    beta5 = betas['beta5']
    beta6 = betas['beta6']
    beta7 = betas['beta7']
    dJglytotal = (betamax*(1+(beta1*Glc))*beta2*Glc*ATPi)/(1+(beta3*ATPi)+((1+(beta4*ATPi))*beta5*Glc)+((1+(beta6*ATPi))*beta7*Glc*Glc))
    return dJglytotal

def find_fPDH(outros,table3,NADHm,Cam):
  NADtot = outros['NADtot']
  u1 = table3['u1']
  u2 = table3['u2']
  KCa2 = table3['KCa2']
  NADm = NADtot-NADHm
  fPDH = 1/(1+(u2*((1+(u1*(1/(np.power((1+(Cam/KCa2)),2)))))*((NADHm/NADm)+1))))
  return fPDH

def find_Jred(outros,betas,table3,NADHm,ATPi,Cam,Glc):
    Jredbasal = table3['Jredbasal']
    dJglytotal = find_dJglytotal(betas,ATPi,Glc)
    fPDH = find_fPDH(outros,table3,NADHm,Cam)
    Jred = Jredbasal+(7.36*fPDH*dJglytotal)
    return Jred

def find_JpTCA(outros,betas,table3,NADHm,ATPi,Cam,Glc):
    Jredbasal = table3['Jredbasal']
    dJglytotal = find_dJglytotal(betas,ATPi,Glc)
    fPDH = find_fPDH(outros,table3,NADHm,Cam)
    JpTCA = (Jredbasal/3)+(0.84*fPDH*dJglytotal)
    return JpTCA

def find_Jpgly(betas,ATPi,Glc):
    dJglytotal = find_dJglytotal(betas,ATPi,Glc)
    Jpgly = 2*dJglytotal
    return Jpgly 

def find_Jhyd(table3,ATPi,dJhyd):
    khyd = table3['khyd']
    Jhyd = (khyd*ATPi)+dJhyd
    return Jhyd

def find_dJhydss(table3,Glc):
    dJhydmax = table3['dJhydmax']
    kGlc = table3['kGlc']
    nhyd = table3['nhyd']
    dJhydss = dJhydmax/(1+(np.power((kGlc/Glc),nhyd)))
    return dJhydss
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Aplicação das entradas
def tscript(t, u):

    global N_ADP, min    

    if(t < t_Glc):
        user_u[inputs.get("Glc")] = Glc_a
        user_u[inputs.get("oligo")] = oligo_a
        user_u[inputs.get("fccp")] = fccp_a
        N_ADP = 4*min

    if(t >= t_Glc and t < t_oligo):
        user_u[inputs.get("Glc")] = Glc_d

    if(t >= t_oligo and t < t_fccp):
        user_u[inputs.get("oligo")] = oligo_d

    if(t >= t_fccp):
        user_u[inputs.get("fccp")] = fccp_d

    return u

#==============================================================================================================================================================

# TEMPO  
um_min = 60000                                      # Milisegundos em um minuto
min = 100  # não reduzir de 100                     # Quantidade de pontos por minuto de simulação
q_min = 10                                          # Quantidade de minutos da simulação
t_i = 0.0                                           # Tempo inicial
t_f = um_min*q_min                                  # Tempo final
N = min*q_min                                       # Número total de pontos                                       
t = np.linspace(t_i,t_f,N)                          # Vetor tempo (ms)
t_n = t*(1/um_min)                                  # Vetor tempo (min)
dt = 10000                                          # Intervalo de tempo da rampa

# ENTRADAS --------------------------------------------------------------------------------------------------------------------------------------------------------------
inputs = {"Glc":0,"oligo":1,"fccp":2}

# Parâmetros Glc
user_GLc = 2.8

Glc_a = 2.8
Glc_d = 7.0
t_Glc = 120000                                  # (ms)
inclinacao_Glc = (Glc_d-Glc_a)/dt
constante_Glc = Glc_a-(inclinacao_Glc*t_Glc) 
# Parâmetros ADP
t_ADP = 4                                       # (min)

# Parâmetros oligo
user_oligo = 1.0

oligo_a = 1.0
oligo_adi = 0.2																	# 13/08			range_oligo = 0 até 0.2 uL 
oligo_d = 1-(4.7*oligo_adi)											# 13/08 (0.06)
t_oligo = 360000                                # (ms)
inclinacao_oligo = (oligo_d-oligo_a)/dt
constante_oligo = oligo_a-(inclinacao_oligo*t_oligo)

# Parâmetros fccp
user_fccp = 1.0

fccp_a = 1.0
fccp_adi = 1.0																    # 13/08			range_fccp = 0 até 1.0 uL
fccp_d = 1+(6.5*fccp_adi)														# 13/08 (7.5)
t_fccp = 480000                                 # (ms)
inclinacao_fccp = (fccp_d-fccp_a)/dt
constante_fccp = fccp_a-(inclinacao_fccp*t_fccp)

# Aplicação de ADP
N_ADP = N

# Inicialmente
u = {};
u[inputs.get("Glc")] = Glc_a
u[inputs.get("oligo")] = oligo_a
u[inputs.get("fccp")] = fccp_a
user_u = {};
user_u[inputs.get("Glc")] = Glc_a
user_u[inputs.get("oligo")] = oligo_a
user_u[inputs.get("fccp")] = fccp_a
old_u = {};
old_u[inputs.get("Glc")] = Glc_a
old_u[inputs.get("oligo")] = oligo_a
old_u[inputs.get("fccp")] = fccp_a
old_t = {};
old_t[inputs.get("Glc")] = 0
old_t[inputs.get("oligo")] = 0
old_t[inputs.get("fccp")] = 0
inclinacao_u = {};
inclinacao_u[inputs.get("Glc")] = 0
inclinacao_u[inputs.get("oligo")] = 0
inclinacao_u[inputs.get("fccp")] = 0

# Condições Iniciais 
NADHm = 0.02
ADPm = 8.2
Cam = 0.004
dp = 150   
ADPi = 0.6
dJhyd = 0.7                                            
z0 = [NADHm,ADPm,Cam,dp,ADPi,dJhyd]
y0 = [0, 0, 0, 0]

enable = 0;
i = 2

# Armazenamento de soluções
NADHm = np.empty_like(t)
ADPm = np.empty_like(t)
Cam = np.empty_like(t)
dp = np.empty_like(t)
ADPi = np.empty_like(t)
dJhyd = np.empty_like(t)
Jo = np.empty_like(t)
JHres = np.empty_like(t)
JHleak = np.empty_like(t)
JpF1 = np.empty_like(t)
JHF1 = np.empty_like(t)
Juni = np.empty_like(t)
JNaCa = np.empty_like(t)
Jant = np.empty_like(t)
Jred = np.empty_like(t)
JpTCA = np.empty_like(t)
Jpgly = np.empty_like(t)
Jhyd = np.empty_like(t)
dJhydss = np.empty_like(t)

# Armazenamento de condições iniciais
NADHm[0] = z0[0]
ADPm[0] = z0[1]
Cam[0] = z0[2]
dp[0] = z0[3]
ADPi[0] = z0[4]
dJhyd[0] = z0[5]
Jo[0] = find_Jo(outros,tabela1,tabela2,NADHm[0],dp[0])
JHres[0] = find_JHres(outros,tabela1,tabela2,NADHm[0],dp[0])
JHleak[0] = find_JHleak(outros,tabela1,dp[0],fccp_a)
JpF1[0] = find_JpF1(outros,tabela1,tabela3,ADPm[0],dp[0],oligo_a)
JHF1[0] = find_JHF1(outros,tabela1,tabela3,ADPm[0],dp[0],oligo_a)
Juni[0] = find_Juni(outros,tabela1,tabela5,0.2,dp[0])
JNaCa[0] = find_JNaCa(outros,tabela1,tabela5,Cam[0],dp[0])
Jant[0] = find_Jant(outros,tabela1,tabela4,ADPm[0],ADPi[0],dp[0])
Jred[0] = find_Jred(outros,betas,table3,NADHm[0],((outros['Atoti'])-ADPi[0]),Cam[0],Glc_a)
JpTCA[0] = find_JpTCA(outros,betas,table3,NADHm[0],((outros['Atoti'])-ADPi[0]),Cam[0],Glc_a)
Jpgly[0] = find_Jpgly(betas,((outros['Atoti'])-ADPi[0]),Glc_a)
Jhyd[0] = find_Jhyd(table3,((outros['Atoti'])-ADPi[0]),dJhyd[0])
dJhydss[0] = find_dJhydss(table3,Glc_a)

# Para os gráficos 
# ---------- Durante toda a corrida
Ca_c = 0.2                                              
# ---------- Após a adição de ADP (segundo laço for)
ADP_c = 2.0
ATP_c = 8.0
tot = ADP_c+ATP_c
d_Jhyd = 0.0

# INTERFACE GRÁFICA
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Page(tk.Frame):
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
    def show(self):
        self.lift()

class Page1(Page):
    def __init__(self, *args, **kwargs):
      Page.__init__(self, *args, **kwargs)
      
      canvas = FigureCanvasTkAgg(figure, self)
      canvas.draw()
      #canvas.grid(row = 10, column = 10)
      canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
      
      #toolbar = NavigationToolbar2Tk(canvas, self)
      #toolbar.update()
      #canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

      variablescontainer = tk.Frame(self)
      variablescontainer.pack(side="top", anchor = "ne")
        
      espaco1 = tk.Label(variablescontainer, text="", width = 1, height = 1)
      espaco2 = tk.Label(variablescontainer, text="", width = 1, height = 1)
      espaco3 = tk.Label(variablescontainer, text="", width = 1, height = 1)
      espaco4 = tk.Label(variablescontainer, text="", width = 1, height = 1)
      
      espaco1.grid(row = 1, column = 1)
      espaco2.grid(row = 1, column = 1000)
      espaco3.grid(row = 1000, column = 1)
      espaco4.grid(row = 1000, column = 1000)
      
      espaco_label = tk.Label(variablescontainer, text="", width = 1, height = 1)
      espaco_label.grid(row = 0, column = 1)
      
      espaco_entry = tk.Label(variablescontainer, text="", width = 1, height = 1)
      espaco_entry.grid(row = 0, column = 15)
      
      espaco_button1 = tk.Label(variablescontainer, text="", width = 1, height = 1)
      espaco_button1.grid(row = 0, column = 25)
      
      #espaco_button2 = tk.Label(variablescontainer, text="", width = 1, height = 1)
      #espaco_button2.grid(row = 0, column = 35)
      
      espaco_button3 = tk.Label(variablescontainer, text="", width = 1, height = 1)
      espaco_button3.grid(row = 0, column = 45)
      
      self.label1 = tk.Label(variablescontainer, text = str(u[inputs.get("Glc")]) + " mM")
      self.label1.grid(row = 10, column = 0, sticky = 'E')
      
      self.entry1 = tk.Entry(variablescontainer)
      self.entry1.grid(row=10,column=10)
      #self.entry1.grid(row=1,column=1)

      btn1 = tk.Button(variablescontainer, text="Glc", command=self.muda_par1, width = 15)
      btn1.grid(row=10,column=20)

      self.label2 = tk.Label(variablescontainer, text = "0.0 uL")
      self.label2.grid(row = 20, column = 0, sticky = 'E')

      self.entry2 = tk.Entry(variablescontainer)
      self.entry2.grid(row=20,column=10)

      btn2 = tk.Button(variablescontainer, text="oligo", command=self.muda_par2, width = 15)
      btn2.grid(row=20,column=20)

      self.label3 = tk.Label(variablescontainer, text = "0.0 uL")
      self.label3.grid(row = 30, column = 0, sticky = 'E')
      
      self.entry3 = tk.Entry(variablescontainer)
      self.entry3.grid(row=30,column=10)

      btn3 = tk.Button(variablescontainer, text="fccp", command=self.muda_par3, width = 15)
      btn3.grid(row=30,column=20)

      self.btnADP = tk.Button(variablescontainer, text="Adicionar ADP", command=self.add_ADP, width = 15, bg = 'yellow')
      self.btnADP.grid(row=40,column=20)

      self.btnenable = tk.Button(variablescontainer, text="Iniciar", command=self.muda_enable, width = 10, height = 3, bg = 'green', fg = 'black')
      self.btnenable.grid(row=10,column=30, rowspan = 20)
      self.btnenable.bind("<Return>", self.muda_enable)
      
      #self.apply = tk.Button(variablescontainer, text="apply all", command=self.apply2all, width = 10, height = 3)
      #self.apply.grid(row=10,column=40, rowspan = 20)
      
      self.reset = tk.Button(variablescontainer, text="reset", command=self.resetride, width = 10, height = 3)
      self.reset.grid(row=10,column=50, rowspan = 20)
      

    def muda_par1(self):
      global user_u 

      entry = float(self.entry1.get());

      if (entry >= 7.0):
        entry = 7.0;

      if (entry > user_u[inputs.get("Glc")]):
        user_u[inputs.get("Glc")] = entry
        self.label1["text"] = str(round(user_u[inputs.get("Glc")], 2)) + " mM"

    def muda_par2(self):
      global user_u

      entry = float(self.entry2.get());

      if (entry >= 0.2):
        entry = 0.2;

      if (1-(4.7*entry) < user_u[inputs.get("oligo")]):
        user_u[inputs.get("oligo")] = 1-(4.7*entry)
        self.label2["text"] = str(entry) + " uL"

    def muda_par3(self):
      global user_u 

      entry = float(self.entry3.get());

      if (entry >= 1.0):
        entry = 1.0;

      if (1+(6.5*entry) > user_u[inputs.get("fccp")]):
        user_u[inputs.get("fccp")] = 1+(6.5*entry)
        self.label3["text"] = str(entry) + " uL"

    def muda_enable(self, event=None):
      global enable

      if enable == 1:
        enable = 0;
        self.btnenable["text"] = "Retomar";
        self.btnenable["bg"] = 'green'
        self.btnenable["fg"] = 'black'
      else:
        enable = 1;
        self.btnenable["text"] = "Pausar";
        self.btnenable["bg"] = 'red'
        self.btnenable["fg"] = 'white'
        
    def apply2all(self):
      global u 

      self.muda_par1()
      self.muda_par2()
      self.muda_par3()
      
    def add_ADP(self):
        global N_ADP, i
        
        if N_ADP > i:
            N_ADP = i;   
            self.btnADP["text"] = "ADP Adicionado!";
            self.btnADP["bg"] = 'gray'
            
    def remove_ADP(self):
        global N_ADP, N
        
        N_ADP = N;
        self.btnADP["text"] = "Adicionar ADP";
        self.btnADP["bg"] = 'yellow'
            
    def resetride(self):
        global i, z0, y0, u, N_ADP, old_t, inclinacao_u, N, enable
        
        enable = 0;

        i = 2
        NADHm = 0.02
        ADPm = 8.2
        Cam = 0.004
        dp = 150   
        ADPi = 0.6
        dJhyd = 0.7
        z0 = [NADHm,ADPm,Cam,dp,ADPi,dJhyd]
        y0 = [0, 0, 0, 0]
        self.remove_ADP();
        
        self.label1["text"] = str(Glc_a) + " mM"
        self.label2["text"] = "0.0 uL"
        self.label3["text"] = "0.0 uL"

        self.btnenable["text"] = "Iniciar";
        self.btnenable["bg"] = 'green'
        self.btnenable["fg"] = 'black'
        
        u[inputs.get("Glc")] = user_u[inputs.get("Glc")] = old_u[inputs.get("Glc")] = Glc_a
        u[inputs.get("oligo")] = user_u[inputs.get("oligo")] = old_u[inputs.get("oligo")] = oligo_a
        u[inputs.get("fccp")] = user_u[inputs.get("fccp")] = old_u[inputs.get("fccp")] = fccp_a
        inclinacao_u = dict.fromkeys(inclinacao_u, 0)
        old_t = dict.fromkeys(old_t, 0)

class Page2(Page):
   def __init__(self, *args, **kwargs):
      Page.__init__(self, *args, **kwargs)
      #label = tk.Label(self, text="This is page 2")
      #label.pack(side="top", fill="both", expand=True)

      textcontainer = tk.Frame(self)
      textcontainer.pack(side="top", fill="both", expand=True)

      content = "Falha ao abrir arquivo \"info.txt\""
      text_area = st.ScrolledText(textcontainer, font = (17), wrap = "word", borderwidth = 3, highlightbackground="grey", relief=tk.SUNKEN) 
        
      text_area.pack(side="top", fill="both", expand=True, padx = (7, 3), pady = (3, 7)) 

      if os.path.exists('info.txt') == True:

        file = open("info.txt","r", encoding="utf8")
        content = file.read() 
        file.close()

      # Inserting Text which is read only 
      text_area.insert(tk.INSERT, content) 


class Page3(Page):
   def __init__(self, *args, **kwargs):
      Page.__init__(self, *args, **kwargs)

      textcontainer = tk.Frame(self)
      textcontainer.pack(side="top", fill="both", expand=True)

      content = "Falha ao abrir arquivo \"tutorial.txt\""
      text_area = st.ScrolledText(textcontainer, font = (17), wrap = "word", borderwidth = 3, highlightbackground="grey", relief=tk.SUNKEN) 
        
      text_area.pack(side="top", fill="both", expand=True, padx = (7, 3), pady = (3, 7)) 

      if os.path.exists('tutorial.txt') == True:

        file = open("tutorial.txt","r", encoding="utf8")
        content = file.read() 
        file.close()

      # Inserting Text which is read only 
      text_area.insert(tk.INSERT, content)

class Page4(Page):
   def __init__(self, *args, **kwargs):
      Page.__init__(self, *args, **kwargs)

      textcontainer = tk.Frame(self)
      textcontainer.pack(side="top", fill="both", expand=True)

      content = "Falha ao abrir arquivo \"about.txt\""
      text_area = st.ScrolledText(textcontainer, font = (17), wrap = "word", borderwidth = 3, highlightbackground="grey", relief=tk.SUNKEN) 
        
      text_area.pack(side="top", fill="both", expand=True, padx = (7, 3), pady = (3, 7)) 

      if os.path.exists('about.txt') == True:

        file = open("about.txt","r", encoding="utf8")
        content = file.read() 
        file.close()

      # Inserting Text which is read only 
      text_area.insert(tk.INSERT, content)

class MainView(tk.Frame):
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        self.master.title("Simulação de Oximetria")
        p1 = Page1(self)
        p2 = Page2(self)
        p3 = Page3(self)
        p4 = Page4(self)

        buttonframe = tk.Frame(self)
        container = tk.Frame(self)
        buttonframe.pack(side="top", fill="x", expand=False)
        container.pack(side="top", fill="both", expand=True)

        p1.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p2.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p3.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p4.place(in_=container, x=0, y=0, relwidth=1, relheight=1)

        b1 = tk.Button(buttonframe, text="Simulador", command=p1.lift)
        b2 = tk.Button(buttonframe, text="Informações", command=p2.lift)
        b3 = tk.Button(buttonframe, text="Ajuda", command=p3.lift)
        b4 = tk.Button(buttonframe, text="Sobre", command=p4.lift)

        b1.pack(side="left")
        b2.pack(side="left")
        b3.pack(side="left")
        b4.pack(side="left")

        p1.show()



if __name__ == "__main__":

  figure = Figure(figsize=(5, 4), dpi=100)
  figureplot = figure.add_subplot(1, 1, 1)

  LARGE_FONT= ("Verdana", 12)
  style.use("ggplot") #outros estilos em plt.style.available

  root = tk.Tk()
  main = MainView(root)
  main.pack(side="top", fill="both", expand=True)
  root.wm_geometry("900x650+200+30")

  ani = animation.FuncAnimation(figure, animate, interval=200)

  if os.path.exists('cell.ico') == True:
    root.iconbitmap(r'cell.ico') #icone de <div>Icons made by <a href="https://www.flaticon.com/authors/surang" title="surang">surang</a> from <a href="https://www.flaticon.com/" title="Flaticon">www.flaticon.com</a></div>
  root.mainloop()
