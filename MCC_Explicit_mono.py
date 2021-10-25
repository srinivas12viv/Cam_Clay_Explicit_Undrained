import numpy as np
import matplotlib.pyplot as plt
import math

print ("This program is for predicting soil's behavior using MCC model in a Undrained TRIAXIAL TEST")

global p_p_o,M,Lamda,k,T,N

p_p_o=float(input("Enter the maximum Isotropic consolidation stress (Po')(kPa):"))

M,Lamda,k,T,N=np.loadtxt('Mat_para.txt',skiprows=1,unpack=True)  #Material Parameters

v=0.25    #Assumed Poisson Ratio

print ("Please define the initial stress state in p-q space-")

p_p_i=float(input("Enter the initial mean effective confining stress on the sample (kPa):"))
q_i=0

print ("Running Triaxial Undrained Test")

SV= N-(Lamda*np.log(p_p_o))+(k*np.log(p_p_o/p_p_i))
eo=SV-1

Strain_level=float(input("Please Enter the strain level you wish to plot the results(%):"))
Strain_incr=float(input("Enter the Strain Increament(%):"))
de=Strain_incr*0.01
iter_no=Strain_level/Strain_incr

#Memory Allocation
De=np.zeros([6,6])
dfds=np.zeros([6,1])
dfdep=np.zeros([6,1])
u=0                         # Pore Water Pressure     
dStrain=np.zeros([6,1])     # Increamental Strain
epsV=0                      # Volumetic Strain
epsD=0                      # Deviatoric Strain


#Tolerance Value set
FTOL=10e-6  #Yield Surface Tolerance Value
STOL=10e-6
LTOL=10e-6
EPS=10e-16

pres_iter=0
S=np.array([[p_p_i],[p_p_i],[p_p_i],[0],[0],[0]])  #Stress (S) =[s11,s22,s33,t12,t13,t23]
Strain=np.array([[0],[0],[0],[0],[0],[0]])
p=(S[0]+2*S[2])/3                                  #Mean Stress
q=S[0]-S[2]	                                   #Deviatoric Stress
MYP=p_p_o					   #Maximum Yield isotropic Pressure
fai=(p*p)+((q*q)/(M*M))-(p_p_o*p)                  #Yield Surface
VR=eo		                                   #Void Ratio
	                                   


f=open('Res.txt','w')
f.write("p' q ep eq u VR \n")
print("%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f" %(p,q,epsV,epsD,u,VR),file=f)
f.close

def pla_sub(S,De,M,eo,Lamda,k,dT,dStrain_Pla):
  dfds=np.zeros([6,1])
  p=(S[0]+2*S[2])/3                                  #Mean Stress
  q=S[0]-S[2]
  MYP=p+(q*q)/(M*M*p)
  for m in range(0,len(dfds)):
    if m<3:
      dfds[m]= ((M*M)*(((2*p)-MYP)/3))+(3*(S[m]-p))
    if m>=3:
      dfds[int(m)]=0 
    
  Kp=((M**4)*(1+eo)*(MYP)*(p)*(2*p-MYP))/(Lamda-k)
  
  dfdsT=np.matrix.transpose(dfds)
  Num1=np.matmul(De,dfds)
  Num2=np.matmul(dfdsT,De)
  Numerator=np.matmul(Num1,Num2)
  Den1 = np.matmul(dfdsT,De)
  Denominator = np.matmul(Den1,dfds) + Kp
  Plastic = Numerator/Denominator
  D = De-Plastic
  den=dT*dStrain_Pla
  dS=np.matmul(D,den) #Plastic Stress Increament
  S_temp=S+dS
  p_temp=(S_temp[0]+2*S_temp[2])/3                
  q_temp=S_temp[0]-S_temp[2]
  MYP_temp=p_temp+((q_temp*q_temp)/(M*M*p_temp))
  dk=MYP_temp-MYP     #Internal Variable Increament   
  
  return dS,dk


while pres_iter<iter_no:

  BM=SV*p/k                 #Bulk Modulus
  SM=(3*BM*(1-2*v))/(2*(1+v))  #Shear Modulus
   
  #Define De
  for m in range(0, len(De)):
    for n in range(0,len(De)):
      if m<3:
        if m==n:
          De[m,n]= BM+((4/3)*SM)
        elif n<3:
          De[m,n]=BM-((2/3)*SM)
      else:
        if m==n:
          De[m,n]=SM
  

  dStrain=np.array([[de],[-de/2],[-de/2],[0],[0],[0]])
  dStress=np.matmul(De,dStrain)
  S_temp=S+dStress                                #Temporary Stress to check elastic or plastic behavior
  p_temp=(S_temp[0]+2*S_temp[2])/3                #Temporary Mean Stress
  q_temp=S_temp[0]-S_temp[2]	                  #Temporary Deviatoric Stress
  
  fai_temp=(p_temp*p_temp)+((q_temp*q_temp)/(M*M))-(MYP*p_temp) 
  
  if fai_temp<=FTOL:
    S=S_temp
    p=p_temp
    q=q_temp
    fai = fai_temp   

  if fai<=FTOL and fai_temp>FTOL:
    
    #Calculation of alpha using Modified Regular Falsi Intersection Scheme using Negative Plastic Multilier

    ao=0.  #alpha_o
    a1=1.  #alpha-1
    NSUB=10  #Typically set to 10
    fsave=fai
    fnew=fai_temp
    fo=fai
    f1=fnew
    i=0
    j=0  
    
    while i<3:
      da=(a1-ao)/NSUB      
      while j<NSUB:  
        a_temp=ao+da
        S_temp=S+(a_temp*dStress)                                
        p_temp=(S_temp[0]+2*S_temp[2])/3                
        q_temp=S_temp[0]-S_temp[2]
        fnew=(p_temp*p_temp)+((q_temp*q_temp)/(M*M))-(MYP*p_temp)
        j=j+1
        if fnew>FTOL:
          a1=a_temp
          if fo<-FTOL:
            f1=fnew
          break
        
        else:
          ao=a_temp
          fo=fnew 
      i=i+1   
   
    a_temp=a1-((a1-ao)*(f1/(f1-fo)))
    S_temp=S+(a_temp*dStress)                                
    p_temp=(S_temp[0]+2*S_temp[2])/3                
    q_temp=S_temp[0]-S_temp[2]	                  
  
    fnew=(p_temp*p_temp)+((q_temp*q_temp)/(M*M))-(MYP*p_temp)

    while abs(fnew) > FTOL:  
    
      if (fo/abs(fo))==(-1*(fnew/abs(fnew))):
        a1=a_temp
        f1=fnew
        if (fnew/abs(fnew))==(fsave/abs(fsave)):
          fo=fo/2
      else:
        ao=a_temp
        fo=fnew
        if (fnew/abs(fnew))==(fsave/abs(fsave)):
          f1=f1/2
      fsave=fnew
      a_temp=a1-((a1-ao)*(f1/(f1-fo)))
      S_temp=S+(a_temp*dStress)                                
      p_temp=(S_temp[0]+2*S_temp[2])/3                
      q_temp=S_temp[0]-S_temp[2]	                  
      fnew=(p_temp*p_temp)+((q_temp*q_temp)/(M*M))-(MYP*p_temp)
      

    a=a_temp
    dStress_ela= a*dStress
    De_inv=np.linalg.inv(De)
    dStrain_ela=np.matmul(De_inv,dStress_ela) 
    dStrain_Pla=dStrain-dStrain_ela            #Plastic Part Strain


    S=S+(a*dStress)                            #Stress State just before plastic behavior
    p=(S[0]+2*S[2])/3                           
    q=S[0]-S[2]	                                
    fai=(p*p)+((q*q)/(M*M))-(MYP*p)
    
    T=0
    dT=1 

    while T<1:
      dS1,dk1=pla_sub(S,De,M,eo,Lamda,k,dT,dStrain_Pla)
      S_temp=S+dS1
      dS2,dk2=pla_sub(S_temp,De,M,eo,Lamda,k,dT,dStrain_Pla)

      S_temp=S+(0.5*(dS1+dS2))
      MYP_temp=MYP+(0.5*(dk1+dk2))
      
      #Relative Error
      R1=abs(np.linalg.norm(dS2-dS1)/(2*np.linalg.norm(S_temp)))
      R2=abs((dk2-dk1)/MYP_temp)
      Rtdt=max(R1,R2,EPS)	#EPS tolerance-10e-16
      dtmin=10e-4
      while Rtdt>STOL:
        q=max(0.9*((STOL/Rtdt)**0.5),0.1)
        dT=max(q*dT,dtmin)
        dS1,dk1=pla_sub(S,De,M,eo,Lamda,k,dT,dStrain_Pla)
        dS2,dk2=pla_sub(S_temp,De,M,eo,Lamda,k,dT,dStrain_Pla)
        S_temp=S+(0.5*(dS1+dS2))
        MYP_temp=MYP+(0.5*(dk1+dk2))
        #Relative Error
        R1=abs(np.linalg.norm(dS2-dS1)/(2*np.linalg.norm(S_temp)))
        R2=abs((dk2-dk1)/MYP_temp)
        Rtdt=max(R1,R2,EPS)	#EPS tolerance-10e-16

      T=T+dT
      S=S_temp
      MYP=MYP_temp
      p=(S[0]+2*S[2])/3                                    #Mean Stress
      q=S[0]-S[2]	                                   #Deviatoric Stress
      fai=(p*p)+((q*q)/(M*M))-(MYP*p)
      
      #DRIFT_Yield_Surface
      S_temp=S
      MYP_temp=MYP
      p_temp=p               
      q_temp=q	                  
      fo=fai
      dfds=np.zeros([6,1])
     
      while abs(fo)>FTOL:
        for m in range(0,len(dfds)):
          if m<3:
            dfds[m]= ((M*M)*(((2*p_temp)-MYP_temp)/3))+(3*(S_temp[m]-p_temp))
          else:
            dfds[int(m)]=0
        Kp=((M**4)*(1+eo)*(MYP_temp)*(p_temp)*(2*p_temp-MYP_temp))/(Lamda-k)
        dfdsT=np.matrix.transpose(dfds)
        Den1 = np.matmul(dfdsT,De)
        Denominator = np.matmul(Den1,dfds) + Kp
        dlamda=fo/Denominator
        S_temp1 = S_temp - (dlamda*np.matmul(De,dfds))
        p_temp1=(S_temp1[0]+2*S_temp1[2])/3                
        q_temp1=S_temp1[0]-S_temp1[2]
        MYP_temp1=p_temp1+(q_temp1*q_temp1)/(M*M*p_temp1)
        f_temp1=(p_temp1*p_temp1)+((q_temp1*q_temp1)/(M*M))-(MYP_temp1*p_temp1)
        
        if abs(f_temp1)>abs(fo):
          dfdsT=np.matrix.transpose(dfds)
          Den = np.matmul(dfdsT,dfds)
          dlamda=fo/Den
          S_temp1 = S_temp - (dlamda*dfds)
          p_temp1=(S_temp[0]+2*S_temp[2])/3                
          q_temp1=S_temp[0]-S_temp[2]
          MYP_temp1=p_temp1+(q_temp1*q_temp1)/(M*M*p_temp1)
          f_temp1=(p_temp1*p_temp1)+((q_temp1*q_temp1)/(M*M))-(MYP_temp1*p_temp1)
         
        
        S_temp=S_temp1
        p_temp=p_temp1
        q_temp=q_temp1
        fo=f_temp1
        MYP_temp=p_temp+(q_temp*q_temp)/(M*M*p_temp)
       
      S=S_temp
      MYP=MYP_temp
      p=(S[0]+2*S[2])/3                                    #Mean Stress
      q=S[0]-S[2]	                                   #Deviatoric Stress
      fai=(p*p)+((q*q)/(M*M))-(MYP*p)

      q=max(0.9*((STOL/Rtdt)**0.5),1.1)
      dT=q*dT
      dT=max(dT,dtmin)
      dT=min(dT, 1-dT)


  Strain=Strain+dStrain               #Strain
  depsV = dStrain[0] + dStrain[1] + dStrain [2] # Increamental Volumetric Strain
  depsD = 0.6667*(dStrain[0] - dStrain[2])     # Increamental Deviatoric Strain
  SV= N-(Lamda*np.log(MYP))+(k*np.log(MYP/p))
  e=SV-1

  #Next Cycle
  pres_iter+=1
  p=(S[0]+2*S[2])/3               #Mean Stress
  q=S[0]-S[2]	                  #Deviatoric Stress
  u=p_p_i+(q/3)-p
  fai=(p*p)+((q*q)/(M*M))-(MYP*p) #Yield Surface
  VR=e		                  #Void Ratio
  epsV=epsV+depsV
  epsD=epsD+depsD
  
  f1=open('Res.txt','a')
  print("%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f" %(p,q,epsV,epsD,u,VR),file=f1)
  plt.figure(1)
  plt.plot(p,q,'o',color='black')
  plt.figure(2)
  plt.plot(Strain[0],q, 'o', color='black')

f1.close()
plt.figure(1)
plt.xlabel("Mean Effective stress (kPa)")
plt.ylabel("Deviatoric stress (kPa)")

plt.figure(2)
plt.xlabel("Axial strain (%)")
plt.ylabel("Deviatoric stress (kPa)")
plt.show()

