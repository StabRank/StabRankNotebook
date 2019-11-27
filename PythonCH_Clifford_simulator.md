```python
import numpy as np
import pandas as pd
import random
import math
import itertools
from copy import deepcopy
```


```python
def CliffordUpdate(psi,gate,qubits,n):
    """psi is a stabilizer state in CH form
    gate is 'S','CZ','CX',or 'H'
    qubits=[q] or qubits=[q1 q2] are the qubits the gate acts on
    n is the total number of qubits
    psi is a dict containing np.arrays/integers {n,s,r,Uc and p}"""
    
    if psi['n']!=n:
        raise AssertionError("Inconsistent number of qubits, n ", n, psi['n'])
    else:
        if gate=='S':
            psi['Uc']=LeftMultS(psi['Uc'],qubits,n)
        elif gate=='CZ':
            psi['Uc']=LeftMultCZ(psi['Uc'],qubits[0],qubits[1],n)
        elif gate=='CX':
            psi['Uc']=LeftMultCNOT(psi['Uc'],qubits[0],qubits[1],n)
        elif gate=='X':
            P=psi['Uc'][qubits[0],:].copy() #Pauli U_c^{\dagger} X_j U_c
            # Now conjugate by H(r)

            numy=np.mod(np.count_nonzero(P[0:n]*P[n:2*n]*psi['r']),2)#Count number of times Y-->-Y when applying H(r)
            P[2*n]=np.mod(P[2*n]+2*numy,4)#adjust phase accordingly

            inds=np.nonzero(psi['r'])[0] #Apply H(r): swap X and Z at appropriate indices

            for j in inds:
                holder=P[j].copy()
                P[j]=P[n+j].copy()
                P[n+j]=holder

            # Now compute P|s>
            v1=P[0:n]
            w1=P[n:2*n]
            alpha=np.mod(P[2*n]+np.sum(v1*w1*(1+2*psi['s']))+2*np.sum((1-v1)*w1*psi['s']),4)
            psi['p']=np.mod(psi['p']+2*alpha,8)
            psi['s']=np.mod(psi['s']+P[0:n],2)



        if gate=='H':

        #Takes a stabilizer state psi in CH form applies the single-qubit H gate to the
        #specified qubit


        #Define Paulis
        # P1=H(r) U_c ^{\dagger} X_qubit U_c H(r)
        # = [v1 | w1 | b1] (v1 is n bits, w1 is n bits, and b1 in the set 0,1,2,3)
        # P2=H(r) U_c ^{\dagger} Z_qubit U_c H(r)
        # = [v2 | w2 | b2]
        # Then the state after applying H is
        #        |phi>= w^p U_c H(r) (P1|s>+P2|s>)/sqrt(2)
            qubit=qubits[0]

            P1=psi['Uc'][qubit,:].copy()

            numy=np.mod(np.count_nonzero(P1[0:n]*P1[n:2*n]*psi['r']),2)
            P1[2*n]=np.mod(P1[2*n]+2*numy,4)

            P2=psi['Uc'][n+qubit,:].copy()
            numy2=np.mod(np.count_nonzero(P2[0:n]*P2[n:2*n]*psi['r']),2)
            P2[2*n]=np.mod(P2[2*n]+2*numy2,4)

            inds=np.nonzero(psi['r'])[0] #Apply H(r): swap X and Z at appropriate indices
            for j in inds:
                holder=P1[j].copy()
                P1[j]=P1[n+j].copy()
                P1[n+j]=holder

                holder=P2[j].copy()
                P2[j]=P2[n+j].copy()
                P2[n+j]=holder

            b1=P1[2*n]
            v1=P1[0:n]
            w1=P1[n:2*n]

            b2=P2[2*n]
            v2=P2[0:n]
            w2=P2[n:2*n]

    #         print(type(b1))
    #         print(type(v1))
    #         print(type(w1))
    #         print(type(psi['s']))
            alpha=np.mod(b1+np.sum(v1*w1*(1+2*psi['s']))+2*np.sum((1-v1)*w1*psi['s']),4)
            beta=np.mod(b2+np.sum(v2*w2*(1+2*psi['s']))+2*np.sum((1-v2)*w2*psi['s']),4)

            s1=np.mod(psi['s']+v1,2)
    #         print('s1=',s1)
            s2=np.mod(psi['s']+v2,2)
    #         print('s2=',s2)
            a=np.mod(beta-alpha,4)
    #         print('a=',a)

            y=np.mod(s1+s2,2)
    #         print(y)
    #         print("=y")

            if np.count_nonzero(y)==0:
                if a==1:
                    psi['p']=np.mod(psi['p']+2*alpha+1,8)
                else:
                    psi['p']=np.mod(psi['p']+2*alpha-1,8)

                psi['s']=s1

            if np.count_nonzero(y)>0:
                #Two cases
                klist=np.nonzero(np.mod(psi['r']+1,2)*y)[0]#Indices where r_k=0 and y_k=1
    #             print(klist)
    #             print("=klist")
                mlist=np.nonzero(psi['r']*y)[0]#Indices where r_k=1 and y_k=1
    #             print(mlist)
    #             print("=mlist")

                if klist.size!=0:# if klist is nonempty
                    k=klist[0]
                    psi['p']=np.mod(psi['p']+2*alpha+2*s1[k]*a,8)
                    psi['r'][k]=1
                    if s1[k]==0:
                        psi['s']=s1.copy()
                    else:
                        psi['s']=s2.copy()


                    for j in klist[1::]: #All except k=klist[0]
                        psi['Uc']=RightMultCNOT(psi['Uc'],k,j,n)

                    for j in mlist:
                        psi['Uc']=RightMultCZ(psi['Uc'],k,j,n)

                    if s1[k]==0:
    #                     print("a=")
    #                     print(a)
                        for j in range(1,a+1):
    #                         print("j=")
    #                         print(j)
                            psi['Uc']=RightMultS(psi['Uc'],k,n)
                    else:
                        a2=np.mod(-a,4)
    #                     print("a2=")
    #                     print(a2)
                        for j in range(1,a2+1):
    #                         print("j=")
    #                         print(j)
                            psi['Uc']=RightMultS(psi['Uc'],k,n)
    #                         print("psi['Uc']=")
    #                         print(psi['Uc'])


                if klist.size==0:# if klist is empty
                    k=mlist[0]

                    # Update phase p
                    if a==1:
                        psi['p']=np.mod(psi['p']+2*alpha+1,8)
                    if a==3:
                        psi['p']=np.mod(psi['p']+2*alpha-1,8)
                    if (a==0 or a==2):
                        psi['p']=np.mod(psi['p']+2*alpha,8)


                    B=np.nonzero(y)[0]#0-indexed in python

                    # Update Uc
                    if s1[k]==1:
                        psi['Uc']=RightMultS(psi['Uc'],k,n)#Apply Z_k gate
                        psi['Uc']=RightMultS(psi['Uc'],k,n)
                    for j in B:
    #                     print('B=',B)
    #                     print('j=',j)
    #                     print('k=',k)
                        if j != k:
                            psi['Uc']=RightMultCNOT(psi['Uc'],j,k,n)

                    if (a==1 or a==3):
                        psi['Uc']=RightMultS(psi['Uc'],k,n)


                    # Update r
    #                 print('r_in=',psi['r'])
    #                 print('k=',k)
                    if (a==0 or a==2):
    #                     print('r_to_update=',psi['r'])
    #                     print('k-th=',[k])
    #                     print('update r')
                        psi['r'][k]=np.mod(psi['r'][k]+1,2)
    #                 print('r_out=',psi['r'])


                    # Update s
                    psi['s']=s1

                    if (s1[k]==1 and (a==1 or a==2) or (s1[k]==0 and (a==0 or a==3))):
                        pass
                    else:
                        psi['s'][k]=np.mod(psi['s'][k]+1,2)
            
            
    return psi
```


```python
def AmplitudeUpdates(Pauli, psi, j, n):
    """Compute amplitude 
    amp=w^p <0|U_c^{\dagger} X_j U_c Pauli*H(r)|s>  where j is a bit index
    Here eps=0,1
    m=-n,...0
    p =0, 1,2, ...,7 
    w=exp(i*pi/4)

    Note <x|psi>=<0|U_c^{\dagger} X(x) U_c H(r)|s>w^p


    X=[x zeros(1,n) 0];
    P=PauliImage(X,psi.Uc,n);"""

    P = MultiplyPauli(psi['Uc'][j,:], Pauli, n)
    v1 = P[0:n]
    eps = 1
    if np.count_nonzero(np.mod(v1*(1-psi['r'])+psi['s']*(1-psi['r']),2)) > 0:    
        #v1 and s should only differ where r is 1
        eps = 0
        m = 0
        p = 0
        amp = 0
    else:
        p = np.mod(psi['p']-2*np.count_nonzero(P[n:2*n]*P[0:n])+2*P[2*n], 8)     
        #Number of Y paulis in P is nnz(P(n+1:2*n).*P(1:n)), each contributes -i
        m = np.count_nonzero(psi['r'])

        # Now compute sign(<v1|H(r)|s>)
        p = np.mod(p+4*np.mod(np.sum(v1*psi['s']*psi['r']), 2), 8)    
        #Look at places where r=v1=s=1. Each contributes a minus sign

        amp = 2**(-m/2)*eps*np.exp(1j*p*np.pi/4)
    return amp
        

```


```python
def Amplitude(x, psi, n):
    """Compute amplitude amp=<x|psi>=2^(-m/2)*eps* w^{p)  where x is an n-bit string and psi is a
    stabilizer state in CH form

    Here eps=0,1
    m=-n,...0
    p =0, 1,2, ...,7 
    w=exp(i*pi/4)

    Note <x|psi>=<0|U_c^{\dagger} X(x) U_c H(r)|s>w^p"""

    X = np.hstack((x,np.zeros(n,dtype=int),0))
    P = PauliImage(X, psi['Uc'], n)
    v1 = P[0:n]
    eps = 1
    if np.count_nonzero(np.mod(v1*(1-psi['r'])+psi['s']*(1-psi['r']),2)) > 0:    
        #v1 and s should only differ where r is 1
        eps = 0
        m = 0
        p = 0
        amp = 0
    else:
        p = np.mod(psi['p']-2*np.count_nonzero(P[n:2*n]*P[0:n])+2*P[2*n], 8) 
        #Number of Y paulis in P is nnz(P(n+1:2*n).*P(1:n)), each contributes -i
        m = np.count_nonzero(psi['r'])

        # Now compute sign(<v1|H(r)|s>)
        p = np.mod(p+4*np.mod(np.sum(v1*psi['s']*psi['r']), 2), 8)    
        #Look at places where r=v1=s=1. Each contributes a minus sign

        amp = 2**(-m/2)*eps*np.exp(1j*p*np.pi/4)
    return amp
```


```python
def CompBasisVector(z):
    """constructs a computational basis vector in C-H form
    |psi>= w^p U_c H(r)|s>
    where
    s is a computational basis vector
    H(r) is hadamards applied to qubits in r 
    U_c is a Clifford such that U_c|0>=|0> (stored by its tableau). 
    tableau format: first n columns are U^{\dagger} X_j U plus phase bit, last
    n columns are U^{\dagger} Z_j U plus phase bit
    p is an integer,  w=exp(i\pi/4).
    Here z vector of length n,  and basis is either 0 (z basis) or 1 (x basis)"""

    if np.ndim(z)!=1:
        z=np.squeeze(z) #Make z into a 1-d array by removing singletons

    n = np.size(z)
    psi={
        'n':n,
        's':np.asarray(z),
        #'s':z,
        'r':np.zeros(n,dtype=int),
        'Uc':np.concatenate((np.eye(2*n,dtype=int), np.zeros([2*n,1],dtype=int)), axis=1),
        'p':0,
        'c':np.exp(0*1j*np.pi)
    }
    return psi

```


```python
def XBasisVector(z):
    """constructs an X-basis vector in C-H form
    |psi>= w^p U_c H(r)|s>
    where
    s is a computational basis vector
    H(r) is hadamards applied to qubits in r 
    U_c is a Clifford such that U_c|0>=|0> (stored by its tableau). 
    tableau format: first n columns are U^{\dagger} X_j U plus phase bit, last
    n columns are U^{\dagger} Z_j U plus phase bit
    p is an integer,  w=exp(i\pi/4).
    Here z vector of length n,  and basis is either 0 (z basis) or 1 (x basis)"""

    if np.ndim(z)!=1:
        z=np.squeeze(z) #Make z into a 1-d array by removing singletons

    n = np.size(z)
    psi={
        'n':n,
        's':np.asarray(z),
        #'s':z,
        'r':np.ones(n,dtype=int),
        'Uc':np.concatenate((np.eye(2*n,dtype=int), np.zeros([2*n,1],dtype=int)), axis=1),
        'p':0,
        'c':np.exp(0*1j*np.pi)
    }
    return psi

```


```python
def RightMultCNOT(tableau, control, target, n):
    """Takes a tableau for a Clifford U and outputs the tableau for
    U*CNOT_{control,target}
     XI --> XX
     YI --> YX
     IZ --> ZZ
     IY-->  ZY 
     XX --> XI
     ZZ--> IZ
     YY--> -XZ
     XY--> YZ
     YX--> YI
     XZ--> -YY
     YZ--> XY
     ZY-->IY
    tableau format: first n columns are U^{\dagger} X_j U plus phase bit, 
    last n columns are U^{\dagger} Z_j U plus phase bit"""

    #Look at the columns of the tableau corresponding to qubits control/target

    for j in range(2*n):
        a = np.array([tableau[j, control], tableau[j, control + n], tableau[j, target], tableau[j, target + n]],dtype=int)
        b = 1*a[0] + 2*a[1] + 4*a[2] + 8*a[3]
        if b == 1:        # [1 0 0 0] % XI-->XX
            tableau[j, target] = 1
        elif b == 3:        #[1 1 0 0] % YI-->YX
            tableau[j, target] = 1
        elif b == 8:        # [0 0 0 1] %IZ-->ZZ
            tableau[j, control + n] = 1
        elif b == 12:        #[0 0 1 1] %IY-->ZY
            tableau[j, control + n] = 1
        elif b == 5:        # [1 0 1 0] % XX-->XI
            tableau[j, target] = 0
        elif b == 10:        #[0 1 0 1] %ZZ-->IZ
            tableau[j, control + n] = 0
        elif b == 15:        #[1 1 1 1] %YY-->-XZ
            tableau[j, control + n] = 0
            tableau[j, target] = 0
            tableau[j, -1] = np.mod(tableau[j,-1] + 2, 4)
        elif b == 13:        #[1 0 1 1] %XY-->YZ
            tableau[j, control + n] = 1
            tableau[j, target] = 0
        elif b == 7:        # [1 1 1 0] %YX-->YI
            tableau[j, target] = 0
        elif b == 9:        # [1 0 0 1] %XZ-->-YY
            tableau[j, control + n] = 1
            tableau[j, target] = 1
            tableau[j, -1] = np.mod(tableau[j, -1] + 2, 4)
        elif b == 11:        #[1 1 0 1] %YZ-->XY
            tableau[j, control + n] = 0
            tableau[j, target] = 1
        elif b == 14:        #[0 1 1 1] %ZY-->IY
            tableau[j, control + n] = 0
        else:
            pass
        
    return tableau
```


```python
def RightMultCZ(tableau, control, target, n):
    """ Takes a tableau for a Clifford U and outputs the tableau for            
 U*CZ_{control,target}                                                   
 XI --> XZ                                                               
 YI --> YZ                                                               
 IY-->  ZY                                                               
 IX--> ZX                                                                
                                                                         
 XX --> YY                                                               
 YY--> XX                                                                
                                                                         
 XY--> -YX                                                               
 YX--> -XY                                                               
                                                                         
 XZ--> XI                                                                
 ZX-->IX                                                                 
 YZ--> YI                                                                
 ZY-->IY                                                                 
"""

    #Look at the columns of the tableau corresponding to qubits control/target

    for j in range(2*n):
        a = np.array([tableau[j, control], tableau[j, control + n], tableau[j, target], tableau[j, target + n]],dtype=int)
        b = 1*a[0] + 2*a[1] + 4*a[2] + 8*a[3]
        if b == 1:        #  [1 0 0 0] % XI-->XZ
            tableau[j,n + target] = 1
        elif b == 3:        #[1 1 0 0] % YI-->YZ
            tableau[j,n + target] = 1
        elif b == 12:        # [0 0 1 1] %IY-->ZY
            tableau[j, control + n] = 1
        elif b == 4:        #[0 0 1 0] %IX-->ZX
            tableau[j, control + n] = 1
        elif b == 5:        # [1 0 1 0] % XX--> YY
            tableau[j, control + n] = 1
            tableau[j, target + n] = 1
        elif b == 15:        #[1 1 1 1] %YY-->XX
            tableau[j, control + n] = 0
            tableau[j, target + n] = 0
        elif b == 13:        #[1 0 1 1] %XY-->-YX
            tableau[j, control + n] = 1
            tableau[j, target + n] = 0
            tableau[j, -1] = np.mod(tableau[j, -1] + 2, 4)
        elif b == 7:        #[1 1 1 0] %YX-->-XY
            tableau[j, control + n] = 0
            tableau[j, target + n] = 1
            tableau[j, -1] = np.mod(tableau[j, -1] + 2, 4)
        elif b == 9:        #[1 0 0 1] %XZ-->XI
            tableau[j, target + n] = 0
        elif b == 6:        #[0 1 1 0] %ZX-->IX
            tableau[j, control + n] = 0
        elif b == 11:        #[1 1 0 1] %YZ-->YI
            tableau[j, target + n] = 0
        elif b == 14:        #[0 1 1 1] %ZY-->IY
            tableau[j, control + n] = 0
        else:
            pass
        
    return tableau
```


```python
def LeftMultCNOT(tableau, control, target, n):
    """Takes a tableau for a Clifford U and outputs the tableau for
    CNOT_{control,target}*U

    X_control --> X_control X_target
    Z_target --> Z_control Z_target"""

    newtableau = tableau.copy()

    
    Xprod = np.zeros(2 * n + 1,dtype=int)
    Xprod[control] = 1
    Xprod[target] = 1 #Pauli X_control*X_target

    Zprod = np.zeros( 2 * n + 1,dtype=int)
    Zprod[n + control] = 1
    Zprod[n + target] = 1 #Pauli Z_control*Z_target

    newtableau[control, :] = PauliImage(Xprod, tableau, n)
    newtableau[n + target, :] = PauliImage(Zprod, tableau, n)
    return newtableau
```


```python
def LeftMultCZ(tableau, control, target, n):
    """ Takes a tableau for a Clifford U and outputs the tableau for
    CZ_{control,target}*U

    X_control --> X_control Z_target
    X_target --> Z_control X_target"""

    newtableau = tableau.copy()

    XZ = np.zeros(2 * n + 1,dtype=int)
    XZ[control] = 1
    XZ[target + n] = 1 #Pauli X_control*Z_target

    ZX = np.zeros( 2 * n + 1,dtype=int)
    ZX[n + control] = 1
    ZX[target] = 1 #Pauli Z_control*X_target


    newtableau[control, :] = PauliImage(XZ, tableau, n) #X_control --> U^{\dagger} X_control Z_target U
    newtableau[target, :] = PauliImage(ZX, tableau, n) #X_target -->U^{\dagger} Z_control X_target U
    return newtableau
```


```python
def LeftMultS(tableau, qubit, n):
    """ Takes a stabilizer tableau for a Clifford U and left multiplies the single-qubit S gate to the
    specified qubit : U--> S_qubit*U

    S^{\dagger} X S=-Y
    S^{\dagger} Z S=Z 

    So must replace row qubit of tableau with image of -Y under U"""
    newtableau = tableau.copy()
    minusYj = np.zeros(2 * n + 1,dtype=int)
    minusYj[2*n] = 2 #phase is -1
    minusYj[qubit] = 1
#     print(minusYj)
#     print(type(_))
#     print(qubit)
    minusYj[qubit[0] + n] = 1
    newtableau[qubit, :] = PauliImage(minusYj, tableau, n)# outputs -U^{\dagger} Y_qubit U^{\dagger}
    return newtableau
```


```python
def RightMultS(tableau, qubit, n):
    """ Takes a stabilizer tableau for a Clifford U and left multiplies the single-qubit S gate to the
    specified qubit : U--> U*S

    Look at the columns of the tableau corresponding to qubit"""  



    for j in range(2*n):
        if tableau[j, qubit] == 1:        #X_qubit--> -Y_qubit; %Y_qubit--> X_qubit
#             print("j=")
#             print(j)
            tableau[j, qubit + n] = np.mod(tableau[j, qubit + n] + 1, 2)
#             print("tableau[j, qubit]=")
#             print(tableau[j, qubit])
            tableau[j, 2*n] = np.mod(tableau[j, 2*n] + 2 * tableau[j, qubit + n], 4)        #Flip phase if X-->-Y
#             print("tableau[j, 2*n]=")
#             print(tableau[j, 2*n])
            
    return tableau
```


```python
def PauliImage(P, tableau, n):
    """Takes as input a tableau for a Clifford U_c and a Pauli P and outputs
    the pauli Uc^{\dagger} P U_c"""


    Pout = np.zeros(2 * n + 1,dtype=int)
    Pout[-1] = P[-1].copy()
    for j in range(n):
        if P[j] == 1:        #Multiply by U_c X_j U_c
            Pout = MultiplyPauli(Pout, tableau[j,:], n)

        if P[n + j] == 1:        #Multiply by U_c Z_j U_c
            Pout = MultiplyPauli(Pout, tableau[n + j,:], n)

        if P[j] == 1 and P[n + j] == 1:        #Correct phase X_j*Z_j*(+i)=Y_j
            Pout[-1] = np.mod(Pout[-1] + 1, 4)
    return Pout
```


```python
def MultiplyPauli(P1, P2, n):
    #Multiplies Paulis P3=P1*P2 as in Eq. 15.17 of Kitaev/Vyalyi/Shen book

    P3 = np.mod(P1 + P2, 2)#Gets everything correct except phase bit
    tau1 = 0
    tau2 = 0
    tau3 = 0
    kap = 0

    for j in range(n):
        tau1 = tau1 + P1[j] * P1[n + j]
        tau2 = tau2 + P2[j] * P2[n + j]
        tau3 = tau3 + P3[j] * P3[n + j]
        kap = np.mod(kap + P2[j] * P1[j + n], 2)

    P3[-1] = np.mod(P1[-1] + P2[-1] + tau1 + tau2 - tau3 + 2 * kap, 4)
    return P3

```


```python
def bin_array(num, m):
    """Convert a positive integer num into an m-bit bit vector"""
    return np.array(list(np.binary_repr(num).zfill(m))).astype(np.int8)
```


```python
def Indicator(nparray, *length):
    if not length:
        n_values = np.max(np.asarray(nparray))+1
        mat=np.eye(n_values,dtype=int)[nparray]
    else:
        n_values = length[0]+1
        mat=np.eye(n_values,dtype=int)[nparray]
    return mat  
```


```python
def U_H(psi,i,n):
    i+=1 #zero-indexing
    h=np.asarray([[1,1],[1,-1]])*2**-0.5
    H=np.kron(np.kron(np.eye(2**(i-1)),h),np.eye(2**(n-i)))
    #print(H)
    return H@psi

def U_X(psi,i,n):
    i+=1 #zero-indexing
    x=np.asarray([[0,1],[1,0]])
    X=np.kron(np.kron(np.eye(2**(i-1)),x),np.eye(2**(n-i)))
    #print(X)
    return X@psi

def U_S(psi,i,n):
    i+=1 #zero-indexing
    s=np.asarray([[1,0],[0,1j]],dtype=complex)
    S=np.kron(np.kron(np.eye(2**(i-1)),s),np.eye(2**(n-i)))
    #print(S)
    return S@psi


def U_CX(psi,i,j,n):
    i+=1 #zero-indexing
    j+=1 #zero-indexing
    # Makes local projector P_ij
    swap=np.asarray([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    cnot=np.asarray([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    if (j<i):
        cnot=swap@cnot@swap
    
    i2=min([i,j])    
    j2=max([i,j])
    j=j2
    i=i2
    
    CNOT=np.kron(np.kron(np.eye(2**(i-1)),cnot),np.eye(2**(n-i-1)))
    
    for k in range(i+1,j):
        swk=np.kron(np.kron(np.eye(2**(k-1)),swap),np.eye(2**(n-k-1)))
        CNOT=swk@CNOT@swk
    
    return CNOT@psi

def U_CZ(psi,i,j,n):
    i+=1 #zero-indexing
    j+=1 #zero-indexing
    # Makes local projector P_ij
    swap=np.asarray([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    cz=np.asarray([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])
    if (j<i):
        cz=swap@cz@swap
    
    i2=min([i,j])    
    j2=max([i,j])
    j=j2
    i=i2
#     print('i,j=',i,j)
    
    CZ=np.kron(np.kron(np.eye(2**(i-1)),cz),np.eye(2**(n-i-1)))
    
#     print(list(range(i+1,j)))
    for k in range(i+1,j):
        swk=np.kron(np.kron(np.eye(2**(k-1)),swap),np.eye(2**(n-k-1)))
        CZ=swk@CZ@swk
#     print(CZ)
    
    return CZ@psi
```


```python
def UnitaryUpdate(psi,gate,qubits,n):
    """psi is a stabilizer state in CH form
    gate is 'S','CZ','CX',or 'H'
    qubits=[q] or qubits=[q1 q2] are the qubits the gate acts on
    n is the total number of qubits
    psi is a dict containing np.arrays/integers {n,s,r,Uc and p}"""
    if gate=='S':
        psi=U_S(psi,qubits,n)
    elif gate=='CZ':
        psi=U_CZ(psi,qubits[0],qubits[1],n)
    elif gate=='CX':
        psi=U_CX(psi,qubits[0],qubits[1],n)
    elif gate=='X':
        psi=U_X(psi,qubits,n)
    elif gate=='H':
        psi=U_H(psi,qubits,n)        
    
    return psi
```


```python
def compUnitaryUpdate(psi,gate,qubits,n):
    """psi is a stabilizer state in CH form
    gate is 'S','CZ','CX',or 'H'
    qubits=[q] or qubits=[q1 q2] are the qubits the gate acts on
    n is the total number of qubits
    psi is a dict containing np.arrays/integers {n,s,r,Uc and p}"""
    if gate=='S':
        psi=U_S(psi,qubits[0],n)
    elif gate=='CZ':
        psi=U_CZ(psi,qubits[0],qubits[1],n)
    elif gate=='CX':
        psi=U_CX(psi,qubits[0],qubits[1],n)
    elif gate=='X':
        psi=U_X(psi,qubits[0],n)
    elif gate=='H':
        psi=U_H(psi,qubits[0],n)        
    
    return psi
```


```python
def fullvec(psi):
    "Write out psi in the computational basis using normal lexic ordering e.g. 000,001,...,111"
    
    n=psi['n']
    psivec=np.zeros((2**n,1),dtype=complex)
    for j in range(2**n):
        psivec[j]=Amplitude(bin_array(j, n),psi,n)
    
    if 'c' in psi:
        psivec=psi['c']*psivec
        
    return psivec
    
```


```python
def normalize(v):
    return v / np.linalg.norm(v)
```


```python
n=3
#Make |100>. 
# Here the leftmost bit is the least significant
# The qubits are labeled 0,1,2,3,...n-1 from left to right
psi=CompBasisVector([1, 0, 0]); 


# Apply CNOT_{12} gate (1 is target, 2 is control)

psi=CliffordUpdate(psi,'CX',[0, 1],n);
# Apply CNOT_{23} gate 
psi=CliffordUpdate(psi,'CX',[1, 2],n);
#State is now |111>

#Apply Hadamards to qubits 2,3
 
psi=CliffordUpdate(psi,'H',[1],n);   
psi=CliffordUpdate(psi,'H',[2],n); 

# State is |1-->

#Apply CZ_{12} and CZ_{13}
psi=CliffordUpdate(psi,'CZ',[0, 1],n);
psi=CliffordUpdate(psi,'CZ',[0, 2],n);

# State is |1++>

# Apply S gate to 1st qubit
   
psi=CliffordUpdate(psi,'S',[0],n);   
# State is i*|1++>



#Note that the ordering of bits is opposite to the tensor product
#ordering. The state q is equal to 
v=(1/2)*1j*np.kron(np.kron(np.asarray([[0],[1]]),np.asarray([[1],[1]])),np.asarray([[1],[1]]))
print(v)
print(fullvec(psi))
                         
#Check q-v=0
# norm(q-v)
```


```python
def update_and_compare(psi_in,phi_in,gtype,gqubits,n):
    psi = CliffordUpdate(psi_in, gtype, gqubits, n)
    phi = compUnitaryUpdate(phi_in, gtype, gqubits, n)
    return psi, phi, np.linalg.norm(fullvec(psi) - phi)
    
```


```python
def displayboth(psi,phi):
    print(np.linalg.norm(fullvec(psi) - phi) )
    print(np.squeeze(np.asarray([fullvec(psi),phi]),axis=2).T)
```


```python
def randomstabilizer(n):
    psi=XBasisVector(np.zeros((n),dtype=int)); # Psi will be the state computed by the CH Clifford simulator
    phi=np.ones((2**n,1),dtype=complex)*2**(-n/2)

    m=100;   #Number of randomly chosen gates
    r=np.random.randint(5, size=m);
    qubits=[];
    nrm=0;

    for j in range(m):
        if r[j] == 0:
            q = np.random.randint(n)
#             print('S',[q])
            qubits = np.append(qubits, q)
            psi, phi, nrm = update_and_compare(psi,phi,'S',[q],n)
#             print('S',[q],nrm)

        if r[j] == 1 or r[j] == 2:
            collision = 0
            while collision == 0:
                q1 = np.random.randint(n)
                q2 = np.random.randint(n)
                if q1 != q2:
                    collision = 1

            qubits = np.append(qubits, [q1, q2])
    #         print(qubits)
            if r[j] == 1:
#                 print('CZ', [q1, q2])
                psi, phi, nrm = update_and_compare(psi,phi,'CZ', [q1, q2],n)
#                 print('CZ', [q1, q2],nrm)

            if r[j] == 2:
#                 print('CX', [q1, q2])
                psi, phi, nrm = update_and_compare(psi,phi,'CX', [q1, q2],n)
#                 print('CX', [q1, q2],nrm)

        if r[j] == 3:
            q = np.random.randint(n)
#             print('H',[q])
            qubits = np.append(qubits, q)
            psi, phi, nrm = update_and_compare(psi,phi,'H',[q],n)
#             print('H',[q],nrm)

        if r[j] == 4:
            q = np.random.randint(n)
#             print('X',[q])
            qubits = np.append(qubits, q)
            psi, phi, nrm = update_and_compare(psi,phi,'X',[q],n)
#             print('X',[q],nrm)

    q=fullvec(psi)

    return psi, phi
```


```python
psi, phi = randomstabilizer(4)
displayboth(psi,phi)
```


```python
def EquatorialA(n):
    "Randomly generate the A matrix corresponding to |phi_A> in Eq 61 "
    
    A=np.zeros((n,n),dtype=np.int8)
    offdiag=np.random.randint(2,size=int(n*(n-1)/2))
    A[np.triu_indices(n, 1)]=offdiag
    A=A+A.T
    np.fill_diagonal(A,np.random.randint(4,size=n))
    return A
```


```python
def ind2vec(ind, N=None):
    ind = np.asarray(ind)
    if N is None: 
        N = ind.max() + 1
    return (np.arange(N) == ind[:,None]).astype(int)
```


```python
def binvec2dec(vec):
    return np.dot((2**np.arange(vec.shape[0], dtype = np.uint64)[::-1]),vec).astype(np.int)
```


```python
def Equatorialfullvec(A):
    """Write out psi in the computational basis using normal lexic ordering e.g. 000,001,...,111. 
    See Sec IV C of BBCCGH"""  
    
    n=A.shape[0]
    psivec=np.zeros((2**n,1),dtype=complex)
    for j in range(2**n):
        x=bin_array(j,n)
        psivec=psivec+1j**((x@A)@x.T)*ind2vec([binvec2dec(x)],2**n).T
    psivec=2**(-n/2)*psivec
        
    return psivec
```


```python
def InnerPsiA(psi,A):
    """Implementation of <phi|phi_A> as described in Lemma 3, Sec IV C of BBCCGH"""  
    
    n=psi['n']
    r=psi['r']
    s=psi['s']
    gamma=psi['Uc'][0:n,-1]
    F=psi['Uc'][0:n,0:n]
    M=psi['Uc'][0:n,n:2*n]
    G=psi['Uc'][n:2*n,n:2*n]
    J=np.mod(M@F.T+np.diag(gamma),4)
    mask=~np.eye(J.shape[0],dtype=bool)
    J[mask]=np.mod(J[mask],2)
    K=G.T@(A+J)@G
    wtr=np.count_nonzero(r)
    term1=2**(-(n+wtr)/2)
    term2=1j**((s@K)@s.T)
    term3=(-1)**np.dot(s,r)
    phase=np.exp(-psi['p']*1j*np.pi/4)*psi['c'].conj()
    locs=np.nonzero(r)[0]
    sumtot=0
    for j in range(2**(wtr)):
        x=np.zeros(n,dtype=int)
        x[locs]=bin_array(j, wtr)#n-bit strings x satisfying x_j <= r_j
        sumtot=sumtot+1j**((x@K)@x.T+2*x@(s+s@K))
    return phase*term1*term2*term3*sumtot
```


```python
def checkInnerPsiA(psi,A):
    aa=abs(InnerPsiA(psi,A))
    bb= abs(np.vdot(fullvec(psi),Equatorialfullvec(A)))
    return print([aa, bb, aa-bb])
```


```python
def etaA(SUBSET,Amat):
    n=Amat.shape[0]
    olap=sum(list(map(lambda x: InnerPsiA(SUBSET[x],Amat), range(len(SUBSET)))))
    return 2**n*abs(olap)**2
```


```python
def eta_estimate(SUBSET,n,Lsamples):
    return np.mean([etaA(SUBSET,EquatorialA(n)) for j in range(Lsamples)])
```


```python
def explicitetaA(SUBSET,Amat):
    tot=0
    for j in range(len(SUBSET)):
        tot=tot+2**6*np.vdot(fullvec(SUBSET[j]),Equatorialfullvec(Amat))
    return tot
```


```python
def HiddenShiftCircuitBuilderCH(numtof,numqub,*hidstring):
    """Build the circuit for implimenting Hidden Shift for Bent functions"""
    n=numqub
    if not hidstring:
        s=np.random.randint(2, size=n)
    else:
        s=np.asarray(hidstring)[0]
    
    q3=np.random.choice(int(n/2), 3, replace=False).copy()
    Og=pd.DataFrame([{'H':[q3[2].copy()]},{'Toff':q3.copy()},{'H':[q3[2].copy()]}])
    Og2=pd.DataFrame([{'H':[int(n/2)+q3[2].copy()]},{'Toff':int(n/2)+q3.copy()},{'H':[int(n/2)+q3[2].copy()]}])
    for j in range(numtof-1):
        for i in range(10):# The number of random Clifford gates Z or C-Z
            r=np.random.randint(2)
            if r==1:
                #print('random Z')
                q=np.random.randint(n/2)
                Og=Og.append(pd.DataFrame([{'S':[q]},{'S':[q]}]),ignore_index=True) # Z gate applied to qubit q
                Og2=Og2.append(pd.DataFrame([{'S':[int(n/2)+q]},{'S':[int(n/2)+q]}]),ignore_index=True)
            elif r==0:
                #print('random CZ')
                q2=np.random.choice(int(60/2), 2, replace=False).copy()
                Og=Og.append(pd.DataFrame([{'CZ':q2}]),ignore_index=True) # C-Z gate between qubits in q2
                Og2=Og2.append(pd.DataFrame([{'CZ':int(n/2)+q2}]),ignore_index=True)
        q3=np.random.choice(int(60/2), 3, replace=False).copy()
        Og=Og.append(pd.DataFrame([{'H':[q3[2]]},{'Toff':q3},{'H':[q3[2]]}]),ignore_index=True)
        Og2=Og2.append(pd.DataFrame([{'H':[int(n/2)+q3[2]]},{'Toff':int(n/2)+q3},{'H':[int(n/2)+q3[1]]}]),ignore_index=True)

    Of=Og.copy()
    Oftil=pd.DataFrame()
    for j in range(int(n/2)):
        Of=Of.append(pd.DataFrame([{'CZ':[j,j+int(n/2)]}]),ignore_index=True)
        Oftil=Oftil.append(pd.DataFrame([{'CZ':[j,j+int(n/2)]}]),ignore_index=True)
    Oftil=Oftil.append(Og2)

    # This portion needs checking
    Xs=pd.DataFrame()
    FT=pd.DataFrame()
    for i in range(n):
        FT=FT.append(pd.DataFrame([{'H':[i]}]))
        if s[i]==1:
            Xs=Xs.append(pd.DataFrame(np.asarray([[[i],np.nan],[np.nan,[i]],[np.nan,[i]],[[i],np.nan]],dtype=object),columns=['H','S']),ignore_index=True)

    U=FT.append([Of,Xs,FT,Oftil,FT],ignore_index=True)
    return U, s 
```


```python
def SamplesToTake(circ,delta):
    numtof=0
    numt=0
    if 'CCZ' in circ.columns:
        numtof=numtof+circ.count()['CCZ']
    if 'Toff' in circ.columns:
        numtof=numtof+circ.count()['Toff']
    if 'T' in circ.columns:
        numt=numt+circ.count()['T']
    
    print('numtof',numtof)
    print('numt',numt)
    c1T=np.cos(np.pi/8)**(-numt)
    c1Toff=(3/4)**(-numtof)
    ksamps=int(round((c1T*c1Toff/delta)**2))
    #2**(-m/2)*eps*np.exp(1j*p*np.pi/4)
    return ksamps
```


```python
def CreateSubsetStates(circ,ksamps,psi_in):
    cliffordgates={'H','S','CZ','CX'}


    pathind=0

    SUBSET=[]
    for path in range(ksamps):
        #print('path=',path)
        #psi=XBasisVector([0,0,0])#create the |+++> state
        psi=deepcopy(psi_in)
        for ind, row in circ.iterrows():
            tmp=row.dropna()
            gtype=tmp.index[0]
#             print(gtype)
            gqubits=np.asarray(tmp.values[0])
            if gtype in cliffordgates:
                psi=CliffordUpdate(psi, gtype, gqubits, n)
            elif gtype=='CCZ':
                x=random.choices(list(range(0,8)), weights=np.ones(8)/6, k=1)[0]
                #print('CCZ','x',x)
                if x==0: # Identity
                    pass
                elif x==1: # CZ_{12}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]],n)
                elif x==2: # CZ_{13}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]],n)
                elif x==3: # CZ_{12,13}Z_{1}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                elif x==4: # CZ_{23}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                elif x==5: # CZ_{23,12}Z_{2}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                elif x==6: # CZ_{13,23}Z_{3}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                elif x==7: # -CZ_{12,13,23}Z_{1,2,3}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                    psi['p']=np.mod(psi['p']+4,8) # This inserts the minus sign

            elif gtype=='Toff':
                psi=CliffordUpdate(psi, 'H', gqubits[[2]], n)
                x=random.choices(list(range(0,8)), weights=np.ones(8)/6, k=1)[0]
#                 x=feynmanpaths[pathind]
#                 print('Toff','x',x)
                pathind+=1
                if x==0: # Identity
                    pass
                elif x==1: # CZ_{12}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]],n)
                elif x==2: # CZ_{13}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]],n)
                elif x==3: # CZ_{12,13}Z_{1}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                elif x==4: # CZ_{23}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                elif x==5: # CZ_{23,12}Z_{2}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                elif x==6: # CZ_{13,23}Z_{3}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                elif x==7: # -CZ_{12,13,23}Z_{1,2,3}
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,1]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[0,2]], n)
                    psi=CliffordUpdate(psi, 'CZ', gqubits[[1,2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[0]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[1]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                    psi=CliffordUpdate(psi, 'S', gqubits[[2]], n)
                    psi['p']=np.mod(psi['p']+4,8) # This inserts the minus sign
#                     psi['c']=(-1)*psi['c']
                psi=CliffordUpdate(psi, 'H', gqubits[[2]], n)

            elif gtype=='T':
                x=random.choices(list(range(0,2)), weights=np.ones(2)/2, k=1)[0]
                # T = c_0 Id + c_1 S where 
                # c_0 = 0.5*(1+1j)*(np.exp(1j*math.pi/4)-1j)
                # c_1 = -0.5*(1+1j)*(np.exp(1j*math.pi/4)-1)
                phi_0=np.angle(0.5*(1+1j)*(np.exp(1j*math.pi/4)-1j))
                phi_1=np.angle(-0.5*(1+1j)*(np.exp(1j*math.pi/4)-1))
                #print('T','x',x)
                if x==0: # Identity
                    pass
                    psi['c']=np.exp(1j*phi_0)*psi['c']
                elif x==1: # S
                    psi=CliffordUpdate(psi, 'S', gqubits, n)
                    psi['c']=np.exp(1j*phi_1)*psi['c']

        SUBSET.append(psi)

    return SUBSET
```


```python
(circ,shift)=HiddenShiftCircuitBuilderCH(1,6,[0,0,0,1,1,1])
print(shift)
# (circ,shift)=HiddenShiftCircuitBuilderCH(1,6)
# print(shift)

SamplesToTake(circ,0.1)
```


```python
n=shift.shape[0]
def SumOverCliffords(circ):
    print('n=',n)
    ksamps=SamplesToTake(circ,0.1)
    print(ksamps)
    return CreateSubsetStates(circ,ksamps,CompBasisVector(np.asarray(np.zeros(n),dtype=int)))
```


```python
n=3
circ=pd.DataFrame([{'CCZ':[0,1,2]}])
ksamps=SamplesToTake(circ,0.1)
print('ksamps',ksamps)

SUBSET=CreateSubsetStates(circ,ksamps,XBasisVector([0,0,0]))

print(normalize(sum(list(map(fullvec,SUBSET)))))
```


```python
circ=pd.DataFrame([{'Toff':[0,1,2]}])
ksamps=SamplesToTake(circ,0.1)
print('ksamps',ksamps)

SUBSET=CreateSubsetStates(circ,ksamps,XBasisVector([0,0,0]))
# SUBSET, phiSUBSET =DebugSubsetStates(circ,ksamps,CompBasisVector([1,1,1]))

print(normalize(sum(list(map(fullvec,SUBSET)))))
```


```python
n=3
circ=pd.DataFrame([{'CCZ':[0,1,2]},{'Toff':[0,1,2]},{'T':[0]},{'T':[1]},{'T':[2]}])
ksamps=SamplesToTake(circ,0.05)
print('ksamps',ksamps)

SUBSET=CreateSubsetStates(circ,ksamps,XBasisVector([0,0,0]))
# SUBSET, phiSUBSET =DebugSubsetStates(circ,ksamps,CompBasisVector([1,1,1]))

print(normalize(sum(list(map(fullvec,SUBSET)))))
```


```python
(circ,shift)=HiddenShiftCircuitBuilderCH(1,6,[0,0,0,1,1,1])
print(shift)
# (circ,shift)=HiddenShiftCircuitBuilderCH(1,6)
# print(shift)

SamplesToTake(circ,0.1)
```


```python
n=shift.shape[0]
def SumOverCliffords(circ):
    print('n=',n)
    ksamps=SamplesToTake(circ,0.1)
    print(ksamps)
    return CreateSubsetStates(circ,ksamps,CompBasisVector(np.asarray(np.zeros(n),dtype=int)))
```


```python
SUBSET=SumOverCliffords(circ)
unnorm=sum(list(map(fullvec,SUBSET)))
answer=normalize(unnorm)
print(np.sum(np.round(answer)-fullvec(CompBasisVector(shift))))
```


```python
A6= EquatorialA(6)
print(A6)
phiA=Equatorialfullvec(A6)
print('<phi_A|phi_A>=',np.vdot(phiA,phiA))
print('<phi_A|psi>=',np.vdot(unnorm,phiA))
print('<psi|psi>=nrm_sq=',np.vdot(unnorm,unnorm))
print('eta_estimate=',eta_estimate(SUBSET,6,100))
```


```python

    
```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```
