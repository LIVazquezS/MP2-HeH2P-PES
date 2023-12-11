'''

An implementation of python to obtain the coefficients of equation 2 and 3 in Phys.Chem.Chem.Phys. 2019, 21, 24976
It uses the packages lmfit
Author: Luis Itza Vazquez-Salazar
email: luisitza.vazquezsalazar@unibas.ch
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Model, Parameter

def r_hat(rab,ri,req):
    return rab + ri*np.exp(-(rab-req))

def vlong(ad,aq,ao,bddq,gd,ri,req,q,rab):
    rhat = r_hat(rab,ri,req)
    t1 = (ad*q**2)/(2*rhat**4)
    t2 = (aq*q**2)/(2*rhat**6)
    t3 = (ao*q**2)/(2*rhat**8)
    t4= (bddq*q**3)/(6*rhat**7)
    t5 = (gd*q**4)/(24*rhat**8)
    vl = -t1-t2-t3-t4-t5
    return vl

def first_term(c0,aab,rab):
    factor = np.exp(-(aab*rab))
    fst = (c0*factor)/rab
    return fst

def rho(bab,rab):
    e_term = np.exp(-bab*rab)
    rho_ab = rab*e_term
    return rho_ab

def sum_term(ci,rab,bab):
    sum = 0
    for i in range(len(ci)):
        s = ci[i]*(rho(bab,rab))**(i+1)
        sum += s
    return sum

def vab(rab,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,He=True):
    # Order of parameters
    # C vector contains all parameters to be optimized
    # C[0] is C0 for first term
    # C[1] is aab
    # from C[2] to C[11] there are the 10 terms of sum
    # C[12] is bab


    # Some scalar parameters

    if He:
        # For He
        ad = 1.384
        aq = 2.275
        ao = 10.620
        bddq = 20.41
        gd = 37.56
        req = 1.4588686
        ri = 8.0
        q = 1.0
    else:
        # For H2+, from table 1
       ad = 4.5
       aq = 15.0
       ao = 131.25
       bddq = 159.75
       gd = 1333.125
       ri = 10.0
       req = 1.9973041778
       q = 1.0

    c_array = [c2,c3,c4,c5,c6,c7,c8,c9,c10,c11]
    vab2 = first_term(c0,c1,rab) + sum_term(c_array,rab,c12)+vlong(ad,aq,ao,bddq,gd,ri,req,q,rab)
    return vab2

def residuals(c,rab,eab):
    model = vab(rab,c)
    return model - eab

def validation(dict,rab):
    e = vab(rab,dict['c0'],dict['c1'],dict['c2'],
            dict['c3'],dict['c4'],dict['c5'],dict['c6'],dict['c7'],dict['c8'],dict['c9']
            ,dict['c10'],dict['c11'],dict['c12'])
    return e

def main(data,He=True):
    df = pd.read_csv(data)
    rab = df['rab'].to_numpy()
    #Zero for He is the atomization energy
    if He:
        eab = df['E'].to_numpy() + 2.898161153989
    else:
        #Zero in case of H2+
        eab = df['E'].to_numpy() + 0.499994784584
    model1 = Model(vab)
    #Initial guess of the parameters, taken from the paper for FCI.
    if He:
        db_p = [3.2233398865576910, 2.8140238576362075, -3.3885062843889191, 86.578369632144941,
                -2929.2064910602976, 73876.513409500127, -1280675.3543930338, 14775779.204776313,
                -110729495.67500649, 512524888.74876255, -1315519328.5011687, 1401896260.3220801,
                2.1737434370604722]
    else:
        db_p = [1.0185817913777047, 1.6426016158795942, -0.7368835063554769,5.1706704100481042,
             -83.346600494014567,874.72395010350556,-5621.3191308852929,23054.48904746825,
             -60955.753058022863,100818.87103555651,-94989.664045708429,38921.331029687601,
               0.99809657513307815]

    params = model1.make_params()
    for i,j in enumerate(model1.param_names):
        params[j] = Parameter(name=j,value=db_p[i],vary=True)
    params['c0'] = Parameter(name='c0',min=0.000001,value=3.2233398865576910,vary=True)
    result = model1.fit(eab,params,rab=rab)
    print(result.fit_report())
    coeff = result.best_values

    e_model = validation(coeff,rab)
    diff = eab-e_model
    rmse = np.sqrt(np.mean(np.square(e_model-eab)))
    print(rmse)
    print(validation(coeff,2.1))
    fig,ax = plt.subplots(1,2)
    ax[0].scatter(eab,e_model)
    # ax[0].set_xlim(-2.90,-2.898)
    # ax[0].set_ylim(-2.90,-2.898)
    ax[1].plot(rab,eab,'r')
    ax[1].scatter(rab,result.eval())
    ax[1].set_xlim(0,10)
    plt.show()

data = 'h2p.csv'
main(data,He=False)

