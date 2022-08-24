# Author: Allan Moreira de Carvalho
# Created: 23.08.2022

from solvers import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import pickle

def Area(x, xt):
    # Area function as defined in Yeom 
    ret = np.zeros_like(x)
    for i,xi in enumerate(x):
        if xi < xt:
            ret[i] = 2.5 + 3.0*(x[i]/xt - 1.5)*(x[i]/xt)**2.0
        elif xi >= xt:
            ret[i] = 3.5 - (x[i]/xt)*(6.0 - 4.5*(x[i]/xt) + (x[i]/xt)**2.0)
    return ret

def interpolate(yp, N_pts):
    i = np.arange(len(yp))
    interp_i = np.linspace(0, i.max(), N_pts)

    yi = interp1d(i, yp, kind='cubic')(interp_i)

    return yi

def Fspline(sp):
    # Setup nozzle shape 
    nozzleL = [0.0,10.0]
    #xt = 5.0

    # Boundary conditions
    T0in = 291.3
    pr = 0.65
    p0in = 104074.0
    Min = 0.01
    
    # Fluid properties
    gamma = 1.4
    R = 287.0
    
    # eulerQ1D setup
    itmaxQ1D = 100000
    itprintQ1D = 1000
    CFLQ1D = 0.1
    tschemeQ1D = 'RK4'
    fschemeQ1D = 'AUSM'
    dttypeQ1D = 'Global'
    dimQ1D = 'Dimensionless'

    noz = nozzle()

    noz.setBC(p0in, T0in, Min, pr*p0in)
    noz.setFluid(R, gamma)

    # Q1D setup
    noz.setX(np.linspace(nozzleL[0], nozzleL[1], nQ1D))

    s = interpolate(sp, nQ1D)
    noz.setS(s)
    
    noz.setQ1DSolver(itmaxQ1D, itprintQ1D, CFLQ1D, tolQ1D, tschemeQ1D, fschemeQ1D, dttypeQ1D, dimQ1D)
    
    noz.setupQ1D(workdir, 'setupQ1D.txt')
    
    noz.setSolverQ1D('./eulerQ1D')
    noz.solveQ1D()

    p = noz.getQ1D('outputs/p.txt')

    return p, noz

def Fxt(xt):
    # Setup nozzle shape 
    nozzleL = [0.0,10.0]
    #xt = 5.0

    # Boundary conditions
    T0in = 291.3
    pr = 0.65
    p0in = 104074.0
    Min = 0.01
    
    # Fluid properties
    gamma = 1.4
    R = 287.0
    
    # eulerQ1D setup
    itmaxQ1D = 100000
    itprintQ1D = 1000
    CFLQ1D = 0.1
    tschemeQ1D = 'RK4'
    fschemeQ1D = 'AUSM'
    dttypeQ1D = 'Global'
    dimQ1D = 'Dimensionless'

    noz = nozzle()

    noz.setBC(p0in, T0in, Min, pr*p0in)
    noz.setFluid(R, gamma)

    # Q1D setup
    noz.setX(np.linspace(nozzleL[0], nozzleL[1], nQ1D))
    noz.setS(Area(noz.xn, xt))
    
    noz.setQ1DSolver(itmaxQ1D, itprintQ1D, CFLQ1D, tolQ1D, tschemeQ1D, fschemeQ1D, dttypeQ1D, dimQ1D)
    
    noz.setupQ1D(workdir, 'setupQ1D.txt')
    
    noz.setSolverQ1D('./eulerQ1D')
    noz.solveQ1D()

    p = noz.getQ1D('outputs/p.txt')

    return p, noz

def J(p, pt):
    N = np.size(pt)
    #return np.sum((p-pt)**2.0)
    #return np.sqrt(np.sum(((p-pt)/pt)**2.0))
    return 1/N*np.sqrt(np.sum((p-pt)**2))
    #return np.sum(((p-pt)/pt)**2.0)
    
def objF(xt, p_target):
    if interiorOnly:
        xt = np.append(np.append(st[0],xt),st[-1])
    
    p_guess, noz_guess = F(xt)
    J_ = J(p_guess, p_target)
    return J_ 

def callB1(xk):
    if plot:
        plt.figure()
        s = interpolate(xk, nQ1D) 
        plt.plot(xi,s)
        plt.plot(xi,st)
        plt.show()
    xhist.append(xk)

def callB(xk):
    if plot:
        plt.figure()
        x=np.linspace(0,10, nQ1D)
        if np.size(xk)==1:
            s = Area(x,xk[0])
        else:
            if interiorOnly:
                xk = np.append(np.append(st[0],xk),st[-1])
            s = interpolate(xk, nQ1D) 
        plt.plot(x,s)
        plt.plot(x,st)
        plt.show()
    xhist.append(xk)

def callBres(xk, res):
    if plot:
        plt.figure()
        x=np.linspace(0,10,nQ1D)
        if np.size(xk)==1:
            s = Area(np.linspace(0,10,nQ1D),xk[0])
        else:
            if interiorOnly:
                xk = np.append(np.append(st[0],xk),st[-1])
            s = interpolate(xk, nQ1D) 
        plt.plot(x,s)
        plt.plot(x,st)
        plt.show()
    xhist.append(xk)

def getSpacedElements(array, numElems = 4):
    out = array[np.round(np.linspace(0, len(array)-1, numElems)).astype(int)]
    return out

def opt(methods, p_target, init, Fun, bounds=None, plots=False, savedict=None, filename='DefaultNameXt.pickle'):
    global F
    global plot
    
    plot = plots
    F = Fun
    global xhist
    if savedict is not None:
        res = savedict
    else:
        res = {}

    for met in methods:
        call = callB
        xhist = []
        if met in ['trust-constr']:
            call = callBres
            res[met] = minimize(fun=objF, x0=init, args=(p_target), bounds=bounds, method=met, options={'disp': True}, callback=call, tol=tol)
        else:
            res[met] = minimize(fun=objF, x0=init, args=(p_target), bounds=bounds, method=met, options={'disp': True}, callback=call, tol=tol)
        res[met]['xhist'] = xhist
        
        if savedict is not None:
            file_to_store = open(filename,'wb')
            pickle.dump(res, file_to_store)
            file_to_store.close()

        if plot:
            print(res[met])
            p_opt, noz_opt = F(res[met].x)
            p_ini, noz_ini = F(init)
            
            x = noz_ini.getQ1D('outputs/x.txt')

            plt.plot(x, p_target, label='target')
            plt.plot(x, p_ini, label='initial guess')
            plt.plot(x, p_opt, ls=':',label='optimization result')

            plt.legend()
            plt.ylabel('pressure')
            plt.xlabel('x')
            plt.title(met)
            plt.show()
        
    return res

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--polynomial', help='polynomial parametrization', default=False)   
    parser.add_argument('-s','--spline', help='spline parametrization', default=False)   
    parser.add_argument('-sbi','--splineBetterInit', help='better spline initialization', default=False)  
    parser.add_argument('-si','--splineInit', help='spline only for the interior points', default=False)  
    parser.add_argument('-t', '--target', help='xt for the target geometry', default=5.5) 
    parser.add_argument('-m', '--method', help='Available: Nelder-Mead,Powell,CG,BFGS,L-BFGS-B,TNC,COBYLA,SLSQP,trust-constr,all', default='all') 
    parser.add_argument('-pl', '--plot', help='plot iterative process and final results', default=False) 
    parser.add_argument('-wd', '--workdir', help='working directory for q1D solver results', default='./results/')
    parser.add_argument('-tol', '--tolerance', help='stop criteria for optimization', default=1e-8) 
    parser.add_argument('-r', '--resultfile', help='result filename', default='resDefaut.pickle') 
    
    args = parser.parse_args()
    if args.method == 'all':
        methods = ['Nelder-Mead','Powell','CG','BFGS','L-BFGS-B','TNC','COBYLA','SLSQP','trust-constr']
    else:
        methods = args.method.split(',')

    global workdir 
    workdir = args.workdir

    global nQ1d
    nQ1D = 100
    global tolQ1D
    tolQ1D = 1e-8

    global tol
    tol = float(args.tolerance)

    plot = args.plot
    xtt = float(args.target)
 
    # target pressure
    p_target, noz_target = Fxt(xt=xtt)
    spt = getSpacedElements(noz_target.Sn,8)
    global st
    st = noz_target.Sn

    global interiorOnly

    if args.polynomial:    
        # initial guess
        xti = 7.5

        # polynomial parametrization optimization
        bounds=[(5.0,10.0)]
        res = {}
        resXt = opt(methods, p_target=p_target, bounds=bounds, init=xti, Fun=Fxt, plots=plot, savedict=res, filename=args.resultfile)

    if args.spline:
        interiorOnly = False 
        # use a uniform distribution for the initial guess
        spi = np.ones(8)
        bounds = [(0.0,2.5) for val in range(8)]

        bounds = [(val,val) for val in spt]
        for i,val in enumerate(bounds):
            if i not in [0,7]:
                bounds[i] = (0.0,2.5)

        resUniform = {}
        opt(methods, p_target=p_target, init=spi, Fun=Fspline, bounds=bounds, plots=plot, savedict=resUniform, filename=args.resultfile)
        
    if args.splineBetterInit: 
        interiorOnly = False   
        # use a better initial guess for spline
        xti = 7.5
        si = Area(np.linspace(0,10,nQ1D), xti)
        #spi=np.linspace(si[0],si[-1], 8)
        spi = getSpacedElements(si,8)

        bounds = [(val,val) for val in spt]
        for i,val in enumerate(bounds):
            if i not in [0,7]:
                bounds[i] = (0.0,2.5)
    
        resBetterInit = {}
        opt(methods, p_target=p_target, init=spi, Fun=Fspline, bounds=bounds, plots=plot, savedict=resBetterInit, filename=args.resultfile)

    if args.splineInit: 
        interiorOnly = True   
        # use a better initial guess for spline
        xti = 7.5
        si = Area(np.linspace(0,10,nQ1D), xti)
        #spi=np.linspace(si[0],si[-1], 8)
        spi = getSpacedElements(si,8)

        bounds = [(val-0.01,val+0.01) for val in spt]
        for i,val in enumerate(bounds):
            if i not in [0,7]:
                bounds[i] = (0.0,2.5)
    
        resBetterInit = {}
        opt(methods, p_target=p_target, init=spi[1:-1], Fun=Fspline, bounds=bounds[1:-1], plots=plot, savedict=resBetterInit, filename=args.resultfile)
