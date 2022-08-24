# Author: Allan Moreira de Carvalho
# Created: 23.08.2022

import numpy as np
import os
import subprocess
import shutil
import timeit
import pandas as pd

class gmsh:
    pcount = 0
    lcount = 0
    loopcount = 0
    scount = 0
    phycount = 0
    gfile = ''

    def addPoint(self, x, y, z):
        self.pcount += 1
        self.gfile += 'Point('+str(self.pcount)+') = {'+str(x)+','+str(y)+','+str(z)+', cl__1};\n'
        return self.pcount

    def addLine(self, points):
        self.lcount += 1
        self.gfile += 'Line('+str(self.lcount)+') = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.lcount

    def addSpline(self, points):
        self.lcount += 1
        self.gfile += 'Spline('+str(self.lcount)+') = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.lcount

    def addLoop(self, lines):
        self.loopcount += 1
        self.gfile += 'Line Loop('+str(self.loopcount)+') = {'
        for i,l in enumerate(lines):
            self.gfile += str(l)
            if i < np.size(lines)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.loopcount

    def addPlaneSurface(self, loops):
        self.scount += 1
        self.gfile += 'Plane Surface('+str(self.scount)+') = {'
        for i,l in enumerate(loops):
            self.gfile += str(l)
            if i < np.size(loops)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.scount
            
    def transfiniteLine(self, line, points, progression=1):
        self.gfile += 'Transfinite Line {' + str(line) + '} = ' + str(points) + ' Using Progression ' + str(progression) + ';\n'

    def transfiniteSurface(self, surface, points):
        self.gfile += 'Transfinite Surface {' + str(surface) + '} = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
    def recombineSurface(self, surface):
        self.gfile += 'Recombine Surface {' + str(surface) + '};\n'

    def addPhysicalLine(self, boundName, line):
        self.phycount +=1
        self.gfile += 'Physical Line(\"' + boundName +'\") = {' + str(line) + '};\n'

    def addPhysicalSurface(self, surface):
        self.phycount +=1
        self.gfile += 'Physical Surface(' + str(self.phycount) +') = {' + str(surface) + '};\n'
        
    def writeLine(self, s):
        self.gfile += s

    def print(self):
        print(self.gfile)

    def writeFile(self, filename):
        f = open(filename,'w')
        f.write(self.gfile)
        f.close()

class nozzle:
    def setX(self, x):
        self.xn = x

    def setXSU2(self, xi, xf, Nn):
        self.Nx = Nn
        self.xnSU2 = np.linspace(xi, xf, Nn)

    def setNySU2(self, Ny):
        self.Ny = Ny

    def setS(self, s):
        self.Sn = s

    def setSSU2(self, f, *args):
        self.SnSU2 = f(self.xnSU2, *args)

    def setBC(self, p0in, T0in, Min, pb):
        self.p0in = p0in
        self.T0in = T0in
        self.pb = pb
        self.Min = Min

    def setFluid(self, R, gamma):
        self.R = R
        self.gamma = gamma

    def setQ1DSolver(self, itamx, itprint, CFL, tol,
                  tscheme, fscheme, dttype, dim):
        self.itmax = itamx
        self.itprint = itprint
        self.CFL = CFL
        self.tol = tol
        self.tscheme = tscheme
        self.fscheme = fscheme
        self.dttype = dttype
        self.dim = dim
    
    def setSolverQ1D(self, solver):
        self.solver = solver

    def setupQ1D(self, filepath, filename):
        self.filepath = filepath
        self.filename = filename
        infilepath = filepath + 'inputs/'
        self.q1dinfilepath = infilepath

        self.setupfile = infilepath + filename
        os.makedirs(infilepath, exist_ok=True)

        outfilepath = filepath + 'outputs/'
        self.q1doutfilepath = outfilepath
        
        os.makedirs(outfilepath, exist_ok=True)
        
        np.savetxt(infilepath+'xn.txt', self.xn)
        np.savetxt(infilepath+'Sn.txt', self.Sn)

        f = open(infilepath + filename,'w')
        f.write('# Domain x-coordinates at cell faces (no need for ghost cells)\n')
        f.write(infilepath + 'xn.txt' + "\n")
        f.write('# Area distributuin at cell faces (no need for ghost cells)\n')
        f.write(infilepath + 'Sn.txt' + "\n")
        f.write('# Inlet Total Pressure [Pa]\n')
        f.write(str(self.p0in) + "\n")
        f.write('# Inlet Total Temperature [K]\n')
        f.write(str(self.T0in) + "\n")
        f.write('# Inlet Mach Number\n')
        f.write(str(self.Min) + "\n")
        f.write('# Outlet Static Pressure [Pa]\n')
        f.write(str(self.pb) + "\n")
        f.write('# Gas constant [J/kg/K]\n')
        f.write(str(self.R) + "\n")
        f.write('# Specific heat ratio\n')
        f.write(str(self.gamma) + "\n")
        f.write('# Maximum number of iterations \n')
        f.write(str(self.itmax) + "\n")
        f.write('# Interval to print iterations \n')
        f.write(str(self.itprint) + "\n")
        f.write('# CFL number \n')
        f.write(str(self.CFL) + "\n")
        f.write('# Convergence criteria \n')
        f.write(str(self.tol) + "\n")
        f.write('# Time integration scheme (Euler, RK4) \n')
        f.write(str(self.tscheme) + "\n")
        f.write('# Flux scheme (Basic, Roe, Split, SplitTVD, VanLeer, LaxFriedrichs, StegerWarming, AUSM) \n')
        f.write(str(self.fscheme) + "\n")
        f.write('# Timestep type (Global, Local) \n')
        f.write(str(self.dttype) + "\n")
        f.write('# Form to be solved (Dimensional, Dimensionless) \n')
        f.write(str(self.dim) + "\n")
        f.write('# Set output directory \n')
        f.write( outfilepath + "\n")
        f.close()
        
    def getQ1D(self, propfile):
        return np.loadtxt(self.filepath + propfile)

    def solveQ1D(self):
        self.setupQ1D(self.filepath, self.filename)
        self.runtimeQ1D = timeit.default_timer()
        print([self.solver, self.setupfile])
        p = subprocess.Popen([self.solver, self.setupfile])
        p.wait()
        self.runtimeQ1D = timeit.default_timer() - self.runtimeQ1D
        print('runtime: ' + str(self.runtimeQ1D) + 's')
        #input()
    
    def readQ1DOutput(self, prop):
        try:
            return np.loadtxt(self.q1doutfilepath+prop)
        except:
            print('file {} cant be reached')

    @property
    def c(self):
        return self.readQ1DOutput('c.txt')

    @property
    def e(self):
        return self.readQ1DOutput('e.txt')

    @property
    def M(self):
        return self.readQ1DOutput('M.txt')

    @property
    def p(self):
        return self.readQ1DOutput('p.txt')

    @property
    def rese(self):
        return self.readQ1DOutput('rese.txt')
    
    @property
    def resrho(self):
        return self.readQ1DOutput('resrho.txt')

    @property
    def resrhou(self):
        return self.readQ1DOutput('resrhou.txt')

    @property
    def rho(self):
        return self.readQ1DOutput('rho.txt')

    @property
    def S(self):
        return self.readQ1DOutput('S.txt')

    @property
    def T(self):
        return self.readQ1DOutput('T.txt')
    
    @property
    def u(self):
        return self.readQ1DOutput('u.txt')

    @property
    def x(self):
        return self.readQ1DOutput('x.txt')

    def genGmshNozzle(self, outputpath, fname):
        L = self.xnSU2[-1]-self.xnSU2[0]
        x = self.xnSU2
        r = self.SnSU2/2
        Nx = self.Nx
        Ny = self.Ny
        inflationRate = 1.03

        geo = gmsh()

        geo.writeLine('cl__1 = 1;\n')
        p1 = geo.addPoint(0, 0, 0)
        p2 = geo.addPoint(L, 0, 0)
        # write wall countour
        wp = []
        for xi,ri in zip(x,r):
            wp.append(  geo.addPoint(xi, ri, 0) )

        l1 = geo.addLine([p1,p2])
        geo.transfiniteLine(l1, Nx, 1)

        l2 = geo.addLine([wp[0],p1])
        geo.transfiniteLine(l2, Ny, inflationRate)

        l3 = geo.addLine([wp[-1],p2])
        geo.transfiniteLine(l3, Ny, inflationRate)

        l4 = geo.addSpline(wp)
        geo.transfiniteLine(l4, Nx, 1)
                        
        loop1 = geo.addLoop([l4,l3,-l1,-l2])
        surface1 = geo.addPlaneSurface([loop1])

        geo.transfiniteSurface(surface1, [wp[0],wp[-1],p2,p1])
        geo.recombineSurface(surface1)

        # add boundary
        geo.addPhysicalLine("WALL",l4)
        geo.addPhysicalLine("INFLOW",l2)
        geo.addPhysicalLine("OUTFLOW",l3)
        geo.addPhysicalLine("SYMMETRY",l1)

        geo.addPhysicalSurface(surface1)

        self.geo = geo

        os.makedirs(outputpath, exist_ok=True)

        fname = fname.split('.cfg')[0]
        geofilename = outputpath + fname + '.geo'

        geo.writeFile(geofilename)
        
        p = subprocess.Popen(['./gmsh', geofilename,'-0','-2','-format','su2','-o', outputpath +fname+'.su2'])
        p.wait()
        p = subprocess.Popen(['./gmsh', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.msh'])
        p.wait()
        p = subprocess.Popen(['./gmsh', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.vtk'])
        p.wait()

        self.meshfilename = outputpath +fname+'.su2'

    def setupSU2Solver(self, itmax, itprint, CFL, tol, tscheme, fscheme, dim):
        self.itmaxSU2 = itmax
        self.itprintSU2 = itprint
        self.CFLSU2 = CFL
        self.tolSU2 = tol
        self.tschemeSU2 = tscheme
        self.fschemeSU2 = fscheme
        self.dimSU2 = dim

    def solveSU2(self, solver):
        self.runtimeSU2 = timeit.default_timer()
        
        p = subprocess.Popen([solver, self.su2cfgfile])
        p.wait()
        # move results to output file
        os.makedirs(self.su2outfilepath, exist_ok=True)
        try:
            shutil.move("./history.csv", self.su2outfilepath + "history.csv")
            shutil.move("./solution.csv", self.su2outfilepath +  "solution.csv")
            shutil.move("./solution.vtk", self.su2outfilepath +  "solution.vtk")
        except:
            print("No solution files exists")

        self.runtimeSU2 = timeit.default_timer() - self.runtimeSU2
        print('runtime: ' + str(self.runtimeSU2) + 's')
        #input() 

    def setupSU2(self, su2filepath, su2filename):
        self.su2infilepath = su2filepath
        self.su2filename = su2filename
        self.su2infilepath = su2filepath + 'inputs/'
        self.su2outfilepath = su2filepath + 'outputs/'
        self.su2cfgfile = self.su2infilepath + su2filename

        # Generate the mesh files
        self.genGmshNozzle(self.su2infilepath, su2filename)

        # Setup the SU2 .cfg file
        meshfile = self.meshfilename
        solfile = 'solution.dat'

        newcfg = f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% SU2 configuration file                                                       %
% Case description: Non-ideal compressible fluid flow in a converging-         %
%                   diverging supersonic nozzle for siloxane fluid MDM         %
% Author: Alberto Guardone                                                     %
% Institution: Politecnico di Milano                                           %
% Date: 2019.05.03                                                             %
% File Version 6.2.0 "Falcon"                                                  %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%
%
% Physical governing equations (EULER, NAVIER_STOKES,
%                               FEM_EULER, FEM_NAVIER_STOKES, FEM_RANS, FEM_LES,
%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,
%                               POISSON_EQUATION)
SOLVER= RANS
%
% Specify turbulence model (NONE, SA, SA_NEG, SST, SA_E, SA_COMP, SA_E_COMP)
KIND_TURB_MODEL= SST
%
% Mathematical problem (DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT)
MATH_PROBLEM= DIRECT
%
% Restart solution (NO, YES)
RESTART_SOL= NO
%
% System of measurements (SI, US)
% International system of units (SI): ( meters, kilograms, Kelvins,
%                                       Newtons = kg m/s^2, Pascals = N/m^2,
%                                       Density = kg/m^3, Speed = m/s,
%                                       Equiv. Area = m^2 )
% United States customary units (US): ( inches, slug, Rankines, lbf = slug ft/s^2,
%                                       psf = lbf/ft^2, Density = slug/ft^3,
%                                       Speed = ft/s, Equiv. Area = ft^2 )
SYSTEM_MEASUREMENTS= SI
%
% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%
%
% Mach number (non-dimensional, based on the free-stream values)
MACH_NUMBER= 0.01
%
% Angle of attack (degrees, only for compressible flows)
AOA= 0.0
%
% Side-slip angle (degrees, only for compressible flows)
SIDESLIP_ANGLE= 0.0
%
% Init option to choose between Reynolds (default) or thermodynamics quantities
% for initializing the solution (REYNOLDS, TD_CONDITIONS)
INIT_OPTION= TD_CONDITIONS
%
% Free-stream option to choose between density and temperature (default) for
% initializing the solution (TEMPERATURE_FS, DENSITY_FS)
FREESTREAM_OPTION= TEMPERATURE_FS
%
% Free-stream pressure (101325.0 N/m^2, 2116.216 psf by default)
FREESTREAM_PRESSURE = 104074.0

% Free-stream temperature (288.15 K, 518.67 R by default)
FREESTREAM_TEMPERATURE= 291.3
%
% Compressible flow non-dimensionalization (DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,
%                              FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)
REF_DIMENSIONALIZATION= {self.dimSU2}

% ---- IDEAL GAS, POLYTROPIC, VAN DER WAALS AND PENG ROBINSON CONSTANTS -------%
%
% Fluid model (STANDARD_AIR, IDEAL_GAS, VW_GAS, PR_GAS,
%              CONSTANT_DENSITY, INC_IDEAL_GAS, INC_IDEAL_GAS_POLY)
FLUID_MODEL= IDEAL_GAS
%
% Ratio of specific heats (1.4 default and the value is hardcoded
%                          for the model STANDARD_AIR, compressible only)
GAMMA_VALUE= {self.gamma}
%
% Specific gas constant (287.058 J/kg*K default and this value is hardcoded
%                        for the model STANDARD_AIR, compressible only)
GAS_CONSTANT= {self.R}
%
% Critical Temperature (131.00 K by default)
CRITICAL_TEMPERATURE= 131.00
%
% Critical Pressure (3588550.0 N/m^2 by default)
CRITICAL_PRESSURE= 3588550.0
%
% Acentric factor (0.035 (air))
ACENTRIC_FACTOR= 0.035

% --------------------------- VISCOSITY MODEL ---------------------------------%
%
% Viscosity model (SUTHERLAND, CONSTANT_VISCOSITY).
VISCOSITY_MODEL= SUTHERLAND
%
% Sutherland Viscosity Ref (1.716E-5 default value for AIR SI)
MU_REF= 1.716E-5
%
% Sutherland Temperature Ref (273.15 K default value for AIR SI)
MU_T_REF= 273.15
%
% Sutherland constant (110.4 default value for AIR SI)
SUTHERLAND_CONSTANT= 110.4

% --------------------------- THERMAL CONDUCTIVITY MODEL ----------------------%
%
% Conductivity model (CONSTANT_CONDUCTIVITY, CONSTANT_PRANDTL).
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
%
% Laminar Prandtl number (0.72 (air), only for CONSTANT_PRANDTL)
PRANDTL_LAM= 0.72
%
% Turbulent Prandtl number (0.9 (air), only for CONSTANT_PRANDTL)
PRANDTL_TURB= 0.90

% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%
%
% Navier-Stokes (no-slip), constant heat flux wall  marker(s) (NONE = no marker)
% Format: ( marker name, constant heat flux (J/m^2), ... )
MARKER_HEATFLUX= ( WALL, 0.0 )
%
% Symmetry boundary marker(s) (NONE = no marker)
MARKER_SYM= ( SYMMETRY )
%
% Riemann boundary marker(s) (NONE = no marker)
% Format: (marker, data kind flag, list of data)
MARKER_RIEMANN= ( INFLOW, TOTAL_CONDITIONS_PT, {self.p0in}, {self.T0in}, 1.0, 0.0, 0.0, OUTFLOW, STATIC_PRESSURE, {self.pb}, 0.0, 0.0, 0.0, 0.0 )

% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%
%
% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)
NUM_METHOD_GRAD= GREEN_GAUSS
%
% CFL number (initial value for the adaptive CFL number)
CFL_NUMBER= {self.CFLSU2}
%
% Adaptive CFL number (NO, YES)
CFL_ADAPT= YES
%
% Parameters of the adaptive CFL number (factor down, factor up, CFL min value,
%                                        CFL max value )
CFL_ADAPT_PARAM= ( 0.1, 2.0, {self.CFLSU2/100} , {self.CFLSU2*100} )
%
% Maximum Delta Time in local time stepping simulations
MAX_DELTA_TIME= 1E6

% ----------- SLOPE LIMITER AND DISSIPATION SENSOR DEFINITION -----------------%
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the flow equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_FLOW= YES
%
% Slope limiter (NONE, VENKATAKRISHNAN, VENKATAKRISHNAN_WANG,
%                BARTH_JESPERSEN, VAN_ALBADA_EDGE)
SLOPE_LIMITER_FLOW= NONE
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the turbulence equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_TURB= NO

% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%
%
% Linear solver or smoother for implicit formulations (BCGSTAB, FGMRES, SMOOTHER_JACOBI,
%                                                      SMOOTHER_ILU, SMOOTHER_LUSGS,
%                                                      SMOOTHER_LINELET)
LINEAR_SOLVER= FGMRES
%
% Preconditioner of the Krylov linear solver (ILU, LU_SGS, LINELET, JACOBI)
LINEAR_SOLVER_PREC= ILU
%
% Linael solver ILU preconditioner fill-in level (0 by default)
LINEAR_SOLVER_ILU_FILL_IN= 0
%
% Minimum error of the linear solver for implicit formulations
LINEAR_SOLVER_ERROR= 1E-6
%
% Max number of iterations of the linear solver for the implicit formulation
LINEAR_SOLVER_ITER= 10

% -------------------------- MULTIGRID PARAMETERS -----------------------------%
%
% Multi-grid levels (0 = no multi-grid)
MGLEVEL= 0

% -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------%
%
% Convective numerical method (JST, LAX-FRIEDRICH, CUSP, ROE, AUSM, AUSMPLUSUP, AUSMPLUSUP2, HLLC,
%                              TURKEL_PREC, MSW, FDS)
CONV_NUM_METHOD_FLOW= {self.fschemeSU2}
%
% Entropy fix coefficient (0.0 implies no entropy fixing, 1.0 implies scalar
%                          artificial dissipation)
ENTROPY_FIX_COEFF= 0.1
%
% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT, EULER_EXPLICIT)
TIME_DISCRE_FLOW= {self.tschemeSU2}

% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%
%
% Convective numerical method (SCALAR_UPWIND)
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
%
% Time discretization (EULER_IMPLICIT)
TIME_DISCRE_TURB= EULER_IMPLICIT
%
% Reduction factor of the CFL coefficient in the turbulence problem
CFL_REDUCTION_TURB= 1.0

% --------------------------- CONVERGENCE PARAMETERS --------------------------%
%
% Number of total iterations
ITER= {self.itmaxSU2}
%
% Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= {self.tolSU2}
%
% Start convergence criteria at iteration number
CONV_STARTITER= 10

% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%
%
% Mesh input file
MESH_FILENAME= {self.meshfilename}
%
% Mesh input file format (SU2, CGNS)
MESH_FORMAT= SU2
%
% Mesh output file
MESH_OUT_FILENAME= mesh_out.su2
%
% Restart flow input file
SOLUTION_FILENAME= solution_flow.dat
%
% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW, PARAVIEW_BINARY,
%                     FIELDVIEW, FIELDVIEW_BINARY)
OUTPUT_FILES = (CSV, PARAVIEW_ASCII)

TABULAR_FORMAT= CSV
%
% Output file convergence history (w/o extension)
CONV_FILENAME= history
%
% Output file restart flow
RESTART_FILENAME= solution.csv
%
% Output file flow (w/o extension) variables
VOLUME_FILENAME= solution
%
% Output file surface flow coefficient (w/o extension)
SURFACE_FILENAME= surface_flow
%
% Writing solution file frequency
OUTPUT_WRT_FREQ= {self.itmaxSU2}
%
% Screen output
SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, RMS_TKE, RMS_DISSIPATION, LIFT, DRAG)
% Writing frequency for screen output
SCREEN_WRT_FREQ_INNER= {self.itprintSU2}
%
SCREEN_WRT_FREQ_OUTER= {self.itprintSU2}
%
SCREEN_WRT_FREQ_TIME= {self.itprintSU2}
"""
        cfgfile = open(self.su2cfgfile, 'w')
        cfgfile.write(newcfg)
        cfgfile.close()