#JRP
import ROOT
import numpy as np
import math
from array import array
import sys
from collections import defaultdict




# --------------------- user config  --------------------- #

# -------------------------------------------------------- #

# avg estimated angle error
# (about 16 degrees between nearest neighbors: 
# 16/2 = approximately 8 degrees to the edge of the crystal)
angle_err = 8.0

if len(sys.argv) < 4:
    print("Proper usage is:  python3 GammaAngularAnalysis.py <analysis type (ng/gg)> <starting run#> <ending run#>") 
    exit(-1)
analysis_type = str(sys.argv[1])
startRun = int(str(sys.argv[2]))
endRun = int(str(sys.argv[3]))


# -------------------------------- functions  -------------------------------- # 
def ang_distrib_func_Legendre(x,params):
        theta_deg = x[0]
        cos_theta = ROOT.TMath.Cos(theta_deg * ROOT.TMath.DegToRad())  # Convert to radians

        #Legendre polynomials
        P0 = 1
        P1 = cos_theta
        P2 = (1/2) * (3 * cos_theta**2 - 1)  # P2(cos θ)
        P3 = (1/2) * (5* cos_theta**3 - 3*cos_theta)
        P4 = (1/8) * (35 * cos_theta**4 - 30 * cos_theta**2 + 3)  # P4(cos θ)

        # Angular distribution function W(θ) = ∑_k (a_k * P_k)
        a0 = params[0]
        a1 = params[1]
        a2 = params[2]
        a3 = params[3]
        a4 = params[4]

        return (a0*P0) +(a1* P1) + (a2 * P2) + (a3 * P3) + (a4 * P4)

def ang_distrib_func(x,params):
        theta_deg = x[0]
        cos_theta = ROOT.TMath.Cos(theta_deg * ROOT.TMath.DegToRad())  # Convert to radians

        #Legendre polynomials
        C0 = 1
        C1 = cos_theta
        C2 = cos_theta**2
        C3 = cos_theta**3
        C4 = cos_theta**4

        # Angular distribution function W(θ) = ∑_k (a_k * P_k)
        a0 = params[0]
        a1 = params[1]
        a2 = params[2]
        a3 = params[3]
        a4 = params[4]

        return (a0*C0) +(a1* C1) + (a2 * C2) + (a3 * C3) + (a4 * C4)


 # -------------------------------------  (n,γ) ANGULAR DISTRIBUTION ANALYSIS ------------------------------------- #
if analysis_type == 'ng':
    if startRun == endRun:
        print(f"Beginning (n,γ) angular distribution analysis on run {startRun}")
        
    else: 
        # sum hists?
        print(f"Beginning (n,γ) angular distribution analysis on runs {startRun}-{endRun}")


   
    

 # -------------------------------------  γ-γ ANGULAR CORRELATION ANALYSIS ------------------------------------- #
elif analysis_type == 'gg':
    if startRun == endRun:
        print(f"Beginning γ-γ angular correlation analysis on run {startRun}")
    else:
        # sum hists?
        print(f"Beginning γ-γ angular correlation analysis on runs {startRun}-{endRun}")


else:
    print("The analysis type must be 'ng' or 'gg'. ")
    print("Proper usage is:  python3 GammaAngularAnalysis.py <analysis type (ng/gg)> <starting run#> <ending run#>") 
    exit(-1)