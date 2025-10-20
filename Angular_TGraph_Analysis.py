#JRP
import ROOT
import yaml
import numpy as np
import math
from collections import defaultdict
from array import array
import sys
import os
from pathlib import Path



if len(sys.argv) < 3:
    print("Proper usage is:  python3 AngularAnalysis.py startingRun# endingRun#") 
    exit(-1)
startRun = int(str(sys.argv[1]))
endRun = int(str(sys.argv[2]))


# ------------------------------- yaml file parameters ------------------------------- #

# *************** begin user input *************** #
yaml_file = 'AngularAnalysis_config.yml' 
init_params = [1, 1/8, 1/24]
#see this yaml file to edit analysis configurations
# *************** end user input *************** #

# open yaml file
with open(yaml_file,'r') as file:
    info = yaml.safe_load(file)


# parameters defined in yaml file
sample_name = info['sample_name']
force_sum_hists = info['force_sum_hists']
do_ngAnalysis = info['do_ngAnalysis']
do_ggAnalysis = info['do_ggAnalysis']
running_a_check = info['running_a_check']
stage1_hist_filepath = info['stage1_hist_filepath']
stage1_hist_prefix = info['stage1_hist_prefix']
stage1_hist_suffix = info['stage1_hist_suffix']
analysis_output_file_prefix = info['analysis_output_file_prefix'] 
#stage1_cfg_filepath = info['stage1_cfg_filepath']
summedhists_outputfile = info['summedhists_outputfile_prefix'] + f'{startRun}_{endRun}.root'
fit_limit_low = float(info['fit_limits_deg']['low'])
fit_limit_high = float(info['fit_limits_deg']['high'])
Esum_cut_low = info['Esum_cuts'][f'{sample_name}_low'] #list
Esum_cut_high = info['Esum_cuts'][f'{sample_name}_high'] #list
angle_precision = info['angle_precision']
angle_err = info['angle_err']
theta = info['detector_spherical_coords']['theta'] 
phi = info['detector_spherical_coords']['phi'] 



# ---------------------------------- define functions  ---------------------------------- #
# angular distribution W in terms of Legendre polynomials
def fitfcos_Legendre(x,params):
    theta_deg = x[0]
    cos_theta = ROOT.TMath.Cos(theta_deg * ROOT.TMath.DegToRad())  # Convert to radians

    #Legendre polynomials
    P2 = 0.5 * (3 * cos_theta**2 - 1)  # P2(cos θ)
    P4 = 0.125 * (35 * cos_theta**4 - 30 * cos_theta**2 + 3)  # P4(cos θ)

    # Angular distribution: W(θ) = 1 + a2 * P2 + a4 * P4
    a2 = params[0]
    a4 = params[1]

    return 1 + (a2 * P2) + (a4 * P4)


# angular distribution fit fcostion in terms of cos(θ)
def fitFunc_cos(x, params):
    theta_deg = x[0]  # Angle in degrees
    cos_theta = ROOT.TMath.Cos(theta_deg * ROOT.TMath.DegToRad())  # Convert to radians

    #return params[3] * (params[0] + (params[1] * (cos_theta)**2) + (params[2] * (cos_theta)**4))
    return params[0] + (params[1] * (cos_theta)**2) + (params[2] * (cos_theta)**4)



def GetAngle(det1, det2, theta, phi):
    theta_1_rad = math.radians(theta[det1])
    theta_2_rad = math.radians(theta[det2])
    delta_phi_rad = math.radians(abs(phi[det1] - phi[det2]))

    cos_angle_btwn = (
        math.cos(theta_1_rad) * math.cos(theta_2_rad)
        + math.sin(theta_1_rad) * math.sin(theta_2_rad) * math.cos(delta_phi_rad)
    )

    # Clamp for numerical safety
    cos_angle_btwn = max(-1.0, min(1.0, cos_angle_btwn))
    angle_btwn_deg = math.degrees(math.acos(cos_angle_btwn))
    return angle_btwn_deg

def compute_detector_pair_angle_counts(theta, phi, angle_precision):
    angle_counts = defaultdict(int)
    dead = {76, 86}   # 0-based IDs
    live = [True]*162
    for d in dead: live[d] = False
    for i in range(162):
        if not live[i]: continue
        for j in range(i + 1, 162): #i+1 to avoid duplicating every pair and self-pairs
            if not live[j]: continue
            angle = GetAngle(i, j, theta, phi)
            rounded_angle = round(angle, angle_precision)
            angle_counts[rounded_angle] += 1

    # Convert to sorted list of tuples like [(angle, count), ...]
    return sorted(angle_counts.items())

def compute_norm_scaling_values(theta, phi, angle_precision):
    angle_counts = compute_detector_pair_angle_counts(theta, phi, angle_precision)
    scaling_values = {}
    for i in range(len(angle_counts)):
        angle = angle_counts[i][0]
        num_pairs = angle_counts[i][1]
        scaling_values[angle] = 1 / num_pairs
    return scaling_values





'''def read_NEnResGates_Ang(stage1_cfg_filepath):
    items = []
    with open(stage1_cfg_filepath) as cfg:
        for line in cfg:
            line = line.split("#", 1)[0].strip()
            if line:
                items.extend(line.split())
    
    #find the entry for NEnResGates_Ang
    i = items.index("NEnResGates_Ang")
    n = int(items[i+1])
    values = list(map(float, items[i+2:i+2+2*n]))
    gates = [(values[j], values[j+1]) for j in range(0,len(values), 2)]
    return n, gates '''




# ---------------------------------- get filepaths  ---------------------------------- #

# initialize lists/dictionaries
# -- geometry normalization values -- #
scaling_values = compute_norm_scaling_values(theta, phi, angle_precision) 
angle_SumofWeights = defaultdict(float)
angle_SumofWeightsSqd = defaultdict(float)
counts = defaultdict(float)
weight = defaultdict(float)

filestosum = []
angles_runnum = {}
use_summed_hists = startRun != endRun #if startRun = endRun, single run, don't use summed hists

# summing stage1 histograms

if force_sum_hists:
    for runnum in range(startRun, endRun + 1):
        rootfiletosum = f"{stage1_hist_filepath}{stage1_hist_prefix}{runnum}{stage1_hist_suffix}"
        filestosum.append(rootfiletosum)
else:
    if use_summed_hists and os.path.exists(summedhists_outputfile): # if summed histogram file exists, ask user if they want to use it
        print(f"\nSummed histogram file exists")
    else: # if sum file does not exist, sum stage1 histograms 
        for runnum in range(startRun, endRun + 1):
            rootfiletosum = f"{stage1_hist_filepath}{stage1_hist_prefix}{runnum}{stage1_hist_suffix}"
            filestosum.append(rootfiletosum)

if use_summed_hists and filestosum:
    print("\nSumming histograms . . . ")
    add_hists_command = f"hadd -f {summedhists_outputfile} {' '.join(filestosum)}"
    ROOT.gSystem.Exec(add_hists_command)
        

# print histogram filepath 
if use_summed_hists:
    runRange = f'{startRun}_{endRun}'
    rootfile = ROOT.TFile(summedhists_outputfile)
    print(f"\nStarting analysis on: {summedhists_outputfile}")

else:
    runRange = startRun  
    rootfilepath = f"{stage1_hist_filepath}{stage1_hist_prefix}{startRun}{stage1_hist_suffix}"
    rootfile = ROOT.TFile(rootfilepath)
    print(f"Starting analysis on: {rootfilepath}")
  
# analysis output file
analysis_output_filename = f'{analysis_output_file_prefix}{runRange}.root'
analysis_output_file = ROOT.TFile(analysis_output_filename, "RECREATE")



if running_a_check:
    #the following check confirms that TGraphs are summing when we sum root files
    '''graphTocheck = 'T_ggAngle_Esum_M2_EnGate_0_5e+06eV'
    file = ROOT.TFile(summedhists_outputfile)
    Tgraph_sum = file.Get(graphTocheck)
    Npoints = Tgraph_sum.GetN()
    print("Summed graph has", Npoints, "points")

    for runnum in range(startRun, endRun + 1):
        singlepath = f"{stage1_hist_filepath}{stage1_hist_prefix}{runnum}{stage1_hist_suffix}"
        singlefile = ROOT.TFile(singlepath)
        Tgraph = singlefile.Get(graphTocheck)
        Npoints = Tgraph.GetN()
        print(f"single run {runnum} graph has {Npoints}" )'''
    
    #the following check is to compare txt files of angle_counts
    angle_counts = compute_detector_pair_angle_counts(theta, phi, angle_precision)
    with open("angle_counts_new.txt", "w") as f:
        for angle, count in angle_counts:
            f.write(f"{angle:.{angle_precision}f} {count}\n")
    with open("scaling_values_new.txt", "w") as f:
        for angle, scale in sorted(scaling_values.items(), key=lambda kv: kv[0]):
            f.write(f"{angle:.1f} {scale}\n")




# what type of analysis are you interested in doing? n-g, g-g, or both? see yaml file to change setup
if not do_ngAnalysis and not do_ggAnalysis:
    print(f"No analysis type chosen. See {yaml_file} to change setup.")



# ------------------------------- g-g Angular Correlation Analysis  ------------------------------- #
if do_ggAnalysis:
    # -- get graphs -- #
    gg_graphs = []
    ggAngle_normcounts = ROOT.TGraphErrors()


    #find hists with keywords
    for key in rootfile.GetListOfKeys():
        kname = key.GetName()
        
        if 'T_ggAngle' in kname:
            gg_graphs.append(kname)
            print(kname)


    print("\n")
    for gname in gg_graphs:
        print("graph name:", gname)
        graph = rootfile.Get(gname) #x = ggAngle, y = Esum, z = counts
        Npoints = graph.GetN()
        for i in range(Npoints):
            x = array('d',[0.])
            y = array('d',[0.])
            graph.GetPoint(i,x,y)
            

            x[0] = round(x[0], angle_precision)
            if Esum_cut_low <= y[0] <= Esum_cut_high:
                if x[0] in scaling_values:
                    weight[0] = scaling_values[x[0]]
                #print(x[0], y[0], weight[0])

                angle_SumofWeights[x[0]] +=weight[0] # gets normalized counts for each angle
                angle_SumofWeightsSqd[x[0]] += weight[0]**2  # for Poisson: sum( (σ=weight)^2 )

        with open('sumweights_new.txt', 'w') as f:
            for i, angle in enumerate(sorted(angle_SumofWeights.keys())):
                y = angle_SumofWeights[angle]
                ey = math.sqrt(angle_SumofWeightsSqd[angle])
                print(i, y, ey, file =f)
                ggAngle_normcounts.SetPoint(i, angle, y)
                ggAngle_normcounts.SetPointError(i, angle_err, ey)
                
            #print("point", i, "angle", angle, "scaling value", scaling_values[angle], "counts", counts[angle], "angle_SumofWeights[angle]", y)

        if gname.startswith("T_ggAngle_Esum_"):
            ending = gname[len("T_ggAngle_Esum_"):]
        else: 
            ending = gname
        
        analysis_output_file.cd()
        canvas = ROOT.TCanvas()
        ggAngle_normcounts.SetName(f"ggAngle_normcounts_{ending}")
        ggAngle_normcounts.SetTitle("Normalized g-g angular correlation;Angle (deg);Normalized Counts")
        ggAngle_normcounts.SetMarkerStyle(20)
        #print(f"number of points in {ggAngle_normcounts.GetName()}: {ggAngle_normcounts.GetN()}")
        ggAngle_normcounts.Draw("AP")
        canvas.Update()
        canvas.Write(f"ggAngle_normcounts_{ending}")

        print(f"\ng-g angular analysis output written to: {analysis_output_filename}\n") 

            
        #xmin, xmax = fit_limit_low, fit_limit_high
        #fcos = ROOT.TF1("fcos", fitFunc_cos, fit_limit_low, fit_limit_high, 4)
        fcos = ROOT.TF1("fcos", fitFunc_cos, fit_limit_low, fit_limit_high, 3)
        fcos.SetParNames("a0", "a2", "a4", "Norm")
        fcos.SetParameters(init_params[0], init_params[1], init_params[2], 1)
        fcos.SetLineColor(ROOT.kRed)
        fit = ggAngle_normcounts.Fit(fcos, "R")

        print("a0 =", fcos.GetParameter(0))
        print("a2 =", fcos.GetParameter(1))
        print("a4 =", fcos.GetParameter(2))
        print("Norm =", fcos.GetParameter(3))
    
        analysis_output_file.cd()

        fit_canvas = ROOT.TCanvas(f"fit_ggAngle_normcounts_{ending}", "Angular Correlation Fit")
        ggAngle_normcounts.SetTitle("Angular Correlation;Angle (degrees);Counts (geometry-normalized)")
        ggAngle_normcounts.SetMarkerStyle(20)
        ggAngle_normcounts.Draw("AP")
        fcos.Draw("same")
        fit_canvas.Update()
        fit_canvas.Write(f"fit_ggAngle_normcounts_{ending}")

        # -- Scaling -- #
        new_norm = fcos.GetParameter(3)
        
        for i in range(ggAngle_normcounts.GetN()):
            x = ggAngle_normcounts.GetX()[i]
            y = ggAngle_normcounts.GetY()[i]
            ex = ggAngle_normcounts.GetErrorX(i)
            ey = ggAngle_normcounts.GetErrorY(i)

            # Divide both y and error by new_norm
            ggAngle_normcounts.SetPoint(i, x, y / new_norm)
            ggAngle_normcounts.SetPointError(i, ex, ey / new_norm)

        # Refit after rescaling
        fcos = ROOT.TF1("fcos", fitFunc_cos, fit_limit_low, fit_limit_high, 4)
        fcos.SetParameters(init_params[0],init_params[1], init_params[2])  
        ggAngle_normcounts.Fit(fcos, "R")
        
        

        # Save new canvas
        analysis_output_file.cd()
        canvas = ROOT.TCanvas(f"fit_scaled_ggAngle_normcounts_{ending}", "Scaled Angular Correlation Fit")
        ggAngle_normcounts.SetTitle("Angular Correlation;Angle (degrees);Counts (geometry-normalized and scaled)")
        ggAngle_normcounts.SetMarkerStyle(20)
        ggAngle_normcounts.Draw("AP")
        fcos.Draw("same")
        canvas.Update()
        canvas.Write(f"fit_scaled_ggAngle_normcounts_{ending}")

        # Optional: Print final parameters
        a0_final = fcos.GetParameter(0)
        a2_final = fcos.GetParameter(1)
        a4_final = fcos.GetParameter(2)

        print(f"Final (normalized) fit parameters:")
        print(f"A     = {fcos.GetParameter(3):.3f}")
        print(f"a0 = {a0_final:.3f}")
        print(f"a2 = {a2_final:.3f}")
        print(f"a4 = {a4_final:.3f}")

        chi2 = fcos.GetChisquare()
        ndf = fcos.GetNDF()
        chi2_ndf = chi2 / ndf if ndf != 0 else float('nan')

        print(f"Chi2 = {chi2:.3f}")
        print(f"NDF  = {ndf}")
        print(f"Chi2/NDF = {chi2_ndf:.3f}")

        print(f"Fit complete. Results saved to {analysis_output_file}.")