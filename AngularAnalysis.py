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
import csv
from colorama import Fore, Style, init
import re
from typing import Any, Iterable, Tuple, Dict

'''DANCE_DATA = os.getenv("DANCE_DATA", "/DANCE/DANCE_Analysis")
STAGE0 = os.path.join(DANCE_DATA, "stage0_bin")
STAGE1 = os.path.join(DANCE_DATA, "stage1_root")'''


# Validate input
if len(sys.argv) < 3 or len(sys.argv[1:]) % 2 != 0:
    print("Usage: python3 AngularAnalysis.py <startRun1> <endRun1> [<startRun2> <endRun2> ...]")
    sys.exit(1)

# Parse run ranges in pairs
run_ranges = []
args = sys.argv[1:]
for i in range(0, len(args), 2):
    start = int(args[i])
    end = int(args[i+1])
    run_ranges.append((start, end))

# Determine if we should use summed histograms
if len(run_ranges) > 1:
    run_range_parts = [f"{start}_{end}" for start, end in run_ranges]
elif run_ranges[0][0] != run_ranges[0][1]:
    run_range_parts = [f"{start}_{end}"]
else:
    run_range_parts = [f"{start}"]

run_range_label = "_".join(run_range_parts)





# *************** begin user input *************** #
yaml_file = 'AngularAnalysis_config.yml' 
#see this yaml file to edit analysis configurations
# *************** end user input *************** #


# open yaml file
with open(yaml_file,'r') as file:
    info = yaml.safe_load(file)





# =============================================== FUNCTIONS  ===============================================#
# angular distribution W in terms of Legendre polynomials
def fitFunc_Legendre(x,params):
    theta_deg = x[0]
    cos_theta = ROOT.TMath.Cos(theta_deg * ROOT.TMath.DegToRad())  # Convert to radians

    #Legendre polynomials
    P2 = 0.5 * (3 * cos_theta**2 - 1)  # P2(cos θ)
    P4 = 0.125 * (35 * cos_theta**4 - 30 * cos_theta**2 + 3)  # P4(cos θ)

    # Angular distribution: W(θ) = 1 + a2 * P2 + a4 * P4
    a2 = params[0]
    a4 = params[1]

    return 1 + (a2 * P2) + (a4 * P4)




# angular distribution fit function in terms of angle
def fitFunc_ang(x, params):
    theta_deg = x[0]  # Angle in degrees
    cos_theta = ROOT.TMath.Cos(theta_deg * ROOT.TMath.DegToRad())  # Convert to radians

    #return params[3] * (params[0] + (params[1] * (cos_theta)**2) + (params[2] * (cos_theta)**4))
    return params[0] + (params[1] * (cos_theta)**2) + (params[2] * (cos_theta)**4)
    # a_0 + a_2 * cos^2(theta) + a_4 * cos^4(theta)




# angular distribution fit function in terms of cos(angle)
def fitFunc_cos(x, params):
    cos_theta = x[0]
    return params[0] + (params[1] * (cos_theta)**2) + (params[2] * (cos_theta)**4)




# angular distribution fit function in terms of cos^2(angle)
def fitFunc_cos_sq(x, params):
    cos_sq = x[0]
    return params[0] + (params[1] * cos_sq) + (params[2] * (cos_sq)**2)




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





def compute_detector_pair_angle_counts(theta, phi, gg, ng):
    angle_counts = defaultdict(int)
    dead = {76,86}
    n = len(theta)
    live = [True]*n
    for d in dead: live[d] = False

    if gg == True:
        for i in range(n):
            if not live[i]: continue
            for j in range(i + 1, n): #i+1 to avoid duplicating every pair and self-pairs
                if not live[j]: continue
                angle = GetAngle(i, j, theta, phi)
                angle_int = int(angle) #make it an integer
                #print(f"{angle_int} deg between {i} and {j}")
                angle_counts[angle_int] += 1
                #print(f"angle={angle}   rounded_angle={rounded_angle}    angle_counts={angle_counts[rounded_angle]}")

    
    if ng == True:
        for i in range(n):
            if not live[i]: continue 
            angle = GetAngle(76, i, theta, phi)  # or use 180 - theta[i]
            #rounded_angle = round(angle / angleBin_width) * angleBin_width
            angle_int = int(angle)
            angle_counts[angle_int] += 1

    # Convert to sorted list of tuples like [(angle, count), ...]
    return sorted(angle_counts.items())





def compute_norm_scaling_values(theta, phi, angleBin_width, gg, ng):
    angle_counts = compute_detector_pair_angle_counts(theta, phi, gg, ng)
    NAngleBins = math.ceil(180 / angleBin_width)
    #for i in range(1, NAngleBins+1):
    #    angleBin_center = (i - 0.5) * angleBin_width
    #    print(f"angle bin {i} center: {angleBin_center}") #the index
    
    # Sum counts over each angle bin, then invert to get one scale per bin.
    # Keys are bin centers (deg), rounded to 6 dp.
    counts_by_angle = dict(angle_counts)  # {int_angle: count}

    NAngleBins = int(math.ceil(180.0 / angleBin_width))

    scaling_values = {}

    for b in range(1, NAngleBins + 1):
        L = (b - 1) * angleBin_width         # left edge
        U = b * angleBin_width               # right edge
        U_eff = 180.0 if b == NAngleBins else U  # last bin includes 180

        a_lo = int(math.ceil(L - 1e-12))         # smallest integer angle in bin
        a_hi = int(math.floor(U_eff - 1e-12))    # largest integer angle in bin
        #print(f"a_lo: {a_lo}     a_hi: {a_hi}")
        total_pairs = 0
        for a in range(a_lo, a_hi + 1):
            total_pairs += counts_by_angle.get(a, 0)

        if total_pairs > 0:
            center = round(0.5 * (L + U_eff), 1)     # bin center (deg)
            #print(f"center: {center}")
            scaling_values[center] = 1.0 / total_pairs

    return scaling_values
    '''
    scaling_values = {}
    for i in range(len(angle_counts)):
        angle = angle_counts[i][0]
        num_pairs = angle_counts[i][1]
        if num_pairs > 0.0:
            scaling_values[angle] = 1 / num_pairs
    return scaling_values'''




'''def geom_scaling_bbb(theta, phi, xaxis, gg, ng):
    """
    Return a dict: binx -> 1/num_pairs for that bin (same bins as hist x-axis).
    """

    nbins = xaxis.GetNbins()
    counts_by_binx = [0] * (nbins + 1)  # 1..nbins (0 unused)
    dead = {76, 86}
    n = len(theta)
    live = [i not in dead for i in range(n)]

    if gg:
        for i in range(n):
            if not live[i]: 
                continue
            for j in range(i + 1, n):
                if not live[j]:
                    continue
                ang = int(GetAngle(i, j, theta, phi))  # in degrees
                binx = xaxis.FindBin(ang)         # ROOT decides the bin
                if 1 <= binx <= nbins:
                    counts_by_binx[binx] += 1

    if ng:
        for i in range(n):
            if not live[i]:
                continue
            ang = GetAngle(76, i, theta, phi)
            binx = xaxis.FindBin(ang)
            if 1 <= binx <= nbins:
                counts_by_binx[binx] += 1

    scaling_by_binx = {}
    for b in range(1, nbins + 1):
        c = counts_by_binx[b]
        if c > 0:
            scaling_by_binx[b] = 1.0 / c
    return scaling_by_binx




def aggregate_scaling_to_hist_bins(scaling_values_by_int_angle, xaxis):
    """
    scaling_values_by_int_angle: dict { int_angle : 1/count_at_that_angle }
    xaxis: ROOT TAxis of the histogram x-axis

    Returns:
      scale_by_binx:   dict { binx (1..nbins) : 1 / sum_a (1/scale[a]) }
      scale_by_center: dict { rounded_bin_center : same value }  # convenience
    """
    nbins = xaxis.GetNbins()
    scale_by_binx = {}
    scale_by_center = {}

    for b in range(1, nbins + 1):
        L = xaxis.GetBinLowEdge(b)
        U = xaxis.GetBinUpEdge(b)

        # Half-open [L, U) except make the LAST bin include U (ROOT convention)
        if b < nbins:
            U_eff = U - 1e-9
        else:
            U_eff = U

        # Integer angles produced by your function use int(angle) (truncate → floor for positives).
        # So include integer 'a' that could have come from any angle in this bin:
        a_lo = math.ceil(L)               # smallest integer >= L
        a_hi = math.floor(U_eff)          # largest integer <= U_eff

        # Sum counts in this bin via your scales (count = 1/scale)
        N_b = 0.0
        for a in range(a_lo, a_hi + 1):
            s = scaling_values_by_int_angle.get(a)
            if s:                         # only if that integer angle existed
                N_b += 1.0 / s

        if N_b > 0.0:
            sf_b = 1.0 / N_b
            scale_by_binx[b] = sf_b

            c = round(xaxis.GetBinCenter(b), 6)  # stable float key
            scale_by_center[c] = sf_b
        # else: leave missing → no geometric acceptance in that bin

    return scale_by_binx, scale_by_center


def _bin_edges_deg(xaxis, binx):
    L = xaxis.GetBinLowEdge(binx)
    U = xaxis.GetBinUpEdge(binx)
    return L, U

def bin_avg_theta_deg(xaxis, binx):
    """⟨θ⟩ over [L,U] assuming uniform weight in θ (degrees)."""
    L, U = _bin_edges_deg(xaxis, binx)
    return 0.5 * (L + U)

def bin_avg_cos(xaxis, binx):
    """⟨cos θ⟩ over [L,U] (θ in degrees on axis; integral done in radians)."""
    L, U = _bin_edges_deg(xaxis, binx)
    Lr, Ur = math.radians(L), math.radians(U)
    if Ur == Lr:
        return math.cos(Lr)
    return (math.sin(Ur) - math.sin(Lr)) / (Ur - Lr)

def bin_avg_cos2(xaxis, binx):
    """⟨cos² θ⟩ over [L,U]."""
    L, U = _bin_edges_deg(xaxis, binx)
    Lr, Ur = math.radians(L), math.radians(U)
    if Ur == Lr:
        c = math.cos(Lr)
        return c*c
    #  ⟨cos²θ⟩ = 1/2 + [sin(2U) - sin(2L)] / [4 (U-L)]
    return 0.5 + (math.sin(2*Ur) - math.sin(2*Lr)) / (4*(Ur - Lr))

def bin_xvalue(xaxis, binx, mode="angle"):
    """
    mode: 'angle' -> ⟨θ⟩ (deg), 'cos' -> ⟨cos θ⟩, 'cos2' -> ⟨cos² θ⟩
    Use this to set the x-value for TGraphErrors.
    """
    if mode == "angle":
        return bin_avg_theta_deg(xaxis, binx)
    elif mode == "cos":
        return bin_avg_cos(xaxis, binx)
    elif mode == "cos2":
        return bin_avg_cos2(xaxis, binx)
    raise ValueError("mode must be 'angle', 'cos', or 'cos2'")'''





'''def correct_1deg_then_sum_to_bins(hist2d, Esum_low, Esum_high,
                                  theta, phi, gg, ng,
                                  coarse_width):
    """
    Step 1: Project 2D (θ × Esum) to 1° θ after the Esum gate.
    Step 2: Divide per-degree by geometric counts G_a (your function).
    Step 3: Sum corrected 1° bins into coarse bins of width 'coarse_width'.
    Returns: x (bin centers in deg), y (corrected counts), ey (errors).
    """
    xaxis = hist2d.GetXaxis()
    yaxis = hist2d.GetYaxis()

    # Project to X with error propagation
    ylo = yaxis.FindBin(Esum_low  + 1e-9)
    yhi = yaxis.FindBin(Esum_high - 1e-9)
    #hX_fine = hist2d.ProjectionX(f"{hist2d.GetName()}_px1deg", ylo, yhi, "e")
    hX_fine = hist2d.ProjectionX(f"{hist2d.GetName()}_px1deg", ylo, yhi)
    # Build geometric G_a at 1°
    angle_counts = compute_detector_pair_angle_counts(theta, phi, gg, ng)

    G = np.zeros(181, dtype=float)   # indices 0..180 map to angles 0..180
    for a, c in angle_counts:
        if 0 <= a <= 180:
            G[a] = float(c)
    #print(G[90]) worked :: 202
    
    # Acceptance-correct per degree
    S = np.zeros(181, dtype=float)
    E2 = np.zeros(181, dtype=float)
    for a in range(181):
        m  = hX_fine.GetBinContent(a + 1)
        em = hX_fine.GetBinError(a + 1)
        if G[a] > 0:
            S[a]  = m / G[a]
            E2[a] = (em / G[a])**2
        else:
            # No geometric acceptance at this angle; skip it
            S[a]  = 0.0
            E2[a] = 0.0

    # Coarse edges from 0..180 (last bin includes 180)
    nb = int(math.floor(180.0 / coarse_width))
    edges = [i * coarse_width for i in range(nb + 1)]
    if edges[-1] < 180.0 - 1e-12:
        edges.append(180.0)

    # Sum corrected 1° bins into coarse bins
    x, y, ey = [], [], []
    for b in range(1, len(edges)):
        L, U = edges[b-1], edges[b]
        U_eff = 180.0 if b == len(edges) - 0 else U  # include 180 in last bin
        a_lo = max(0,   int(math.ceil(L  - 1e-12)))
        a_hi = min(180, int(math.floor(U_eff - 1e-12)))
        #print(f"a_lo: {a_lo}  a_hi: {a_hi}")
        y_b  = S[a_lo:a_hi+1].sum()
        ey_b = math.sqrt(E2[a_lo:a_hi+1].sum())
        x_b  = 0.5 * (L + U_eff)
        if y_b > 0:
            x.append(x_b); y.append(y_b); ey.append(ey_b)

    return np.array(x), np.array(y), np.array(ey)'''



def extract_graph_points(gr):
    """
    Return numpy arrays (x, y, ex, ey) from a ROOT.TGraphErrors.
    """
    n = gr.GetN()
    xs, ys, exs, eys = [], [], [], []
    x = array('d', [0.0])
    y = array('d', [0.0])
    for i in range(n):
        gr.GetPoint(i, x, y)        # fills x[0], y[0]
        xs.append(float(x[0]))
        ys.append(float(y[0]))
        exs.append(float(gr.GetErrorX(i)))
        eys.append(float(gr.GetErrorY(i)))
    return np.array(xs), np.array(ys), np.array(exs), np.array(eys)

def export_graph_to_csv(gr, out_path):
    """
    Write angle_deg, norm_counts, ex, ey to CSV (header included).
    """
    xs, ys, exs, eys = extract_graph_points(gr)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["angle_deg", "norm_counts", "ex", "ey"])
        for i in range(len(xs)):
            w.writerow([f"{xs[i]:.6g}", f"{ys[i]:.6g}", f"{exs[i]:.6g}", f"{eys[i]:.6g}"])


def export_graph_to_txt(gr, out_path, header=True, precision=6, delimiter="  "):
    """
    Export graph points to a plain text file (space-delimited by default).
    Columns: x  y  ex  ey
    """
    xs, ys, exs, eys = extract_graph_points(gr)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    hdr = "x  y  ex  ey" if header else ""
    fmt = f"%.{precision}f"
    arr = np.column_stack([xs, ys, exs, eys])

    # np.savetxt doesn’t accept a custom delimiter for headers, which is fine.
    np.savetxt(out_path, arr, header=hdr, comments="" if header else "", fmt=fmt, delimiter=delimiter)


def build_filename_from_yaml(prefix, run_range_label, info, sample_name):
    outID = 'analysis_output_additional_identifiers'
    if not info[outID] or not info[outID]['enabled']:
        # No identifiers at all
        return f"{prefix}{run_range_label}.root"

    join_with     = info[outID]['join_with']
    decimal_sub   = info[outID]['decimal_sub']
    float_fmt_def = info[outID]['float_fmt_default']
    fields        = info[outID]['fields']
    order        = list(fields.keys())

    parts = [f"{prefix}{sample_name}_{run_range_label}"]


    for name in order:
        nameID = info[outID]['fields'][name]
        if not nameID or not info[outID]['fields'][name]['include']:
            continue
        label = info[outID]['fields'][name]['label']
        vals = []
        for v in info[outID]['fields'][name]['values']:
            # if it's a variable name defined in your script, use its value
            if isinstance(v, str) and v in globals():
                vals.append(globals()[v])
            else:
                vals.append(v)  # use literal value from YAML (number or string)

        # Format: remove trailing .0, optionally replace '.' with 'p'
        val_str = "_".join(
            str(float(v)).rstrip("0").rstrip(".") if isinstance(v, (int, float)) else str(v)
            for v in vals
        )
        if decimal_sub:
            val_str = val_str.replace(".", decimal_sub)

        parts.append(f"{label}{val_str}")

    #print(parts)

    return join_with.join(parts) + ".root"

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


# =============================================== PARAMETERS  ===============================================#
sample_name = info['sample_name']
do_ngAnalysis = info['do_ngAnalysis']
do_ggAnalysis = info['do_ggAnalysis']
stage1_hist_filepath = info['stage1_hist_filepath']
stage1_hist_prefix = info['stage1_hist_prefix']
stage1_hist_suffix = info['stage1_hist_suffix']
analysis_output_file_prefix = info['analysis_output_file_prefix']
#stage1_cfg_filepath = info['stage1_cfg_filepath']
summedhists_outputfile_prefix = info['summedhists_outputfile_prefix']
if info.get(f'{sample_name}_init_params') and f'{sample_name}_init_params' in info:
    init_params = info[f'{sample_name}_init_params'] # a list
else:
    init_params = [0,0,0] #default
numParams = int(info['numParams'])

fit_limit_ang_low = float(info['fit_limits_deg']['low'])
fit_limit_ang_high = float(info['fit_limits_deg']['high'])
fit_limit_cos_low = float(info['fit_limits_cos']['low'])
fit_limit_cos_high = float(info['fit_limits_cos']['high'])
fit_limit_cos_sq_low = float(info['fit_limits_cos_sq']['low'])
fit_limit_cos_sq_high = float(info['fit_limits_cos_sq']['high'])
min_accepted_angle = fit_limit_ang_low
max_accepted_angle = fit_limit_ang_high
gg_Esum_cut_low = info['gg_Esum_cuts'][f'{sample_name}_low'] 
gg_Esum_cut_high = info['gg_Esum_cuts'][f'{sample_name}_high'] 
if info.get('ng_Esum_cuts') and f'{sample_name}_low' in info['ng_Esum_cuts']:
    ng_Esum_cut_low = info['ng_Esum_cuts'][f'{sample_name}_low'] 
else:
    ng_Esum_cut_low = None
if info.get('ng_Esum_cuts') and f'{sample_name}_high' in info['ng_Esum_cuts']:
    ng_Esum_cut_high = info['ng_Esum_cuts'][f'{sample_name}_high'] 
else: 
    ng_Esum_cut_high = None
#print("\nng Esum cuts: ", ng_Esum_cut_low, ng_Esum_cut_high)
#print("gg Esum cuts: ", gg_Esum_cut_low, gg_Esum_cut_high)
angleBin_width = info['angleBin_width']
angleRebin_width = info['angleRebin_width']
angle_err = info['angle_err']
theta = info['detector_spherical_coords']['theta'] 
phi = info['detector_spherical_coords']['phi'] 

# =============================================== GET FILEPATHS  ===============================================#

# initialize lists/dictionaries
# -- geometry normalization values -- #
angle_counts_gg = compute_detector_pair_angle_counts(theta, phi, gg = True, ng = False) 
#print("angle_counts_gg: \n", angle_counts_gg)
scaling_values_gg = compute_norm_scaling_values(theta, phi, angleBin_width, gg = True, ng = False) 
scaling_values_ng = compute_norm_scaling_values(theta, phi, angleBin_width, gg = False, ng = True) 
print(f"\n\nscaling_values_gg: \n{scaling_values_gg}")
#print(f"\n\nscaling_values_gg: \n{scaling_values_gg} \n\nscaling_values_ng: \n {scaling_values_ng}")
angle_SumofWeights = defaultdict(float)
angle_SumofWeightsSqd = defaultdict(float)
counts = defaultdict(float)
weight = defaultdict(float)

filestosum = []
angles_runnum = {}

# Determine if we should use summed histograms
if len(run_ranges) > 1:
    use_summed_hists = True
elif run_ranges[0][0] != run_ranges[0][1]:
    use_summed_hists = True
else:
    use_summed_hists = False

# Decide how to get histograms
summedhists_outputfile = f"{summedhists_outputfile_prefix}{sample_name}_{run_range_label}_angleBinWidth{angleBin_width}.root"
if use_summed_hists and os.path.exists(summedhists_outputfile):
    print(f"\nSummed histogram file already exists: {summedhists_outputfile}")
    user_input = input(Fore.CYAN + "Would you like to use this summed histogram file? (y/n): ").strip().lower()

    if user_input == 'y':
        print(Fore.RESET)
        print("Using existing summed histogram.")
    else:
        print(Fore.RESET)
        print("Rebuilding summed histogram from individual run files...")
        for startRun, endRun in run_ranges:
            for runnum in range(startRun, endRun + 1):
                rootfiletosum = f"{stage1_hist_filepath}{stage1_hist_prefix}{runnum}{stage1_hist_suffix}"
                filestosum.append(rootfiletosum)
else:
    for startRun, endRun in run_ranges:
        for runnum in range(startRun, endRun + 1):
            rootfiletosum = f"{stage1_hist_filepath}{stage1_hist_prefix}{runnum}{stage1_hist_suffix}"
            filestosum.append(rootfiletosum)

# If we're summing, run the hadd command
if use_summed_hists and filestosum:
    print("\nSumming histograms . . . ")
    add_hists_command = f"hadd -f {summedhists_outputfile} {' '.join(filestosum)}"
    ROOT.gSystem.Exec(add_hists_command)

# Load the ROOT file
if use_summed_hists:
    rootfile = ROOT.TFile(summedhists_outputfile)
    print(f"\nStarting analysis on: {summedhists_outputfile}")
else:
    rootfilepath = f"{stage1_hist_filepath}{stage1_hist_prefix}{startRun}{stage1_hist_suffix}"
    rootfile = ROOT.TFile(rootfilepath)
    print(f"\nStarting analysis on: {rootfilepath}")

# Create output file
analysis_output_filename = build_filename_from_yaml(analysis_output_file_prefix, run_range_label, info, sample_name)
#print("\nanalysis output file name: ", analysis_output_filename)


#analysis_output_filename = f'{analysis_output_file_prefix}{run_range_label}_angleBinWidth{angleBin_width}_fitRange{fit_limit_ang_low}_{fit_limit_ang_high}.root'
analysis_output_file = ROOT.TFile(analysis_output_filename, "RECREATE")







# =============================================== ANALYSIS  ===============================================#

# what type of analysis are you interested in doing? n-g, g-g, or both? see yaml file to change setup
if not do_ngAnalysis and not do_ggAnalysis:
    print(f"No analysis type chosen. See {yaml_file} to change setup.")


#>>>>>>>>>>>>>>>>>>>>>>>>> g-g Angular Correlation <<<<<<<<<<<<<<<<<<<<<<<<<#
if do_ggAnalysis:

    #-------------------------- GET HISTS --------------------------#
    gg_hists = []

    #find hists with keywords
    for key in rootfile.GetListOfKeys():
        kname = key.GetName()
        name = rootfile.Get(kname)
        if 'TGraph' not in str(type(name)) and kname.startswith("ggAngle"):
            gg_hists.append(name)

        #get most common crystals per cluster
        if kname.startswith("hClusterSize"):
            hClusterSize = name
    maxbin = hClusterSize.GetMaximumBin()
    most_common = hClusterSize.GetXaxis().GetBinCenter(maxbin)
    count = hClusterSize.GetBinContent(maxbin)
    print(f"\nMost common crystals per cluster = {most_common:.0f} (count = {int(count)})\n")


    #--------------------- FILL TGRAPHERRORS ---------------------#
    for hist in gg_hists:
        xaxis = hist.GetXaxis()
        yaxis = hist.GetYaxis()
        #scalefactor_bbb = geom_scaling_bbb(theta, phi, xaxis, gg = True, ng = False) 
        #scale_by_binx, scale_by_center = aggregate_scaling_to_hist_bins(scaling_values_gg, xaxis)
        ybin_lo = yaxis.FindBin(gg_Esum_cut_low + 1e-9)   # tiny offset to avoid edge issues
        ybin_hi = yaxis.FindBin(gg_Esum_cut_high - 1e-9)
        #print(f"Nbins = {hist.GetNbinsX()}")
        x_ang, x_cos, x_cos_sq, y, ey  = [] , [],  [], [], []
        for binx in range(1, hist.GetNbinsX() + 1):
            # ----- bin check ---- #
            binlow = xaxis.GetBinLowEdge(binx)
            bincenter = xaxis.GetBinCenter(binx)
            binhigh = xaxis.GetBinUpEdge(binx)
            #print(f"Bin {binx}: [{binlow:.2f}, {binhigh:.2f}] center={bincenter:.2f}")
            # --------------------- #
            
            angle = xaxis.GetBinCenter(binx)
            #print(f"xaxis.GetBinCenter({binx}) = {angle}")
            cos_theta = math.cos(math.radians(angle))
            sin_theta = math.sin(math.radians(angle))
            cos_theta_sq = cos_theta * cos_theta
            cos_err = abs(sin_theta) * math.radians(angle_err)
            cos_sq_err = 2 * abs(cos_theta * sin_theta) * math.radians(angle_err)

            '''angle = bin_xvalue(xaxis, binx, mode="angle")
            cos_theta = bin_xvalue(xaxis, binx, mode="cos")
            cos_theta_sq = bin_xvalue(xaxis, binx, mode="cos2")
            cos_err = 0.05
            cos_sq_err = 0.05'''


            if angle >= min_accepted_angle and angle <= max_accepted_angle:
                sum_counts = 0.0
                sum_counts_err2 = 0.0
                for biny in range(ybin_lo, ybin_hi + 1):
                    val =  hist.GetBinContent(binx, biny)
                    err = hist.GetBinError(binx,biny)
                    sum_counts += val
                    sum_counts_err2 += err * err

                if sum_counts <= 0:
                    continue
                angle
                scalefactor = scaling_values_gg[angle]
                #scalefactor = scalefactor_bbb.get(binx)
                #scalefactor = scale_by_binx.get(binx)
                if scalefactor is None:
                    continue
                scaled_counts = sum_counts * abs(scalefactor)
                # error propagation: (σ_scaled)^2 = (s^2)*(σ_sum)^2 (+ optional scale term)
                #scaled_counts_err2 = (scalefactor**2) * sum_counts_err2
                scaled_counts_err = math.sqrt(sum_counts_err2) * abs(scalefactor)
                # include uncertainty on the scale factor here if you have it
                '''if angle in scaling_errs:
                    ds = scaling_errs[angle]
                    # add (∂/∂s)^2 term: (sum_counts * ds)^2
                    scaled_counts_err2 += (sum_counts * ds) ** 2'''
                



                #print(f"angle={angle}, counts={int(sum_counts)}, scalefactor={scalefactor:.4f}, scaled_counts={scaled_counts:.3f}, scaled_counts_err={scaled_counts_err:.3f}")
                x_ang.append(angle)
                x_cos.append(cos_theta)
                x_cos_sq.append(cos_theta_sq)
                y.append(scaled_counts)
                ey.append(scaled_counts_err)
                #x_deg, y_corr, y_err = correct_1deg_then_sum_to_bins(hist, gg_Esum_cut_low, gg_Esum_cut_high, theta, phi, gg=True, ng=False, coarse_width=angleRebin_width)

        
        n_ang = len(x_ang)
        ggAngle_normcounts= ROOT.TGraphErrors(n_ang)

        #n_ang_test = len(x_deg)
        ##test_ggAngle_normcounts = ROOT.TGraphErrors(n_ang_test)

        for i in range(n_ang):
            ggAngle_normcounts.SetPoint(i, x_ang[i], y[i])
            ggAngle_normcounts.SetPointError(i, angle_err, ey[i])
        #for i in range(n_ang_test):
            #test_ggAngle_normcounts.SetPoint(i, x_deg[i], y_corr[i])
            #test_ggAngle_normcounts.SetPointError(i, angle_err, y_err[i])
        
        
        n_cos = len(x_cos)
        ggCosAngle_normcounts = ROOT.TGraphErrors(n_cos)
        for i in range(n_cos):
            ggCosAngle_normcounts.SetPoint(i, x_cos[i], y[i])
            ggCosAngle_normcounts.SetPointError(i, cos_err, ey[i])

        n_cos_sq = len(x_cos_sq)
        ggCos2Angle_normcounts = ROOT.TGraphErrors(n_cos_sq)
        for i in range(n_cos_sq):
            ggCos2Angle_normcounts.SetPoint(i, x_cos_sq[i], y[i])
            ggCos2Angle_normcounts.SetPointError(i, cos_sq_err, ey[i])

        ending = hist.GetName()[len("ggAngle_Esum_"):]

        # counts vs angle
        analysis_output_file.cd()
        canvas_ang = ROOT.TCanvas()
        #test_ggAngle_normcounts.SetName(f"#test_ggAngle_normcounts_{ending}")
        ggAngle_normcounts.SetName(f"ggAngle_normcounts_{ending}")
        #test_ggAngle_normcounts.SetTitle("#gamma-#gamma angular correlation; #Delta#sigma_{#gamma#gamma}(deg);Counts (geometry-normalized)")
        ggAngle_normcounts.SetTitle("#gamma-#gamma angular correlation; #Delta#sigma_{#gamma#gamma}(deg);Counts (geometry-normalized)")
        ggAngle_normcounts.SetMarkerStyle(20)
        #print(f"number of points in {ggAngle_normcounts.GetName()}: {ggAngle_normcounts.GetN()}")
        ggAngle_normcounts.Draw("AP")
        #test_ggAngle_normcounts.Draw("AP")
        canvas_ang.Update()
        canvas_ang.Write(f"ggAngle_normcounts_{ending}")
        #canvas_ang.Write(f"#test_ggAngle_normcounts_{ending}")


        
        # counts vs cos(angle)
        analysis_output_file.cd()
        canvas_cos = ROOT.TCanvas()
        ggCosAngle_normcounts.SetName(f"ggCosAngle_normcounts_{ending}")
        ggCosAngle_normcounts.SetTitle("#gamma-#gamma angular correlation;cos(#Delta#sigma_{#gamma#gamma});Counts (geometry-normalized)")
        ggCosAngle_normcounts.SetMarkerStyle(20)
        ggCosAngle_normcounts.Draw("AP")
        canvas_cos.Update()
        canvas_cos.Write(f"ggCosAngle_normcounts_{ending}")
        
        # counts vs cos(angle)
        analysis_output_file.cd()
        canvas_cos_sq = ROOT.TCanvas()
        ggCos2Angle_normcounts.SetName(f"ggCos2Angle_normcounts_{ending}")
        ggCos2Angle_normcounts.SetTitle("#gamma-#gamma angular correlation;cos^{2}(#Delta#sigma_{#gamma#gamma});Counts (geometry-normalized)")
        ggCos2Angle_normcounts.SetMarkerStyle(20)
        ggCos2Angle_normcounts.Draw("AP")
        canvas_cos_sq.Update()
        canvas_cos_sq.Write(f"ggCos2Angle_normcounts_{ending}")

        # export data to csv file
        export_txt = info['export_txt']['enabled']
        if export_txt:
            graph_registry = {
            "ggAngle_normcounts":   ggAngle_normcounts,
            "ggCosAngle_normcounts": ggCosAngle_normcounts,
            "ggCos2Angle_normcounts": ggCos2Angle_normcounts,
            # add more here if needed
            }
            graphs_to_txt = info['export_txt']['graphs_to_txt']
            txt_dir = "AngularAnalysis_ExportTXT"
            for gname in graphs_to_txt:
                gr = graph_registry.get(gname)
                gtitle = gr.GetName()
                if gr is None:
                    print(f"Skipping '{gname}': not found in registry.")
                    continue
                out_file = f"{sample_name}_{run_range_label}_{gtitle}.txt"
                out_path = os.path.join(txt_dir, out_file)
                export_graph_to_txt(gr, out_path)




        #----------------------------- FITTING -----------------------------#
        fitAng = ROOT.TF1("fitAng", fitFunc_ang, fit_limit_ang_low, fit_limit_ang_high, numParams)
        fitAng.SetParNames("a0", "a2", "a4")
        fitAng.SetParameters(init_params[0], init_params[1], init_params[2])
        fitAng.SetLineColor(ROOT.kRed)
        ggAngle_normcounts.Fit(fitAng, "R")
        #test_ggAngle_normcounts.Fit(fitAng, "R")
        print("\nNormalized fit parameters (unscaled)")
        print("a0 =", fitAng.GetParameter(0))
        print("a2 =", fitAng.GetParameter(1))
        print("a4 =", fitAng.GetParameter(2))

        fitCos = ROOT.TF1("fitCos", fitFunc_cos, fit_limit_cos_low, fit_limit_cos_high, numParams)
        fitCos.SetParNames("a0", "a2", "a4")
        fitCos.SetParameters(init_params[0], init_params[1], init_params[2])
        fitCos.SetLineColor(ROOT.kRed)
        ggCosAngle_normcounts.Fit(fitCos, "R")
        print("\nNormalized fit parameters (unscaled)")
        print("a0 =", fitCos.GetParameter(0))
        print("a2 =", fitCos.GetParameter(1))
        print("a4 =", fitCos.GetParameter(2))

        fitCos_sq = ROOT.TF1("fitCos_sq", fitFunc_cos_sq, fit_limit_cos_sq_low, fit_limit_cos_sq_high, numParams)
        fitCos_sq.SetParNames("a0", "a2", "a4")
        fitCos_sq.SetParameters(init_params[0], init_params[1], init_params[2])
        fitCos_sq.SetLineColor(ROOT.kRed)
        ggCos2Angle_normcounts.Fit(fitCos_sq, "R")
        print("\nNormalized fit parameters (unscaled)")
        print("a0 =", fitCos_sq.GetParameter(0))
        print("a2 =", fitCos_sq.GetParameter(1))
        print("a4 =", fitCos_sq.GetParameter(2))

        
        analysis_output_file.cd()
        
        #-------------- SCALING PARAMS AND DRAWNG FIT --------------#

        #creating box to display parameters on canvas
        box = ROOT.TPaveText(0.62, 0.58, 0.89, 0.89, "NDC")
        box.SetFillColor(0)
        box.SetFillStyle(1001)
        box.SetBorderSize(1)
        box.SetTextFont(42)
        box.SetTextSize(0.03)


        # ➤ geometry-normalized counts vs angle
        cfit_ang = ROOT.TCanvas(f"fit_ggAngle_normcounts_{ending}", "g-g Angular Correlation Fit")
        ggAngle_normcounts.SetTitle("#gamma-#gamma Angular Correlation;#Delta#sigma_{#gamma#gamma} (degrees);Counts (geometry normalized)")
        ggAngle_normcounts.GetXaxis().SetRangeUser(fit_limit_ang_low, fit_limit_ang_high)
        ggAngle_normcounts.SetMarkerStyle(20)
        ggAngle_normcounts.Draw("AP")
        fitAng.Draw("same")
        print("\nFinal normalized and scaled fit parameters for counts vs angle:")
        # scaling 
        box.AddText("W(#Delta#sigma)/a_{0}")
        for i in range(fitAng.GetNpar()):
            name = fitAng.GetParName(i) or f"p{i}"
            val  = fitAng.GetParameter(i) / fitAng.GetParameter(0)
            err  = val * math.sqrt((fitAng.GetParError(i)/fitAng.GetParameter(i))**2 + (fitAng.GetParError(0)/fitAng.GetParameter(0))**2)
            print(f"a{i} = {val}, error = {err}")
            box.AddText(f"{name} = {val:.3g} #pm {err:.3g}")
        chi2 = fitAng.GetChisquare()
        ndf  = fitAng.GetNDF()
        if ndf > 0:
            box.AddText(f"#chi^{{2}}/NDF = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
        else:
            box.AddText(f"#chi^{{2}} = {chi2:.1f}")
        box.Draw()
        cfit_ang.Update()
        cfit_ang.Write(f"fit_ggAngle_normcounts_{ending}")
        
        # ➤ geometry-normalized counts vs cos(angle)
        cfit_cos = ROOT.TCanvas(f"fit_ggCosAngle_normcounts_{ending}", "g-g Angular Correlation Fit")
        ggCosAngle_normcounts.SetTitle("#gamma-#gamma Angular Correlation;cos(#Delta#sigma_{#gamma#gamma});Counts (geometry-normalized)")
        ggCosAngle_normcounts.GetXaxis().SetRangeUser(fit_limit_cos_low, fit_limit_cos_high)
        ggCosAngle_normcounts.SetMarkerStyle(20)
        ggCosAngle_normcounts.Draw("AP")
        fitCos.Draw("same")
        print("\nFinal normalized and scaled fit parameters for counts vs cos(angle):")
        # scaling 
        box.Clear()
        box.AddText("W(#Delta#sigma)/a_{0}")
        for i in range(fitCos.GetNpar()):
            name = fitCos.GetParName(i) or f"p{i}"
            val  = fitCos.GetParameter(i) / fitCos.GetParameter(0)
            err  = val * math.sqrt((fitCos.GetParError(i)/fitCos.GetParameter(i))**2 + (fitCos.GetParError(0)/fitCos.GetParameter(0))**2)
            print(f"a{i} = {val}, error = {err}")
            box.AddText(f"{name} = {val:.3g} #pm {err:.3g}")
        chi2 = fitCos.GetChisquare()
        ndf  = fitCos.GetNDF()
        if ndf > 0:
            box.AddText(f"#chi^{{2}}/NDF = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
        else:
            box.AddText(f"#chi^{{2}} = {chi2:.1f}")
        box.Draw()
        cfit_cos.Update()
        cfit_cos.Write(f"fit_ggCosAngle_normcounts_{ending}")

        # ➤ geometry-normalized counts vs cos^2(angle)
        cfit_cos_sq = ROOT.TCanvas(f"fit_ggCos2Angle_normcounts_{ending}", "g-g Angular Correlation Fit")
        ggCos2Angle_normcounts.SetTitle("#gamma-#gamma Angular Correlation;cos^{2}(#Delta#sigma_{#gamma#gamma});Counts (geometry-normalized)")
        ggCos2Angle_normcounts.GetXaxis().SetRangeUser(fit_limit_cos_sq_low, fit_limit_cos_sq_high)
        ggCos2Angle_normcounts.SetMarkerStyle(20)
        ggCos2Angle_normcounts.Draw("AP")
        fitCos_sq.Draw("same")
        print("\nFinal normalized and scaled fit parameters for counts vs cos^2(angle):")
        # scaling 
        box.Clear()
        box.AddText("W(#Delta#sigma)/a_{0}")
        for i in range(fitCos_sq.GetNpar()):
            name = fitCos_sq.GetParName(i) or f"p{i}"
            val  = fitCos_sq.GetParameter(i) / fitCos_sq.GetParameter(0)
            err  = val * math.sqrt((fitCos_sq.GetParError(i)/fitCos_sq.GetParameter(i))**2 + (fitCos_sq.GetParError(0)/fitCos_sq.GetParameter(0))**2)
            print(f"a{i} = {val}, error = {err}")
            box.AddText(f"{name} = {val:.3g} #pm {err:.3g}")
        chi2 = fitCos_sq.GetChisquare()
        ndf  = fitCos_sq.GetNDF()
        if ndf > 0:
            box.AddText(f"#chi^{{2}}/NDF = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
        else:
            box.AddText(f"#chi^{{2}} = {chi2:.1f}")
        box.Draw()
        cfit_cos_sq.Update()
        cfit_cos_sq.Write(f"fit_ggCos2Angle_normcounts_{ending}")


        print(Fore.GREEN 
                + "------------------------------------------" 
                +f"\n            Analysis complete.\n" 
                + "------------------------------------------"  
                +Fore.RESET + " \n\nResults saved to: " 
                + Fore.WHITE + f"{analysis_output_file.GetName()}.")
        
        print(Fore.RESET)

#>>>>>>>>>>>>>>>>>>>>>>>>> n-g Angular Distribution <<<<<<<<<<<<<<<<<<<<<<<<<#
if do_ngAnalysis:

    #-------------------------- GET HISTS --------------------------#
    ng_hists = []


    #find hists with keywords
    for key in rootfile.GetListOfKeys():
        kname = key.GetName()
        name = rootfile.Get(kname)
        if 'TGraph' not in str(type(name)) and kname.startswith("ngAngle"):
            ng_hists.append(name)

        #get most common crystals per cluster
        '''if kname.startswith("hClusterSize"):
            hClusterSize = name
    maxbin = hClusterSize.GetMaximumBin()
    most_common = hClusterSize.GetXaxis().GetBinCenter(maxbin)
    count = hClusterSize.GetBinContent(maxbin)
    print(f"Most common crystals per cluster = {most_common:.0f} (count = {int(count)})")'''



    for hist in ng_hists:
        xaxis = hist.GetXaxis()
        yaxis = hist.GetYaxis()

        if ng_Esum_cut_low:
            ybin_lo = yaxis.FindBin(ng_Esum_cut_low + 1e-9)   # tiny offset to avoid edge issues
        else:
            ybin_lo = 1
        if ng_Esum_cut_high:
            ybin_hi = yaxis.FindBin(ng_Esum_cut_high - 1e-9)
        else:
            ybin_hi = yaxis.GetNbins()

        x_ang, x_cos, x_cos_sq, y, ey  = [] , [],  [], [], []
        for binx in range(1, hist.GetNbinsX() + 1):
            #angle = round(xaxis.GetBinCenter(binx) / angleBin_width) * angleBin_width
            angle = int(xaxis.GetBinCenter(binx))
            cos_theta = math.cos(math.radians(angle))
            sin_theta = math.sin(math.radians(angle))
            cos_theta_sq = cos_theta * cos_theta
            cos_err = abs(sin_theta) * math.radians(angle_err)
            cos_sq_err = 2 * abs(cos_theta * sin_theta) * math.radians(angle_err)
            if angle >= min_accepted_angle and angle <= max_accepted_angle:
                sum_counts = 0.0
                sum_counts_err2 = 0.0
                for biny in range(ybin_lo, ybin_hi + 1):
                    val =  hist.GetBinContent(binx, biny)
                    err = hist.GetBinError(binx,biny)
                    sum_counts += val
                    sum_counts_err2 += err**2

                if sum_counts <= 0:
                    continue

                scalefactor = scaling_values_ng[angle]
                scaled_counts = sum_counts * scalefactor
                # error propagation: (σ_scaled)^2 = (s^2)*(σ_sum)^2 (+ optional scale term)

                scaled_counts_err2 = (scalefactor**2) * sum_counts_err2
                scaled_counts_err = math.sqrt(scaled_counts_err2)
                # include uncertainty on the scale factor here if you have it
                '''if angle in scaling_errs:
                    ds = scaling_errs[angle]
                    # add (∂/∂s)^2 term: (sum_counts * ds)^2
                    scaled_counts_err2 += (sum_counts * ds) ** 2'''
                


                #print(f"angle={angle}, counts={int(sum_counts)}, scalefactor={scalefactor:.4f}, scaled_counts={scaled_counts:.3f}, scaled_counts_err={scaled_counts_err:.3f}")
                x_ang.append(angle)
                x_cos.append(cos_theta)
                x_cos_sq.append(cos_theta_sq)
                y.append(scaled_counts)
                ey.append(scaled_counts_err)
        
        n_ang = len(x_ang)
        ngAngle_normcounts= ROOT.TGraphErrors(n_ang)
        for i in range(n_ang):
            ngAngle_normcounts.SetPoint(i, x_ang[i], y[i])
            ngAngle_normcounts.SetPointError(i, angle_err, ey[i])

        n_cos = len(x_cos)
        ngCosAngle_normcounts = ROOT.TGraphErrors(n_cos)
        for i in range(n_cos):
            ngCosAngle_normcounts.SetPoint(i, x_cos[i], y[i])
            ngCosAngle_normcounts.SetPointError(i, cos_err, ey[i])

        n_cos_sq = len(x_cos_sq)
        ngCos2Angle_normcounts = ROOT.TGraphErrors(n_cos_sq)
        for i in range(n_cos_sq):
            ngCos2Angle_normcounts.SetPoint(i, x_cos_sq[i], y[i])
            ngCos2Angle_normcounts.SetPointError(i, cos_sq_err, ey[i])

        ending = hist.GetName()[len("ngAngle_Esum_"):]

        
        # counts vs angle
        analysis_output_file.cd()
        canvas_ang = ROOT.TCanvas()
        ngAngle_normcounts.SetName(f"ngAngle_normcounts_{ending}")
        ngAngle_normcounts.SetTitle("n-#gamma angular distribution; #Delta#sigma_{n#gamma}(deg);Counts (geometry-normalized)")
        ngAngle_normcounts.SetMarkerStyle(20)
        #print(f"number of points in {ngAngle_normcounts.GetName()}: {ngAngle_normcounts.GetN()}")
        ngAngle_normcounts.Draw("AP")
        canvas_ang.Update()
        canvas_ang.Write(f"ngAngle_normcounts_{ending}")
        
        # counts vs cos(angle)
        analysis_output_file.cd()
        canvas_cos = ROOT.TCanvas()
        ngCosAngle_normcounts.SetName(f"ngCosAngle_normcounts_{ending}")
        ngCosAngle_normcounts.SetTitle("n-#gamma angular distribution;cos(#Delta#sigma_{n#gamma});Counts (geometry-normalized)")
        ngCosAngle_normcounts.SetMarkerStyle(20)
        ngCosAngle_normcounts.Draw("AP")
        canvas_cos.Update()
        canvas_cos.Write(f"ngCosAngle_normcounts_{ending}")
        
        # counts vs cos(angle)
        analysis_output_file.cd()
        canvas_cos_sq = ROOT.TCanvas()
        ngCos2Angle_normcounts.SetName(f"ngCos2Angle_normcounts_{ending}")
        ngCos2Angle_normcounts.SetTitle("n-#gamma angular distribution;cos^{2}(#Delta#sigma_{n#gamma});Counts (geometry-normalized)")
        ngCos2Angle_normcounts.SetMarkerStyle(20)
        ngCos2Angle_normcounts.Draw("AP")
        canvas_cos_sq.Update()
        canvas_cos_sq.Write(f"ngCos2Angle_normcounts_{ending}")

        # export data to csv file
        export_csv = False
        if export_csv == True:
            csv_dir = "ExportCSV"
            csv_name = f"{sample_name}{ngAngle_normcounts.GetName()}.csv"
            csv_path = os.path.join(csv_dir, csv_name)
            export_graph_to_csv(ngAngle_normcounts, csv_path)
            print(f"Wrote CSV: {csv_path}")



        #----------------------------- FITTING -----------------------------#
        fitAng = ROOT.TF1("fitAng", fitFunc_ang, fit_limit_ang_low, fit_limit_ang_high, numParams)
        fitAng.SetParNames("a0", "a2", "a4")
        fitAng.SetParameters(init_params[0], init_params[1], init_params[2])
        fitAng.SetLineColor(ROOT.kRed)
        ngAngle_normcounts.Fit(fitAng, "R")
        print("\nNormalized fit parameters (unscaled)")
        print("a0 =", fitAng.GetParameter(0))
        print("a2 =", fitAng.GetParameter(1))
        print("a4 =", fitAng.GetParameter(2))

        fitCos = ROOT.TF1("fitCos", fitFunc_cos, fit_limit_cos_low, fit_limit_cos_high, numParams)
        fitCos.SetParNames("a0", "a2", "a4")
        fitCos.SetParameters(init_params[0], init_params[1], init_params[2])
        fitCos.SetLineColor(ROOT.kRed)
        ngCosAngle_normcounts.Fit(fitCos, "R")
        print("\nNormalized fit parameters (unscaled)")
        print("a0 =", fitCos.GetParameter(0))
        print("a2 =", fitCos.GetParameter(1))
        print("a4 =", fitCos.GetParameter(2))

        fitCos_sq = ROOT.TF1("fitCos_sq", fitFunc_cos_sq, fit_limit_cos_sq_low, fit_limit_cos_sq_high, numParams)
        fitCos_sq.SetParNames("a0", "a2", "a4")
        fitCos_sq.SetParameters(init_params[0], init_params[1], init_params[2])
        fitCos_sq.SetLineColor(ROOT.kRed)
        ngCos2Angle_normcounts.Fit(fitCos_sq, "R")
        print("\nNormalized fit parameters (unscaled)")
        print("a0 =", fitCos_sq.GetParameter(0))
        print("a2 =", fitCos_sq.GetParameter(1))
        print("a4 =", fitCos_sq.GetParameter(2))

        
        analysis_output_file.cd()
        
        #-------------- SCALING PARAMS AND DRAWNG FIT --------------#

        #creating box to display parameters on canvas
        box = ROOT.TPaveText(0.62, 0.58, 0.89, 0.89, "NDC")
        box.SetFillColor(0)
        box.SetFillStyle(1001)
        box.SetBorderSize(1)
        box.SetTextFont(42)
        box.SetTextSize(0.03)


        # ➤ geometry-normalized counts vs angle
        cfit_ang = ROOT.TCanvas(f"fit_ngAngle_normcounts_{ending}", "n-g Angular Distribution Fit")
        ngAngle_normcounts.SetTitle("n-#gamma Angular Distribution;#Delta#sigma_{n#gamma} (degrees);Counts (geometry normalized)")
        ngAngle_normcounts.SetMarkerStyle(20)
        ngAngle_normcounts.Draw("AP")
        fitAng.Draw("same")
        print("\nFinal normalized and scaled fit parameters for counts vs angle:")
        # scaling 
        box.AddText("W(#Delta#sigma)/a_{0}")
        for i in range(fitAng.GetNpar()):
            name = fitAng.GetParName(i) or f"p{i}"
            val  = fitAng.GetParameter(i) / fitAng.GetParameter(0)
            err  = val * math.sqrt((fitAng.GetParError(i)/fitAng.GetParameter(i))**2 + (fitAng.GetParError(0)/fitAng.GetParameter(0))**2)
            print(f"a{i} = {val}, error = {err}")
            box.AddText(f"{name} = {val:.3g} #pm {err:.3g}")
        chi2 = fitAng.GetChisquare()
        ndf  = fitAng.GetNDF()
        if ndf > 0:
            box.AddText(f"#chi^{{2}}/NDF = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
        else:
            box.AddText(f"#chi^{{2}} = {chi2:.1f}")
        box.Draw()
        cfit_ang.Update()
        cfit_ang.Write(f"fit_ngAngle_normcounts_{ending}")
        
        # ➤ geometry-normalized counts vs cos(angle)
        cfit_cos = ROOT.TCanvas(f"fit_ngCosAngle_normcounts_{ending}", "n-g Angular Distribution Fit")
        ngCosAngle_normcounts.SetTitle("n-#gamma Angular Distribution;cos(#Delta#sigma_{n#gamma}) ;Counts (geometry-normalized)")
        ngCosAngle_normcounts.SetMarkerStyle(20)
        ngCosAngle_normcounts.Draw("AP")
        fitCos.Draw("same")
        print("\nFinal normalized and scaled fit parameters for counts vs cos(angle):")
        # scaling 
        box.Clear()
        box.AddText("W(#Delta#sigma)/a_{0}")
        for i in range(fitCos.GetNpar()):
            name = fitCos.GetParName(i) or f"p{i}"
            val  = fitCos.GetParameter(i) / fitCos.GetParameter(0)
            err  = val * math.sqrt((fitCos.GetParError(i)/fitCos.GetParameter(i))**2 + (fitCos.GetParError(0)/fitCos.GetParameter(0))**2)
            print(f"a{i} = {val}, error = {err}")
            box.AddText(f"{name} = {val:.3g} #pm {err:.3g}")
        chi2 = fitCos.GetChisquare()
        ndf  = fitCos.GetNDF()
        if ndf > 0:
            box.AddText(f"#chi^{{2}}/NDF = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
        else:
            box.AddText(f"#chi^{{2}} = {chi2:.1f}")
        box.Draw()
        cfit_cos.Update()
        cfit_cos.Write(f"fit_ngCosAngle_normcounts_{ending}")

        # ➤ geometry-normalized counts vs cos^2(angle)
        cfit_cos_sq = ROOT.TCanvas(f"fit_ngCos2Angle_normcounts_{ending}", "n-g Angular Distribution Fit")
        ngCos2Angle_normcounts.SetTitle("n-#gamma Angular Distribution;cos^{2}(#Delta#sigma_{n#gamma});Counts (geometry-normalized)")
        ngCos2Angle_normcounts.SetMarkerStyle(20)
        ngCos2Angle_normcounts.Draw("AP")
        fitCos_sq.Draw("same")
        print("\nFinal normalized and scaled fit parameters for counts vs cos^2(angle):")
        # scaling 
        box.Clear()
        box.AddText("W(#Delta#sigma)/a_{0}")
        for i in range(fitCos_sq.GetNpar()):
            name = fitCos_sq.GetParName(i) or f"p{i}"
            val  = fitCos_sq.GetParameter(i) / fitCos_sq.GetParameter(0)
            err  = val * math.sqrt((fitCos_sq.GetParError(i)/fitCos_sq.GetParameter(i))**2 + (fitCos_sq.GetParError(0)/fitCos_sq.GetParameter(0))**2)
            print(f"a{i} = {val}, error = {err}")
            box.AddText(f"{name} = {val:.3g} #pm {err:.3g}")
        chi2 = fitCos_sq.GetChisquare()
        ndf  = fitCos_sq.GetNDF()
        if ndf > 0:
            box.AddText(f"#chi^{{2}}/NDF = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
        else:
            box.AddText(f"#chi^{{2}} = {chi2:.1f}")
        box.Draw()
        cfit_cos_sq.Update()
        cfit_cos_sq.Write(f"fit_ngCos2Angle_normcounts_{ending}")


        print(f"Fit complete. Results saved to {analysis_output_file}.")


  



