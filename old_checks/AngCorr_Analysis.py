#JRP
import ROOT
import numpy as np
import math
from array import array
import sys
from collections import defaultdict



# ****************************** begin user input  ****************************** #

#run settings
make_histograms = True
sumhists = False
fix_params = False
do_fit = True
#esum windows
lowEsum_cut = 2.4
highEsum_cut = 2.7
#initial parameters
init_params = [1, 1/8, 1/24]
#normvalue = 19.43 # arbitrary overall normalization value obtained from fit value p3 , 1.0 for placeholder
#fit limits
fit_limit_low = 30
fit_limit_high = 180
# avg estimated angle error
# about 16 degrees between nearest neighbors: 16/2 = 8 degrees to the edge of the crystal
angle_err = 8.0

# ****************************** end user input  ****************************** #


if len(sys.argv) < 3:
    print("Proper usage is:  python3 AngCorr_Analysis.py startingRun# endingRun#") 
    exit(-1)
startRun = int(str(sys.argv[1]))
endRun = int(str(sys.argv[2]))

if startRun == endRun:
    use_summed_hists = False
else:
    use_summed_hists = True

file_path = 'detector_angles_matrix.txt'
matrix = np.loadtxt(file_path)

def ang_distrib_func(x,params):
    theta_deg = x[0]
    cos_theta = ROOT.TMath.Cos(theta_deg * ROOT.TMath.DegToRad())  # Convert to radians

    #Legendre polynomials
    P2 = 0.5 * (3 * cos_theta**2 - 1)  # P2(cos θ)
    P4 = 0.125 * (35 * cos_theta**4 - 30 * cos_theta**2 + 3)  # P4(cos θ)

    # Angular distribution: W(θ) = norm * (1 + a2 * P2 + a4 * P4)
    a2 = params[0]
    a4 = params[1]

    return 1 + (a2 * P2) + (a4 * P4)


# angular correlation function where x = angle
def ang_corr_func_ang(x, params):
    theta = x[0]  # Angle in degrees
    cos_theta = ROOT.TMath.Cos(theta * ROOT.TMath.DegToRad())  # Convert to radians

    #Legendre polynomials
    P2 = 0.5 * (3 * cos_theta**2 - 1)  # P2(cos θ)
    P4 = 0.125 * (35 * cos_theta**4 - 30 * cos_theta**2 + 3)  # P4(cos θ)

    return params[3] * (params[0] + (params[1] * (cos_theta)**2) + (params[2] * (cos_theta)**4))

def ang_corr_func_ang_noA(x, params):
    theta = x[0]  # Angle in degrees
    cos_theta = ROOT.TMath.Cos(theta * ROOT.TMath.DegToRad())  # Convert to radians

    #Legendre polynomials
    P2 = 0.5 * (3 * cos_theta**2 - 1)  # P2(cos θ)
    P4 = 0.125 * (35 * cos_theta**4 - 30 * cos_theta**2 + 3)  # P4(cos θ)

    return params[0] + (params[1] * (cos_theta)**2) + (params[2] * (cos_theta)**4)


# angular correlation function where x = cos(θ)
def ang_corr_func_cos(x, params):
    P2 = 0.5 * (3 * x[0]**2 - 1)  # P2(cos θ)
    P4 = 0.125 * (35 * x[0]**4 - 30 * x[0]**2 + 3)  # P4(cos θ)

    return params[3] * (params[0] + (params[1] * (x[0])**2) + (params[2] * (x[0])**4))



# here I included detectors 76 and 86 even though they weren't taking data
# a check showed it'd be better with them included
scaling_values = {}
with open("angle_counts.txt", "r") as file:
    for line in file:
        angle_val, num_pairs = map(float, line.strip().split())
        if num_pairs <= 0:
            continue

        scaling_values[angle_val] = 1 / num_pairs
#print(len(scaling_values), scaling_values) #len=140

with open("scaling_values.txt", "w") as f:
    for angle, scale in sorted(scaling_values.items(), key=lambda kv: kv[0]):
        f.write(f"{angle:.1f} {scale}\n")
# initialize lists/dictionaries
filestosum = []
angles_runnum = {}

# read unique angles
with open("unique_angles.txt", "r") as file:
    unique_angles = [float(line.strip()) for line in file]

# go through run range
for runnum in range(startRun, endRun +1):
    
    # get histogram
    histfile = ROOT.TFile(f"/home/jr2514/DANCE/DANCE_Analysis/stage1_root/treeTest/Stage1_Histograms_Run_{runnum}_10ns_CW_0ns_CBT_0ns_DEBT.root")
    # get eventTree
    eventTree = histfile.Get("eventTree")

    # define output for each run
    outfile_name = f"/data2/lansce/jr2514/AngCorr_Analysis_Output/Analysis_Output_{runnum}.root"
    filestosum.append(outfile_name)


#-------------------- make hists --------------------#  

    if make_histograms == True:
        output_file = ROOT.TFile(outfile_name, "RECREATE")
        print(f"Analyzing Run {runnum} ...")
        # From global.h
        NeutronE_From = 0.005  # Replace with the actual starting value
        NeutronE_To = 5e6     # Replace with the actual ending value
        NeutronE_BinsPerDecade = 200  # Replace with the number of bins per decade
        GammaE_From = 0.0 #Gamma energy [MeV] - low limit
        GammaE_To = 20.0 #Gamma energy [MeV] - upper limit
        GammaE_NoOfBins = 200 #Number of bins

        # Initialize the fixed-size list and NEbins
        x = [0.0] * 5000  # Preallocate a list of length 5000 with default value 0.0
        NEbins = 0
        NoOfEnergyBins = GammaE_NoOfBins
        # Loop through logarithmic values
        lx = math.log10(NeutronE_From)
        while lx < math.log10(NeutronE_To):
            x[NEbins] = 10**lx  # Assign the computed value to the list
            lx += 1.0 / NeutronE_BinsPerDecade
            NEbins += 1

        NEbins -= 1  # Adjust NEbins as in the original C++ code


        Mbins = [0.0] * 21  # Preallocate a list for 21 bins
        EtotBins = [0.0] * (NoOfEnergyBins + 1)  # Preallocate a list for EtotBins

        # Fill Mbins
        for i in range(21):
            Mbins[i] = i  # Optionally, use Mbins[i] = 0.5 + 1.0 * i for different bin edges
        DEGamma = (GammaE_To - GammaE_From) / GammaE_NoOfBins

        # Fill EtotBins
        for i in range(NoOfEnergyBins + 1):
            # Use EtotBins[i] = i * (16.0 / 128.0) for the commented-out option
            EtotBins[i] = GammaE_From + i * DEGamma

        x_array = array('d', x)
        Mbins_array = array('d', Mbins)
        EtotBins_array = array('d', EtotBins)

        # Create the 3D histogram
        En_Esum_Mcl = ROOT.TH3F(
            "En_Esum_Mcl",               # Histogram name
            "En_Esum_Mcl",               # Histogram title
            NEbins,                      # Number of bins for neutron energy
            x_array,                     # Bin edges for neutron energy
            NoOfEnergyBins,              # Number of bins for Etot
            EtotBins_array,              # Bin edges for Etot
            len(Mbins) - 1,              # Number of bins for multiplicity
            Mbins_array                  # Bin edges for multiplicity
        )



        angle_counts_raw = ROOT.TH1F("angle_counts_raw", "angle_counts_raw;Angle (degrees);Counts", 1800, 0, 180)
        angle_counts_selected = ROOT.TH1F("angle_counts_selected", "angle_counts_selected;Angle (degrees);Counts", 1800, 0, 180)
        angle_Esum_2D = ROOT.TH2D("angle_Esum_2D", "angle_Esum_2D;Angle (degrees);ESum (MeV);Counts",
            1800, 0, 180,  #Angle: 0 to 180 degrees  
            400, 0, 20) #Esum: 0 to 15 MeV
        angle_Esum_counts = ROOT.TH3F("angle_Esum_counts", "angle_Esum_counts;Angle (degrees);ESum (MeV);Counts",
            1800, 0, 180,  #Angle: 0 to 180 degrees  
            400, 0, 20, #Esum: 0 to 15 MeV
            1000, 0, 45e6) #Counts: 0 to 45e6
        costheta_Esum_counts = ROOT.TH3F("costheta_Esum_counts", "cos(#theta) vs ESum vs Counts;cos(#theta);ESum (MeV);Counts",
            1800, -1, 1,  #Cos(theta): -1 to 1 
            200, 0, 15, #Esum: 0 to 15 MeV
            1000, 0, 45e6) #Counts: 0 to 45e6
        angle_Esum_counts_cut = ROOT.TH3F("angle_Esum_counts_cut", "Angle vs ESum vs Counts;Angle (degrees);ESum (MeV);Counts",
            1800, 0, 180,  #Angle: 0 to 180 degrees  
            200, 0, 15, #Esum: 0 to 15 MeV
            1000, 0, 45e6) #Counts: 0 to 45e6
        angle_detpairs = ROOT.TH1F("angle_detpairs", "Angle vs. Unique Pairs", 1800, 0, 180)
        clustermult_Esum = ROOT.TH2F("clustermult_Esum", "Cluster Multiplicity vs ESum;Cluster Multiplicity;ESum (MeV)", 30, 0, 30, 200, 0, 15)
        crystalmult_Esum = ROOT.TH2F("crystalmult_Esum", "Crystal Multiplicity vs ESum;Crystal Multiplicity;ESum (MeV)", 30, 0, 30, 200, 0, 15)
        crystalmult_clustermult_Esum = ROOT.TH3F("crystalmult_clustermult_Esum", "Crystal Mult. vs Cluster Mult. vs ESum;Crystal Multiplicity;Cluster Multiplicity;ESum (MeV)",
                                                30, 0, 30,
                                                30, 0, 30,
                                                200, 0, 15)
        #En_Esum_Mcl = ROOT.TH3F("En_Esum_Mcl", "Neutron Energy vs ESum vs Cluster Mult.;En (eV);ESum (MeV);Cluster Mult.",
        #                       200, 0.005, 5e6,
        #                      201, 0, 20,
        #                     21, 0, 20)
        
        angle_event_counts_raw = defaultdict(int)
        angle_event_counts_selected = defaultdict(int)
        angle_sum = defaultdict(float)
        angle_sum_err2 = defaultdict(float)  # store squared errors
        graph_index = 0  # counter for TGraph point index
        detector_pair_set = {}
        for event in eventTree:
            eventType = event.eventType
            ESum = event.ESum
            En = event.En
            crystal_mult = event.crystal_mult
            cluster_mult = event.cluster_mult
            tof_corr = event.tof_corr
            crystal_energies = event.crystal_energies
            eventDetectorIDs = event.eventDetectorIDs
            

            
            if eventType == 'DANCE':
                #En_Esum_Mcl.Fill(En, ESum, cluster_mult)
                #crystalmult_Esum.Fill(crystal_mult, ESum)
                #clustermult_Esum.Fill(cluster_mult, ESum)
                #crystalmult_clustermult_Esum.Fill(crystal_mult, cluster_mult, ESum)
                if crystal_mult == 2 and cluster_mult == 2:
                    det1 = eventDetectorIDs[0]
                    det2 = eventDetectorIDs[1]
                    angle = matrix[det1][det2]
                    #angle = matrix[det1-1][det2-1]
                    
                    angle_event_counts_raw[angle] += 1 # get counts for raw data for each angle
                    angle_Esum_2D.Fill(angle,ESum)
                    angle_Esum_counts.Fill(angle, ESum, 1) # fill hist
                    # Only add point to graph if ESum is within the cut window
                    if lowEsum_cut <= ESum <= highEsum_cut:
                        angle_event_counts_selected[angle] += 1 # get counts for energy selected data for each angle
                        scale_factor = scaling_values.get(angle, 1.0) # get scale_factor
                        weight = scale_factor
                        #print("angle", angle, "scale factor", scale_factor)
                        #weight = scale_factor / normvalue # = 1 / (num of detector pairs * arbitrary normalization value)
                        angle_sum[angle] += weight # gets normalized counts for each angle
                        angle_sum_err2[angle] += weight**2  # for Poisson: sum( (σ=weight)^2 )
                        #norm_count_error = scale_factor / normvalue
                        

                    #costheta_Esum_counts.Fill(ROOT.TMath.Cos(angle * ROOT.TMath.DegToRad()), ESum, 1)

                    detector_pair = tuple(sorted([det1, det2]))
                    if detector_pair not in detector_pair_set:
                        detector_pair_set[detector_pair] = angle

                        #if angle< 21.0:
                            #print(det1, det2, angle)
                        angle_detpairs.Fill(angle)
                    
        with open("sumweights.txt", "w") as f:
            angle_Esum_counts_norm_graph = ROOT.TGraphErrors()
            for i, angle in enumerate(sorted(angle_sum.keys())):
                y = angle_sum[angle]
                ey = math.sqrt(angle_sum_err2[angle])  # √(Σw²)
                print(i, y, ey, file = f)
                angle_Esum_counts_norm_graph.SetPoint(i, angle, y)
                angle_Esum_counts_norm_graph.SetPointError(i, angle_err, ey)
        print(angle_sum)
        #print("sum of weights keys size:" , angle_sum.keys().size())
        # normalize to W(90°) = 1
        '''scale_90 = None
        for i in range(angle_Esum_counts_norm_graph.GetN()):
            x = angle_Esum_counts_norm_graph.GetX()[i]
            if abs(x - 90.0) < 0.5:
                scale_90 = angle_Esum_counts_norm_graph.GetY()[i]
                break

        if scale_90:
            for i in range(angle_Esum_counts_norm_graph.GetN()):
                y = angle_Esum_counts_norm_graph.GetY()[i]
                ey = angle_Esum_counts_norm_graph.GetErrorY(i)
                angle_Esum_counts_norm_graph.SetPoint(i, angle_Esum_counts_norm_graph.GetX()[i], y / scale_90)
                angle_Esum_counts_norm_graph.SetPointError(i, angle_err, ey / scale_90)'''

        print("Writing histograms ...")

        #En_Esum_Mcl.Write()
        #angle_counts.Write()
        #angle_detpairs.Write()
        #crystalmult_clustermult_Esum.Write()
        #crystalmult_Esum.Write()
        #clustermult_Esum.Write()

        # write/project angle_Esum_counts
        #angle_Esum_2D.Write()
        #angle_Esum_counts_yx = angle_Esum_counts.Project3D("yx")
        #angle_Esum_counts_yx.Write()
        #angle_Esum_counts_yx_cut = angle_Esum_counts_yx.Clone("angle_Esum_counts_yx_cut")
        #angle_Esum_counts_yx_projy = angle_Esum_counts_yx.ProjectionY("angle_Esum_counts_yx_projY")
        #angle_Esum_counts_yx_projy.Write()
        #angle_Esum_counts_yx_projx = angle_Esum_counts_yx.ProjectionX("angle_Esum_counts_yx_projX")
        #angle_Esum_counts_yx_projx.Write()
        #EsumAxis = angle_Esum_counts_yx_cut.GetYaxis()
        #lowbin = EsumAxis.FindBin(lowEsum_cut)
        #highbin = EsumAxis.FindBin(highEsum_cut)
        #angle_Esum_counts_yx_projx_cut = angle_Esum_counts_yx_cut.ProjectionX("angle_Esum_counts_yx_projX_cut", lowbin, highbin)
       # angle_Esum_counts_yx_projx_cut.Write()

       
        output_file.cd()
        canvas = ROOT.TCanvas("angle_geometrynormalized", "Angle vs Geometry-Normalized Counts")
        angle_Esum_counts_norm_graph.SetTitle("Angular Correlation;Angle (degrees);Normalized Counts")
        angle_Esum_counts_norm_graph.SetMarkerStyle(20)           # Style 20: solid circles
        angle_Esum_counts_norm_graph.Draw("AP")                   # A = axis, P = points
        canvas.Update()
        canvas.Write("Angle_GeometryNormalizedCounts_drawn")
        angle_Esum_counts_norm_graph.SetName("Angle_GeometryNormalizedCounts_graph")
        angle_Esum_counts_norm_graph.Write("Angle_GeometryNormalizedCounts_graph")
        Npoints = angle_Esum_counts_norm_graph.GetN()
        with open(f'normalized_counts_data_{runnum}.txt', "w") as f:
            for i in range(Npoints):
                x = array('d',[0.])
                y = array('d',[0.])
                angle_Esum_counts_norm_graph.GetPoint(i, x, y)
                ex = angle_Esum_counts_norm_graph.GetErrorX(i)
                ey = angle_Esum_counts_norm_graph.GetErrorY(i)
                #print(f"Point {i}: angle = {x[0]} ± {ex}, count = {y[0]} ± {ey}", file =f)
            
            '''# write/project costheta_Esum_counts
            #costheta_Esum_counts.Write()
            costheta_Esum_counts_yx = costheta_Esum_counts.Project3D("yx")
            #costheta_Esum_counts_yx.Write()
            costheta_Esum_counts_yx_cut = costheta_Esum_counts_yx.Clone("costheta_Esum_counts_yx_cut")
            costheta_Esum_counts_yx_projx = costheta_Esum_counts_yx.ProjectionX("costheta_Esum_counts_yx_projX")
            #costheta_Esum_counts_yx_projx.Write()
            EsumAxis = costheta_Esum_counts_yx_cut.GetYaxis()
            lowbin = EsumAxis.FindBin(lowEsum_cut)
            highbin = EsumAxis.FindBin(highEsum_cut)
            costheta_Esum_counts_yx_projx_cut = costheta_Esum_counts_yx_cut.ProjectionX("costheta_Esum_counts_yx_projX_cut", lowbin, highbin)
            costheta_Esum_counts_yx_projx_cut.Write()'''



    
        print(f"Run {runnum}: Analysis Complete")
output_file.Close()
print("Output file closed successfully")
        

if use_summed_hists == True:
    summedhists_outputfile = f"/data2/lansce/jr2514/AngCorr_Analysis_Output/SummedHists_{startRun}_{endRun}.root"


    #-------------------- sum hists --------------------#

    if sumhists == True:
        print("Summing histograms ...")
        add_hists_command = f"hadd -f {summedhists_outputfile} {' '.join(filestosum)}"
        ROOT.gSystem.Exec(add_hists_command) 
        print(f"Summed histograms written to ~/DANCE/AngularCorrelations/AngCorr_Analysis_Output/SummedHists_{startRun}_{endRun}.root")

#-------------------- do fit --------------------#

if do_fit == True:

    if use_summed_hists == True:
        runRange= f'{startRun}_{endRun}'
        hfile = ROOT.TFile(f"/home/jr2514/DANCE/AngularCorrelations/AngCorr_Analysis_Output/SummedHists_{startRun}_{endRun}.root")
    
    else:
        runRange = startRun
        hfile = ROOT.TFile(f"/home/jr2514/DANCE/DANCE_Analysis/stage1_root/Stage1_Histograms_Run_{startRun}_10ns_CW_0ns_CBT_0ns_DEBT.root")
    print(hfile)
    fit_output_name = f"/home/jr2514/DANCE/AngularCorrelations/AngCorr_Analysis_Output/AngCorrFit_{runRange}.root"
    fit_outputfile = ROOT.TFile(fit_output_name, "RECREATE")
    analysis_outfilename = f"/data2/lansce/jr2514/AngCorr_Analysis_Output/Analysis_Output_{runnum}.root"
    analysis_output = ROOT.TFile(analysis_outfilename, "READ")
    

    angle_geometrynormalized_graph = analysis_output.Get("Angle_GeometryNormalizedCounts_graph")
    if angle_geometrynormalized_graph:
        print("successfully obtained hist from analysis")
    else:
        print("analysis hist is null -- not obtained")
    Npoints = angle_geometrynormalized_graph.GetN()
    #print("Npoints: ", angle_geometrynormalized_graph.GetN())
    for i in range(Npoints):
   # for i in range(Npoints):
        x = array('d',[0.])
        y = array('d',[0.])
        angle_geometrynormalized_graph.GetPoint(i, x, y)
        x_val = x[0]
        y_val = y[0]
        angle_geometrynormalized_graph.SetPoint(i, x_val, y_val)
        #angle_geometrynormalized_graph.SetPoint(i, x_val, y_val / normvalue)
        
        ey = angle_geometrynormalized_graph.GetErrorY(i)
        #print("error_y =", ey)
        #print("normvalue=", normvalue)
        #print("error_y / normvalue =",ey/normvalue)
        ex = angle_geometrynormalized_graph.GetErrorX(i)
        angle_geometrynormalized_graph.SetPointError(i, ex, ey)
        #angle_geometrynormalized_graph.SetPointError(i, ex, ey / normvalue)
    fitFunc_ang = ROOT.TF1("fitFunc_ang", ang_corr_func_ang, fit_limit_low, fit_limit_high,4)
    fitFunc_ang.SetParameters(init_params[0], init_params[1],init_params[2])
    if fix_params == True:
        fitFunc_ang.FixParameter(0, init_params[0])
        fitFunc_ang.FixParameter(1, init_params[1])
        fitFunc_ang.FixParameter(2, init_params[2])



    angle_geometrynormalized_graph.Fit(fitFunc_ang, "R")
    fit_outputfile.cd()
    canv = ROOT.TCanvas("angle_geometrynormalized_fit", "Angular Correlation Fit")
    angle_geometrynormalized_graph.SetTitle("Angular Correlation;Angle (degrees);Counts (geometry-normalized)")
    angle_geometrynormalized_graph.SetMarkerStyle(20)
    angle_geometrynormalized_graph.Draw("AP")
    fitFunc_ang.Draw("same")
    canv.Update()
    canv.Write("Angle_GeometryNormalizedCounts_drawn")
    
    
    # check uncertainty
    for i in range(angle_geometrynormalized_graph.GetN()):
        y = angle_geometrynormalized_graph.GetY()[i]
        ey = angle_geometrynormalized_graph.GetErrorY(i)
        #print(f"angle = {angle_geometrynormalized_graph.GetX()[i]:.1f}, y = {y:.3f}, error = {ey:.3f}, % uncertainty = {100 * ey/y:.2f}%")

    
    #---- SCALE COUNTS using Fit Parameter 3 (A)----#
    
    new_norm = fitFunc_ang.GetParameter(3)
    
    for i in range(angle_geometrynormalized_graph.GetN()):
        x = angle_geometrynormalized_graph.GetX()[i]
        y = angle_geometrynormalized_graph.GetY()[i]
        ex = angle_geometrynormalized_graph.GetErrorX(i)
        ey = angle_geometrynormalized_graph.GetErrorY(i)

        # Divide both y and error by new_norm
        angle_geometrynormalized_graph.SetPoint(i, x, y / new_norm)
        angle_geometrynormalized_graph.SetPointError(i, ex, ey / new_norm)

    # Refit after rescaling
    fitFunc_ang = ROOT.TF1("fitFunc_ang", ang_corr_func_ang, fit_limit_low, fit_limit_high,4)
    fitFunc_ang.SetParameters(init_params[0],init_params[1], init_params[2])  
    angle_geometrynormalized_graph.Fit(fitFunc_ang, "R")
    
    

    # Save new canvas
    fit_outputfile.cd()
    canvas = ROOT.TCanvas("angle_geometrynormalized_fit_scaled", "Scaled Angular Correlation Fit")
    angle_geometrynormalized_graph.SetTitle("Angular Correlation;Angle (degrees);Counts (geometry-normalized and scaled)")
    angle_geometrynormalized_graph.SetMarkerStyle(20)
    angle_geometrynormalized_graph.Draw("AP")
    fitFunc_ang.Draw("same")
    canvas.Update()
    canvas.Write("Angle_GeometryNormalizedCounts_graph_Scaled")

    # Optional: Print final parameters
    a0_final = fitFunc_ang.GetParameter(0)
    a2_final = fitFunc_ang.GetParameter(1)
    a4_final = fitFunc_ang.GetParameter(2)

    print(f"Final (normalized) fit parameters:")
    print(f"A     = {fitFunc_ang.GetParameter(3):.3f}")
    print(f"a0 = {a0_final:.3f}")
    print(f"a2 = {a2_final:.3f}")
    print(f"a4 = {a4_final:.3f}")

    chi2 = fitFunc_ang.GetChisquare()
    ndf = fitFunc_ang.GetNDF()
    chi2_ndf = chi2 / ndf if ndf != 0 else float('nan')

    print(f"Chi2 = {chi2:.3f}")
    print(f"NDF  = {ndf}")
    print(f"Chi2/NDF = {chi2_ndf:.3f}")

    print(f"Fit complete. Results saved to {fit_output_name}.")
    

    '''
    hfile = ROOT.TFile(summedhists_outputfile, "READ")

    fit_output_name = f"~/DANCE/AngularCorrelations/AngCorr_Analysis_Output/AngCorrFit_{startRun}_{endRun}.root"
    # Check if the histogram exists
    angle_counts = hfile.Get("angle_Esum_counts_yx_projX_cut").Clone("angle_counts")
    angle_counts_norm = hfile.Get("angle_Esum_counts_norm_yx_projX_cut").Clone("angle_counts_norm")
    cos_theta = hfile.Get("costheta_Esum_counts_yx_projX_cut").Clone("cos_theta")

    if not angle_counts or not cos_theta:
        print("Error: Histogram 'angle_Esum_counts_yx_projX_cut' or 'costheta_Esum_counts_yx_projX_cut' not found in file.")
    else:
        #angle_counts.GetYaxis().SetRangeUser(1000,2000)
        #cos_theta.GetYaxis().SetRangeUser(1000,2000)
        angle_counts.Write()
        cos_theta.Write()

        # Create a new ROOT file to store fit results
        fit_outputfile = ROOT.TFile(fit_output_name, "RECREATE")
        # create TF1 fit function 
        fitFunc_ang = ROOT.TF1("fitFunc_ang", ang_corr_func_ang, 0, 180,4)
        fitFunc_cos = ROOT.TF1("fitFunc_cos", ang_corr_func_cos, -1, 1,4)
        # set parameters (initialized at top of code line 16)
        fitFunc_ang.SetParameters(init_params[0], init_params[1],init_params[2])
        if fix_params == True:
            fitFunc_ang.FixParameter(0, init_params[0])
            fitFunc_ang.FixParameter(1, init_params[1])
            fitFunc_ang.FixParameter(2, init_params[2])  # Initial parameter values



            

    # draw angle_counts
        angle_counts_canvas = ROOT.TCanvas("angle_counts", "", 800, 600)
        angle_counts.SetTitle("counts vs angle")  # Use #theta for theta in LaTeX
        angle_counts.GetXaxis().SetTitle("angle")  # Use #theta for theta in LaTeX
        angle_counts.GetYaxis().SetTitle("counts")  # Use #theta for theta in LaTeX
        
        if do_fit == True:
            fit_outputfile.cd()
            print("fitting angle_counts")
            angle_counts.Fit(fitFunc_ang, "R")
            angle_counts.Draw()
            fitFunc_ang.Draw("same")
            angle_counts_canvas.Write("AngleX_CountsY_fitted")
        else:
            fit_outputfile.cd()
            angle_counts.Draw()
            angle_counts_canvas.Write("AngleX_CountsY")

    # normalize counts for (angle, counts) hist
        
        angles = []
        angles_err = []
        cos_angle = []
        counts = []
        norm_counts = []
        with open("angle_counts.txt", "r") as file:
            for line in file:
                value = line.split()
                angle = float(value[0])
                #print("angle", angle) #all good
                no_of_pairs = float(value[1])
                #print("no of pairs", no_of_pairs) #all good
                scaling_values[angle] = 1 / no_of_pairs
                #print("scaling value (1/no_of_pairs)", scaling_values[angle]) #all good
        print("scaling value of angle 27.9", scaling_values[27.9])     
        max_count = angle_counts.GetMaximum()
        print("max count", max_count)
        max_count_bin = angle_counts.GetMaximumBin()
        print("max count bin", max_count_bin)
        num_scaledcount_bins = angle_counts.GetNbinsY()
        #angle_counts_norm = ROOT.TH2F("angle_scaledcounts", "Angles vs Scaled Counts;Angle (degrees);Scaled Counts", 1800, 0, 180, num_scaledcount_bins, 0, max_count)
        
        #make TGraph
        angle_counts_norm = ROOT.TGraph()
        angle_counts_norm_err = ROOT.TGraphErrors()
        
        #bin_width = angle_counts.GetXaxis().GetBinWidth(1)
        #x_min = angle_counts.GetXaxis().GetXmin()
        for binx in range(1, angle_counts.GetNbinsX() + 1):
            count = angle_counts.GetBinContent(binx)
            if count != 0.0:
                #print("count / bin content:", count)
                
                
                # go through all the unique angles and see what this angle is closest to
                # Define your computed angle
                angle = angle_counts.GetXaxis().GetBinCenter(binx)
                # Find the closest angle
                closest_angle = min(unique_angles, key=lambda x: abs(x - angle))
                #print("angle bin content:", angle, " closest unique angle:", closest_angle)
                angle = closest_angle
                
                scale_factor = scaling_values[angle]
                print("angle:",angle, "scale_factor:",scale_factor)
                angles.append(angle)
                cos_angle.append(ROOT.TMath.Cos(angle * ROOT.TMath.DegToRad()))
                counts.append(int(count))
                norm_counts.append(count*scale_factor)
        with open("normalized_counts_data.txt", "w") as f:
            for i in range(len(angles)):
                angle_counts_norm.SetPoint(i, angles[i], norm_counts[i]/normvalue)
                angle_counts_norm_err.SetPoint(i, angles[i], norm_counts[i]/normvalue)
                angle_counts_norm_err.SetPointError(i, angle_err, 0.0)
                print(f"angle: {angles[i]}, count: {counts[i]}, normalized count: {norm_counts[i]}", file=f)




        fitFunc_ang.SetParameters(init_params[0], init_params[1],init_params[2])
        if fix_params == True:
            fitFunc_ang.FixParameter(0, init_params[0])
            fitFunc_ang.FixParameter(1, init_params[1])
            fitFunc_ang.FixParameter(2, init_params[2])



    # draw angle_counts_norm
        angle_counts_norm_canvas = ROOT.TCanvas("angle_counts_norm","angle_counts_norm")
        angle_counts_norm.SetTitle("x: angle(#theta) y: normalized counts")  # Use #theta for theta in LaTeX
        angle_counts_norm.GetXaxis().SetTitle("angle #theta")  # Use #theta for theta in LaTeX
        angle_counts_norm.GetYaxis().SetTitle("normalized counts")  # Use #theta for theta in LaTeX
        angle_counts_norm.SetMarkerStyle(20)

        angle_counts_norm_err_canvas = ROOT.TCanvas("angle_counts_norm_err","angle_counts_norm_err")
        angle_counts_norm_err.SetTitle("x: angle(#theta) y: normalized counts")  # Use #theta for theta in LaTeX
        angle_counts_norm_err.GetXaxis().SetTitle("angle #theta")  # Use #theta for theta in LaTeX
        angle_counts_norm_err.GetYaxis().SetTitle("normalized counts")  # Use #theta for theta in LaTeX
        angle_counts_norm_err.SetMarkerStyle(20)
        if do_fit == True:
            print("fitting angle_counts_norm")
            angle_counts_norm.Fit(fitFunc_ang, "R")
            fit_outputfile.cd()
            angle_counts_norm.Draw("AP")
            fitFunc_ang.Draw("same")
            angle_counts_norm_canvas.Update()
            angle_counts_norm_canvas.Write("AngleX_NormCountsY_fitted")

            print("fitting angle_counts_norm_err")
            #angle_counts_norm_err.Fit(fitFunc_ang,"R")
            fit_outputfile.cd()
            angle_counts_norm_err.Draw("AP")
            fitFunc_ang.Draw("same")
            angle_counts_norm_err_canvas.Update()
            angle_counts_norm_err_canvas.Write("AngleX_err_NormCountsY_fitted")
        else:
            fit_outputfile.cd()
            angle_counts_norm.Draw("AP")
            angle_counts_norm_canvas.Write("AngleX_NormCountsY")
            angle_counts_norm_err.Draw("AP")
            angle_counts_norm_err_canvas.Write("AngleX_err_NormCountsY")


        
        fitFunc_cos.SetParameters(init_params[0], init_params[1],init_params[2])
        if fix_params == True:
            fitFunc_cos.FixParameter(0, init_params[0])
            fitFunc_cos.FixParameter(1, init_params[1])
            fitFunc_cos.FixParameter(2, init_params[2])

        if do_fit == True:
            print("fitting cos_theta")
            cos_theta.Fit(fitFunc_cos, "R")


    # draw cos_theta_counts
        canvas_costheta_counts = ROOT.TCanvas("costheta_counts", "Counts vs cos(#theta)", 800, 600)
        cos_theta.SetTitle("counts vs cos(#theta)")  # Use #theta for theta in LaTeX
        cos_theta.GetXaxis().SetTitle("cos(#theta)")  # Use #theta for theta in LaTeX
        cos_theta.GetYaxis().SetTitle("counts")  # Use #theta for theta in LaTeX
        cos_theta.Draw("E")
        if do_fit == True:
            fitFunc_cos.Draw("same")
            fit_outputfile.cd()
            canvas_costheta_counts.Write("Cos(theta)X_CountsY_fitted")
        else:
            fit_outputfile.cd()
            canvas_costheta_counts.Write("Cos(theta)X_CountsY")


    # normalize counts for ( cos(theta), counts ) hist 
        num_scaledcount_bins = cos_theta.GetNbinsY() + 1
        max_count = cos_theta.GetMaximum()
        cos_theta_counts_norm = ROOT.TH2F("cos_theta_counts_norm", "x: cos(#theta) y: Scaled Counts;cos(#theta);Normalized Counts", 1800, -1.0, 1.0, num_scaledcount_bins, 0, max_count)
        for i in range(len(cos_angle)):
            cos_theta_counts_norm.Fill(cos_angle[i], norm_counts[i])

    # draw cos_theta_counts_norm
        cos_theta_counts_norm_canvas = ROOT.TCanvas("cos_theta_counts_norm","cos_theta_counts_norm",800, 600)
        cos_theta_counts_norm.SetTitle("x: cos(#theta) y: normalized counts")  # Use #theta for theta in LaTeX
        cos_theta_counts_norm.GetXaxis().SetTitle("cos(#theta)")  # Use #theta for theta in LaTeX
        cos_theta_counts_norm.GetYaxis().SetTitle("normalized counts")  # Use #theta for theta in LaTeX
        cos_theta_counts_norm.SetMarkerStyle(20)
        if do_fit == True:
            fit_outputfile.cd()
            cos_theta_counts_norm.Draw()
            fitFunc_cos.Draw("same")
            cos_theta_counts_norm_canvas.Write("Cos(theta)X_NormCountsY_fitted")
        else:
            fit_outputfile.cd()
            cos_theta_counts_norm.Draw()
            cos_theta_counts_norm_canvas.Write("Cos(theta)X_NormCountsY")

        # Create histograms for w(theta) and g(theta)
        nbins = angle_counts.GetNbinsX()  # Use the same binning as the data
        cos_bins = np.linspace(-1, 1, nbins + 1)  # Define bins for cos(theta)

        h_w_costheta = ROOT.TH1D("h_w_costheta", "w(#theta) vs cos(#theta)", nbins, -1, 1)
        h_g_theta = ROOT.TH1D("h_g_costheta", "g(#theta) {= w(#theta) / w(90)} vs #theta", nbins, 0, 180)

        # Fit the angle_counts histogram
        angle_counts.Fit(fitFunc_ang, "R")

        # Extract the fit parameters
        p0 = fitFunc_ang.GetParameter(0)  # p0 (constant)
        p1 = fitFunc_ang.GetParameter(1)  # p1 (coefficient of cos^2)
        p2 = fitFunc_ang.GetParameter(2)  # p2 (coefficient of cos^4)

        w_90 = p0  # Reference value w(90)

        # Loop over the actual bin centers of angle_counts
        for i in range(1, nbins + 1):
            theta = angle_counts_norm.GetBinCenter(i)  # Get actual theta from histogram
            cos_theta = ROOT.TMath.Cos(theta * ROOT.TMath.DegToRad())  # Convert to cos(theta)

            # Calculate W(theta) using the fitted parameters
            w_theta = p0 + (p1 * cos_theta**2) + (p2 * cos_theta**4)

            # Fill the histograms
            h_w_costheta.Fill(cos_theta, w_theta)
            h_g_theta.Fill(theta, w_theta / w_90)

        # Create and save the w(theta) vs cos(theta) plot
        canvas_w = ROOT.TCanvas("w(theta)_cos(theta)", "w(#theta) vs cos(#theta)", 800, 600)
        h_w_costheta.SetTitle("w(#theta) vs cos(#theta)")
        h_w_costheta.GetXaxis().SetTitle("cos(#theta)")
        h_w_costheta.GetYaxis().SetTitle("w(#theta)")
        h_w_costheta.Draw("E")
        fit_outputfile.cd()
        canvas_w.Write("w(theta)X_cos(theta)Y")

        # Create and save the g(theta) vs theta plot
        canvas_g = ROOT.TCanvas("g(theta)_theta", "g(#theta) {= w(#theta)/w(90)} vs #theta ", 800, 600)
        h_g_theta.SetTitle("g(#theta) vs #theta")
        h_g_theta.GetXaxis().SetTitle("#theta")
        h_g_theta.GetYaxis().SetTitle("g(#theta)")
        h_g_theta.Draw("E")
        fit_outputfile.cd()
        canvas_g.Write("g(theta)X_thetaY")'''



        # Save everything into the new ROOT file





