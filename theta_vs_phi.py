import ROOT
import math
import sys
import numpy as np

detectors = []
theta = []
phi = []
radius = []
shapes = []
red_theta = []
red_phi= []
blue_theta = []
blue_phi = []
green_theta = []
green_phi = []
yellow_theta = []
yellow_phi = []

file_path_angles = '/home/jr2514/DANCE/AngularDistributions/angles.txt'
file_path_shapes = '/home/jr2514/DANCE/AngularDistributions/crystalshapes.txt'

with open(file_path_shapes, 'r') as file:
    for line in file:
        terms = line.split()
        print(terms)
        shapes.append(str(terms[1]))
        detectors.append(int(terms[0]))

output_filename = "theta_vs_phi.root"
output_file = ROOT.TFile(output_filename, "RECREATE") 


with open(file_path_angles, 'r') as file:
    for line in file:
                # Step 3: Split the line into terms
        terms = line.split()  
        theta.append(float(terms[2]))
        phi.append(float(terms[3]))
        radius.append(float(terms[5]))

for i in range(0, len(theta)):
    if shapes[i] == 'b':
        blue_theta.append(theta[i])
        blue_phi.append(phi[i])
    if shapes[i] == 'g':
        green_theta.append(theta[i])
        green_phi.append(phi[i])
    if shapes[i] == 'r':
        red_theta.append(theta[i])
        red_phi.append(phi[i])
    if shapes[i] == 'y':
        yellow_theta.append(theta[i])
        yellow_phi.append(phi[i])


assert len(theta) == len(phi), "Theta and Phi lists must be of the same length."






# Convert lists to numpy arrays
red_theta_arr = np.array(red_theta)
red_phi_arr = np.array(red_phi)
blue_theta_arr = np.array(blue_theta)
blue_phi_arr = np.array(blue_phi)
green_theta_arr = np.array(green_theta)
green_phi_arr = np.array(green_phi)
yellow_theta_arr = np.array(yellow_theta)
yellow_phi_arr = np.array(yellow_phi)
mg = ROOT.TMultiGraph()

g1 = ROOT.TGraph(len(blue_theta_arr), blue_theta_arr, blue_phi_arr)
g1.SetMarkerColor(ROOT.kBlue)
g1.SetMarkerStyle(107)
mg.Add(g1)

g2 = ROOT.TGraph(len(red_theta_arr), red_theta_arr, red_phi_arr)
g2.SetMarkerColor(ROOT.kRed)
g2.SetMarkerStyle(107)
mg.Add(g2)

g3 = ROOT.TGraph(len(green_theta_arr), green_theta_arr, green_phi_arr)
g3.SetMarkerColor(ROOT.kGreen)
g3.SetMarkerStyle(107)
mg.Add(g3)

g4 = ROOT.TGraph(len(yellow_theta_arr), yellow_theta_arr, yellow_phi_arr)
g4.SetMarkerStyle(107)
g4.SetMarkerColor(ROOT.kYellow+2)

mg.Add(g4)
#graph.SetMarkerSize(4.0)
# Set graph title and axis labels
mg.SetTitle("Theta vs Phi; Theta (degrees); Phi (degrees)")



# Draw the graph
c1 = ROOT.TCanvas("c1", "ThetavsPhi")
legend = ROOT.TLegend(0.89, 0.6, 0.98, 0.8)  # Position of legend (x1, y1, x2, y2)
legend.AddEntry(g3, "A", "p")
legend.AddEntry(g2, "B", "p")     # "p" indicates point style
legend.AddEntry(g1, "C", "p")
legend.AddEntry(g4, "D", "p")
#legend.SetFillStyle(3001)
mg.Draw("AP")  # A: axis, P: points
legend.Draw()
c1.Write()

# Update and show the canvas

