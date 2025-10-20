import ROOT
import math
from math import acos, cos, sin, radians, degrees, trunc
import numpy as np

def GetAngle():
    cos_angle_btwn =  cos(radians(theta[i]))*cos(radians(theta[j])) + sin(radians(theta[i]))*sin(radians(theta[j]))*cos(radians(abs(phi[i]-phi[j])))
    cos_angle_btwn = max(-1.0, min(1.0, cos_angle_btwn))
    angle_btwn_deg = round(math.degrees(math.acos(cos_angle_btwn)),1)
    return angle_btwn_deg

theta = []
phi = []
angles = []
unique_angles = []
detector_combo1 = []
detector_combo2 = []
combos = []
detectors = []



theta = (90. , 99.6795 , 105.786 , 90 , 74.2137 , 80.3205 , 90. , 108.001 , 117.284 , 120.001 , 106.458 , 90 , 73.5423 , 59.9987 , 62.716 , 71.9993 , 97.41 , 115.379 , 125.971 , 131.842 , 133.908 , 119.471 , 106.458 , 90 , 73.5423 , 60.529 , 46.0923 , 48.1581 , 54.0294 , 64.6209 , 82.59 , 89.9979 , 105.097 , 121.719 , 134.479 , 144.002 , 149.499 , 148.287 , 133.907 , 120.003 , 105.786 , 90 , 74.2137 , 60.0005 , 46.0934 , 31.7134 , 30.501 , 35.9954 , 45.5211 , 58.281 , 74.9026 , 97.4085 , 115.378 , 134.478 , 150.535 , 161.868 , 164.908 , 149.496 , 131.839 , 117.284 , 99.6795 , 80.3205 , 62.716 , 48.1612 , 30.5037 , 15.0917 , 18.1315 , 29.4654 , 45.5223 , 64.6224 , 82.5915 , 90 , 107.998 , 125.968 , 143.998 , 161.864 , 180 , 161.864 , 144.002 , 125.968 , 108.002 , 90 , 72.0019 , 54.0322 , 36.0019 , 18.1361 , 0 , 18.1361 , 35.9981 , 54.0322 , 71.9981 , 99.6795 , 117.284 , 131.839 , 149.496 , 164.908 , 161.868 , 150.535 , 134.478 , 115.378 , 97.4085 , 82.5915 , 64.6224 , 45.5223 , 29.4654 , 18.1315 , 15.0918 , 30.5037 , 48.1612 , 62.716 , 80.3205 , 90 , 105.786 , 120 , 133.907 , 148.287 , 149.499 , 144.005 , 134.479 , 121.719 , 105.097 , 90.0021 , 74.9026 , 58.281 , 45.5211 , 35.9976 , 30.501 , 31.7134 , 46.0935 , 59.9966 , 74.2137 , 90 , 106.458 , 119.471 , 133.908 , 131.842 , 125.971 , 115.379 , 97.41 , 82.59 , 64.6209 , 54.0294 , 48.1581 , 46.0923 , 60.529 , 73.5423 , 90 , 106.458 , 120.001 , 117.284 , 108.001 , 90 , 71.9993 , 62.716 , 59.9987 , 73.5423 , 90 , 105.786 , 99.6795 , 80.3205 , 74.2137 , 90)
phi = (0. , 346.422 , 5.27057 , 16.6216 , 5.27057 , 346.422 , 331.184 , 333.434 , 350.352 , 10.8129 , 23.9913 , 31.7189 , 23.9913 , 10.8129 , 350.352 , 333.434 , 318.313 , 319.236 , 336.211 , 353.747 , 18.2257 , 31.7205 , 39.4487 , 46.8184 , 39.4487 , 31.7205 , 18.2257 , 353.747 , 336.211 , 319.236 , 318.313 , 301.713 , 301.713 , 301.712 , 315.341 , 333.428 , 359.312 , 31.7231 , 45.2171 , 52.6262 , 58.1694 , 63.44 , 58.1694 , 52.6279 , 45.2171 , 31.7231 , 359.312 , 333.434 , 315.341 , 301.712 , 301.713 , 285.113 , 284.188 , 288.081 , 301.709 , 326.181 , 31.7284 , 64.1331 , 69.6945 , 73.0877 , 77.0176 , 77.0176 , 73.0877 , 69.6945 , 64.1331 , 31.7284 , 326.181 , 301.709 , 288.081 , 284.188 , 285.113 , 272.256 , 270 , 267.213 , 270 , 277.264 , 0 , 97.2639 , 90 , 87.2127 , 90 , 92.2556 , 90 , 87.2126 , 90 , 97.2639 , 0 , 277.264 , 270 , 267.213 , 270 , 257.018 , 253.088 , 249.694 , 244.133 , 211.728 , 146.181 , 121.709 , 108.081 , 104.188 , 105.113 , 105.113 , 104.188 , 108.081 , 121.709 , 146.181 , 211.728 , 244.133 , 249.694 , 253.088 , 257.018 , 243.44 , 238.169 , 232.628 , 225.217 , 211.723 , 179.312 , 153.434 , 135.341 , 121.712 , 121.713 , 121.713 , 121.713 , 121.712 , 135.341 , 153.428 , 179.312 , 211.723 , 225.217 , 232.626 , 238.169 , 226.818 , 219.449 , 211.72 , 198.226 , 173.747 , 156.211 , 139.236 , 138.313 , 138.313 , 139.236 , 156.211 , 173.747 , 198.226 , 211.72 , 219.449 , 211.719 , 203.991 , 190.813 , 170.352 , 153.434 , 151.184 , 153.434 , 170.352 , 190.813 , 203.991 , 196.622 , 185.271 , 166.422 , 166.422 , 185.271 , 180.)

#print(len(theta)) #160
for i in range(0,len(theta)):
    for j in range(0,len(theta)):
        if i == 76 and j == 76:
            print("coordinates of detector 76:",theta[i], phi[i])
        if i == 86 and j == 86:
            print("coordinates of detector 86:",theta[i], phi[i])
        #print("det",i, " and det",j, "  ", degrees(acos(cos(radians(theta[i]))*cos(radians(theta[j])) + sin(radians(theta[i]))*sin(radians(theta[j]))*cos(abs(radians(phi[i]) - radians(phi[j]))))))
        #print("det",i, " and det",j, "  ", cos(radians(theta[i]))*cos(radians(theta[j])) + sin(radians(theta[i]))*sin(radians(theta[j]))*cos(radians(abs(phi[i]-phi[j]))))
        angle_btwn = GetAngle() # cosine of the angle
        
        '''if angle_cos == 1.0:
            angle_cos == 1.0
        if angle_cos < -0.99999:
            angle_cos == -1.0'''
        combo1 = f"{i}"
        combo2 = f"{j}"
        detector_combo1.append(combo1)
        detector_combo2.append(combo2)
        angles.append(angle_btwn)

matrix = np.zeros((162, 162))
#print(len(angles))
for ang in angles:
    if ang not in unique_angles:
        unique_angles.append(ang)   
data = []
#print(unique_angles)
#print(len(unique_angles))

for i in range(0,len(unique_angles)):
    for j in range(0,len(angles)):
        if unique_angles[i] == angles[j] and float(unique_angles[i]) != 0.0:
            entry = [unique_angles[i], detector_combo1[j], detector_combo2[j]]
            data.append(entry)
            matrix[int(detector_combo1[j])][int(detector_combo2[j])] = float(unique_angles[i])  
            #print(unique_angles[i], detector_combo1[j], detector_combo2[j], file=file)
#print(data) 
           

# Iterate over the data
    
angle_counts = {}
detector_pair_set = []
unique_data = []


for entry in data:
    angle, detector1, detector2 = entry
    # include all detectors (76,86) or else normalization won't work properly
    detector_pair = tuple(sorted([detector1, detector2]))
    if detector_pair not in detector_pair_set:
        detector_pair_set.append(detector_pair)
        unique_data.append(entry)
        # Count the number of unique detector pairs for each angle
        angle_counts[angle] = angle_counts.get(angle, 0) + 1

print(detector_pair_set)

with open("angle_counts.txt", "w") as file:
    for angle, counts in sorted(angle_counts.items(), key=lambda x: x[0]):
        file.write(f"{angle} {counts}\n")
   
        
         

# has exactly 160 choose 2 entries
with open("angle_combinations.txt", "w") as file:
    for line in unique_data:
        angle, detector1, detector2 = line
        file.write(f"{angle} {detector1} {detector2}\n")
          

with open("unique_angles.txt", "w") as file:
    for y in unique_angles:
        file.write(f"{y}\n") 


submatrix = matrix[:, :]
with open("detector_angles_matrix.txt","w") as file:
    for row in submatrix:
        row_string = ' '.join(f'{val:.2f}' for val in row)  # Adjust format if needed
        file.write(row_string + '\n')

file_path = 'detector_angles_matrix.txt'
matrix = np.loadtxt(file_path)
# Check symmetry
symmetric = True
with open("nearest_neighbors.txt","w") as file:
    for i in range(1,161):
        for j in range(1,161):
            if matrix[i][j] != matrix[j][i]:
                print(f"Not symmetric at ({i}, {j}): {matrix[i][j]} != {matrix[j][i]}")
                symmetric = False
            if matrix[i][j] < 20.0 and matrix[i][j] > 0.0:
                print(i+1, j+1, matrix[i][j], file=file)


if symmetric:
    print("The matrix is symmetric.")
else:
    print("The matrix is NOT symmetric.")
            
        