# # Delete .DS_Store files
import os
print("Deleting.......")
for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if file.endswith(".DS_Store"):
            print(os.path.join(root, file))
            os.remove( os.path.join(root, file) )



directory = "/Users/Schilder/Desktop/x2k/"

# Filter kinase perturbation data
import re
with open(directory + "Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_COMBINED.txt") as kinase_pert_combo, \
        open(directory + "Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_FILTERED.txt","w") as filteredFile:
    for line in kinase_pert_combo:
        splitLine = re.split(r'\t+', line)
        experiment = splitLine[0]
        condition = splitLine[0].split("_")[1]
        direction = splitLine[0][-3:]
        if condition in ("druginhibition", "knockdown", "knockout", "defectivemutant") and direction == "-dn":
            filteredFile.write(line)
            print("upper")
        elif condition in ("drugactivation", "overexpression") and direction == "-up":
            filteredFile.write(line)
            print("downer")
        elif "mutant" or "activemutant":
            print("Not sure?...")
    kinase_pert_combo.close()
    filteredFile.close()