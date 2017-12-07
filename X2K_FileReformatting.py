directory = "/Users/Schilder/Desktop/x2k/"

# # Delete .DS_Store files
import os
print("Deleting.......")
for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if file.endswith(".DS_Store"):
            print(os.path.join(root, file))
            os.remove( os.path.join(root, file) )

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


# Subset testgmt file for overfitting test
with open(directory + "Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_FILTERED.txt") as testgmt,\
    open(directory + "Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.txt","w") as subset1,\
    open(directory + "Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET2.txt","w") as subset2:
    count = 1
    for line in testgmt:
        if count <= 130:
            subset1.write(line)
        if count > 130:
            subset2.write(line)
        count += 1


