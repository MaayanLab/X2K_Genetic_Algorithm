#RUIIUSNUFSDF TESTrer

directory = "/Users/Schilder/Desktop/x2k/"

# Filter GEO kinase perturbation data
## 1. By simply stripping the 1.0 after each gene
import re
with open(directory+"Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_up.gmt") as GEO_up,\
        open(directory+"Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_up_noVals.txt","w") as noVals_up,\
        open(directory+"Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_down.gmt") as GEO_dn,\
        open(directory+"Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_dn_noVals.txt","w") as noVals_dn:
    for line in GEO_up:
        splitLine = re.split(r'\t+', line)
        experiment = splitLine[0]
        genes = splitLine[1:-1]
        editedGenes=[]
        for gene in genes:
            editedGenes.append(gene.strip(",1.0"))
        noVals_up.write(experiment+"\t\t"+"\t".join(editedGenes)+"\n")
    for line in GEO_dn:
        splitLine = re.split(r'\t+', line)
        experiment = splitLine[0]
        genes = splitLine[1:-1]
        editedGenes=[]
        for gene in genes:
            editedGenes.append(gene.strip(",1.0"))
        noVals_dn.write(experiment+"\t\t"+"\t".join(editedGenes)+"\n")



## 1. By matching perturbation types with up/dn files
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





directory = "/Users/Schilder/Desktop/x2k/"
# Subset testgmt file for overfitting test
def split_testgmt(testPercent, input, output1, output2, writeType="w", up_dn="dn"):
    with open(directory + "Validation/Perturbation_Data/"+input) as testgmt,\
        open(directory + "Validation/Perturbation_Data/"+output1, writeType) as subset1,\
        open(directory + "Validation/Perturbation_Data/"+output2, writeType) as subset2:
        testgmtLines = testgmt.readlines()
        division = round(len(testgmtLines)*testPercent/100)
        count = 1
        for line in testgmtLines:
            if count <= division:
                newLine = line.split("\t")
                newLine[0] = newLine[0]+"_"+up_dn
                subset1.write("\t".join(newLine))
            if count > division:
                newLine = line.split("\t")
                newLine[0] = newLine[0]+"_"+up_dn
                subset2.write("\t".join(newLine))
            count += 1


# Split GEO gene pert data and then combine up/dn files
split_testgmt(testPercent=80, \
              input="GEO/Kinase_Perturbations_from_GEO_up_noVals.txt", \
              output1="GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt", \
              output2="GEO/Kinase_Perturbations_from_GEO_SUBSET.20per.txt", \
              up_dn="up")
split_testgmt(testPercent=80, \
              input="GEO/Kinase_Perturbations_from_GEO_up_noVals.txt", \
              output1="GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt", \
              output2="GEO/Kinase_Perturbations_from_GEO_SUBSET.20per.txt", \
              up_dn="dn", writeType="a")


# Split LINCS L1000 chem pert data and then combine up/dn files
split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_dn_DRH.kinaseInihibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt",\
              up_dn="up")
split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_up_DRH.kinaseInihibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt",\
              up_dn="dn", writeType="a")


directory = "/Users/Schilder/Desktop/x2k/"

# Mix up and down, then write to new file, then get top n lines (to limit computational load on GA)
def reduced_split_testgmt(totalExperiments, testPercent, input, output1, output2, writeType="w", up_dn="dn"):
    with open(directory + "Validation/Perturbation_Data/"+input) as testgmt,\
        open(directory + "Validation/Perturbation_Data/"+output1, writeType) as subset1,\
        open(directory + "Validation/Perturbation_Data/"+output2, writeType) as subset2:
        testgmtLines = testgmt.readlines()
        #division = round(len(testgmtLines)*testPercent/100)
        division = round((totalExperiments*.80)/2)
        count = 1

        for line in testgmtLines:
            if count <= division:
                newLine = line.split("\t")
                newLine[0] = newLine[0]+"_"+up_dn
                subset1.write("\t".join(newLine))
            if count > division:
                newLine = line.split("\t")
                newLine[0] = newLine[0]+"_"+up_dn
                subset2.write("\t".join(newLine))
            count += 1

#directory = "/home/maayanlab/PycharmProjects/X2K_Genetic_Algorithm/" # Change to your directory
directory = "/Users/schilder/Desktop/X2K_Genetic_Algorithm/"

# LINCS L1000 + DRH
reduced_split_testgmt(totalExperiments=570, testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_dn_DRH.kinaseInihibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt",\
              up_dn="dn")
reduced_split_testgmt(totalExperiments=570, testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_up_DRH.kinaseInihibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt",\
              up_dn="up", writeType="a")




# LINCS L1000 + KINOMEscan
reduced_split_testgmt(totalExperiments=570, testPercent=80,\
              input="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.DN_filtered.txt",\
              output1="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.SUBSET1",\
              output2="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.SUBSET1",\
              up_dn="dn")
reduced_split_testgmt(totalExperiments=570, testPercent=80,\
              input="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.UP_filtered.txt",\
              output1="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.SUBSET1.txt",\
              output2="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.SUBSET2.txt",\
              up_dn="up", writeType="a")