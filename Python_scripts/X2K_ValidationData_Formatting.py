# ____________________________________________________________
# ______________X2K Validation Data Reformatting______________
# ____________________________________________________________

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


# Combine L1000_kinasePerturbation UpDn files
def reformatL1000_kinasePert():
    with open("../X2K_Databases/Validation/Perturbation_Data/LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_up.txt") as Lup,\
        open("../X2K_Databases/Validation/Perturbation_Data/LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_down.txt") as Ldn, \
        open("../X2K_Databases/Validation/Perturbation_Data/LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations.txt","w") as combo:
        up = Lup.readlines()
        dn = Ldn.readlines()
        def format(lines, UpDn):
            for line in lines:
                lsplit = line.split("\t")
                gene = lsplit[0].split("_")[0]
                name = "_".join(lsplit[0].split("_")[1:]+[UpDn])
                DEGs=lsplit[2:]
                newDEGS=[]
                for d in DEGs:
                    newDEGS.append(d.replace(",1.0",""))
                combo.write(gene+"\t"+name+"\t"+"\t".join(newDEGS))
        format(up, "up")
        format(dn, "dn")






directory = "/Users/Schilder/Desktop/X2K_Genetic_Algorithm/"
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
              input="GEO/Kinase_Perturbations_from_GEO_dn_noVals.txt", \
              output1="GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt", \
              output2="GEO/Kinase_Perturbations_from_GEO_SUBSET.20per.txt", \
              up_dn="dn", writeType="a")


# Split LINCS L1000 chem pert data and then combine up/dn files
split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_up_DRH.kinaseInihibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET2.txt",\
              up_dn="up")
split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_dn_DRH.kinaseInhibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET2.txt",\
              up_dn="dn", writeType="a")


directory = "/Users/Schilder/Desktop/X2K_Genetic_Algorithm/"

# Mix up and down, then write to new file, then get top n lines (to limit computational load on GA)
def reduced_split_testgmt(testPercent, input, output1, output2, writeType="w", up_dn="dn"):
    with open(directory + "Validation/Perturbation_Data/"+input) as testgmt,\
        open(directory + "Validation/Perturbation_Data/"+output1, writeType) as subset1,\
        open(directory + "Validation/Perturbation_Data/"+output2, writeType) as subset2:
        testgmtLines = testgmt.readlines()
        ## Can put a cap on dataset to limit computational load
        division = round(len(testgmtLines)*testPercent/100)
        #division = round((570*.80)/2)
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

## LINCS L1000 + DRH
reduced_split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_up_DRH.kinaseInhibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET2.txt",\
              up_dn="up")
reduced_split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_dn_DRH.kinaseInhibitors.txt",\
              output1="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInhibitors_SUBSET2.txt",\
              up_dn="dn", writeType="a")

## LINCS L1000 + KinomeScan
reduced_split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.UP_filtered.txt",\
              output1="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET2.txt",\
              up_dn="up")
reduced_split_testgmt(testPercent=80,\
              input="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan.DN_filtered.txt",\
              output1="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET1.txt",\
              output2="LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET2.txt",\
              up_dn="dn", writeType="a")
## LINCS L1000 Kinase Perturbations ****
reduced_split_testgmt(testPercent=80,\
              input="LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_down.txt",\
              output1="LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_SUBSET1.80per.txt",\
              output2="LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_SUBSET2.20per.txt",\
              up_dn="up")
reduced_split_testgmt(testPercent=80,\
              input="LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_up.txt",\
              output1="LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_SUBSET1.80per.txt",\
              output2="LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_SUBSET2.20per.txt",\
              up_dn="dn", writeType="a")


Dir = "X2K_Genetic_Algorithm/Validation/Perturbation_Data/" # Change to your directory

# Combine datasets:
with open(Dir+"LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET1.80per.txt") as L1000_sub1,\
        open(Dir+"LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.20per.txt") as L1000_sub2,\
        open(Dir+"GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt") as GEO_sub1,\
        open(Dir+"GEO/Kinase_Perturbations_from_GEO_SUBSET2.20per.txt") as GEO_sub2,\
        open(Dir+"Combined/GEO-KinasePert_L1000-DRH_SUBSET1-80per.txt", "a") as sub1,\
        open(Dir+"Combined/GEO-KinasePert_L1000-DRH_SUBSET2-20per.txt", "a") as sub2:
    for line in GEO_sub1:
        lineSp = line.split("\t")
        name = lineSp[0].split("_")
        newname = "GEO-KinasePert_" + "_".join([ "-".join(name[1:]), "["+name[0]+"]"])
        newLine = newname+"\t"+ "\t".join(lineSp[1:])
        sub1.write(newLine)
    for line in L1000_sub1:
        lineSp = line.split("\t")
        name = lineSp[0].split("_")
        newName = "L1000-DrugRepurposingHub_" + "-".join( name[0:3]+[name[-1]]) + "_"+ name[3]
        newLine = newName+"\t"+"\t".join(lineSp[1:])
        sub1.write(newLine)
        
    for line in GEO_sub2:
        lineSp = line.split("\t")
        name = lineSp[0].split("_")
        newname = "GEO-KinasePert" + "_".join(["-".join(name[1:]), "[" + name[0] + "]"])
        newLine = newname + "\t" + "\t".join(lineSp[1:])
        sub2.write(newLine)
    for line in L1000_sub2:
        lineSp = line.split("\t")
        name = lineSp[0].split("_")
        newName = "L1000-DrugRepurposingHub_" + "-".join(name[0:3] + [name[-1]]) + "_" + name[3]
        newLine = newName+"\t"+"\t".join(lineSp[1:])
        sub2.write(newLine)


## ------- Scramble DGEs in in testGMT to make a different randomized baseline: ------- ##
with open("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt") as testGMT,\
        open("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per-SCRAMBLED.txt","w") as scrambled:
    GMTlines = testGMT.readlines()
    DEGlist=[]; lenList=[]
    for line in GMTlines:
        DEGs = line.split("\t")[2:]
        DEGs[-1] = DEGs[-1].strip("\n")
        for deg in DEGs:
            DEGlist.append( deg )
        lenList.append( len(DEGs) )

    from random import shuffle
    shuffle(DEGlist)

    import numpy as np
    for i,line in enumerate(GMTlines):
        print(line)
        spltLine = line.split("\t")
        newDEGs = np.random.choice(DEGlist, size=lenList[i], replace=False)
        scrambled.write(spltLine[0]+"\t"+spltLine[1]+"\t"+"\t".join(newDEGs)+"\n")


## ----- Check the frequency of each Kinase in the testGMT ----- ##
with open("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt") as GEOgmt:
    import pandas as pd
    def getKinaseGroupAndFamily(kinaseList):
        geneSyn = pd.read_table("../X2K_Databases/General_Resources/Moshe_mapping/mappingFile_2017.txt", header=None,
                                names=["Gene", "GeneSyn"])
        kinaseGroups = pd.read_csv("../X2K_Databases/General_Resources/Kinase.com/Kinase_Groups_&_Families.csv")
        kgroupDict = dict(zip(kinaseGroups.Name, zip(kinaseGroups.Group, kinaseGroups.Family)))
        for k in kinaseList:
            print(k)
            syns = [k]
            if k in list(geneSyn.Gene):
                newSyns = list(geneSyn[geneSyn.Gene == k].GeneSyn.unique())
                for s in newSyns:
                    syns.append(s)
            overlap = list(set(syns).intersection(set(kgroupDict.keys())))
            if len(overlap) > 0:
                results = kgroupDict[overlap[0]]
            else:
                results = ["NA", "NA"]
        return results


    import collections
    gmt = GEOgmt.readlines()
    kinaseList=[]
    for line in gmt:
        kinaseList.append( line.split('\t')[0].split("_")[0] )
    # Kinase-level
    x = collections.Counter(kinaseList)
    x.most_common()
    # Get Group, Family
    GROUPS=[]; FAMILIES=[]
    for ki in kinaseList:
        group, family = getKinaseGroupAndFamily([ki])
        GROUPS.append(group); FAMILIES.append(family)
    # Group-level
    collections.Counter(GROUPS).most_common()
    # Family-level
    collections.Counter(FAMILIES).most_common()


## ---- Filter GMT to only include one perturbation of each kinase ---- ##
with open("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt") as GEOgmt,\
        open("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per-UniqueKinases.txt","w") as filtGEOgmt:
    GEOlines = GEOgmt.readlines()
    kinaseRecorder=[]
    for line in GEOlines:
        kinase = line.split("")


directory = "/Users/Schilder/Desktop/X2K_Genetic_Algorithm/"

# Mix up and down, then write to new file, then get top n lines (to limit computational load on GA)
def uniqueKinasesGMT(testPercent, input, output1, output2, writeType="w", up_dn="dn"):
    with open(directory + "Validation/Perturbation_Data/"+input) as testgmt,\
        open(directory + "Validation/Perturbation_Data/"+output1, writeType) as subset1,\
        open(directory + "Validation/Perturbation_Data/"+output2, writeType) as subset2:
        testgmtLines = testgmt.readlines()
        #division = round(len(testgmtLines)*testPercent/100)
        division = round((570*.80)/2)
        count = 1
        kinaseRecorder_dn=[]; kinaseRecorder_up=[]

        for line in testgmtLines:
            kinase = line.split("")
            if count <= division and kinase not in kinaseRecorder_dn:
                newLine = line.split("\t")
                newLine[0] = newLine[0]+"_"+up_dn
                subset1.write("\t".join(newLine))
            if count > division:
                newLine = line.split("\t")
                newLine[0] = newLine[0]+"_"+up_dn
                subset2.write("\t".join(newLine))
            count += 1

# Make a fake GMT with random genes
## Matched to GEO dataset
def randomDEGsGMT(testGMT):
    randomGMT = testGMT.strip(".txt")+"-RandomDEGs.txt"
    with open(testGMT) as GEO80, open(randomGMT, "w") as GEO80_random:
        import pandas as pd
        from random import choice
        allGenes = pd.read_table("../X2K_Databases/General_Resources/Moshe_mapping/mappingFile_2017.txt", header=None)[0].str.strip().tolist()
        geo80 = GEO80.readlines()
        for line in geo80:
            splitLine = line.split("\t")
            DEGs = splitLine[2:]
            randomDEGs=[]
            for x in range(len(DEGs)):
                randomDEGs.append( choice(allGenes) )
            GEO80_random.write("\t".join([splitLine[0]]+randomDEGs)+"\n")
randomDEGsGMT("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt")
randomDEGsGMT("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET2.20per.txt")

## One row for each kinase in the kinome, with 300 DEGs each
def randomKinaseGMT(numGenes):
    import pandas as pd
    from random import choice
    kinases = pd.read_table("data/KEA/kea_ranks.txt", header=None)[0].tolist()
    allGenes = pd.read_table("../X2K_Databases/General_Resources/Moshe_mapping/mappingFile_2017.txt", header=None)[0].str.strip().tolist()
    with open("Validation/Perturbation_Data/Random_GMTs/allKinases_randomGenes.gmt","w") as GMT:
        for k in kinases:
            randomDEGs=[]
            for x in range(numGenes):
                randomDEGs.append(str(choice(allGenes)).upper())
            GMT.write(k+"_allKinases-randomGenes"+"\t\t"+"\t".join(randomDEGs)+"\n")
randomKinaseGMT(300)


from Python_scripts.X2K_Pipeline import X2K_fitness
def calculateFitness(population, allDataDF, genCount=1, fitnessMethod='targetAdjustedOverlap', testFitness=False):
    import pandas as pd; import numpy as np
    fitDF=pd.DataFrame(columns=allDataDF.columns)
    for i,indiv in enumerate(population):
        print("\n****** Generation "+str(genCount)+"  ::  Individual " + str(i+1) + " ******")
        # Delete .DS_Store files
        import os
        for root, dirs, files in os.walk(os.getcwd()):
            for file in files:
                if file.endswith(".DS_Store"):
                    os.remove(os.path.join(root, file))
        if testFitness == True:
            new_fitness = sum(map(int, indiv))  # Test fitness
            print('Fake fitness= ' + str(new_fitness))
            newBinary = indiv
            fitDF = fitDF.append(pd.DataFrame(np.column_stack([genCount,indiv,newBinary,fitnessMethod,new_fitness,"NA","NA","NA","NA","NA","NA","NA"]), columns=allDataDF.columns))
        else:
            # Calculate fitness (ONLY if it hasn't been previously calculated)
            if indiv not in allDataDF.newBinary.tolist():
                new_fitness, PPI_size, newBinary, chea_parameters, g2n_string, kea_string, baselineFitness, TKs, PKs = X2K_fitness(indiv,fitnessMethod)
                # Make new Dataframe
                fitDF = fitDF.append( pd.DataFrame(np.column_stack([genCount, indiv, newBinary,fitnessMethod, new_fitness, baselineFitness, \
                                                                    PPI_size, chea_parameters, g2n_string, kea_string, TKs, PKs]), columns=allDataDF.columns) )
                print("PPI Size = " +str(PPI_size))
            else:
                newBinary = indiv
                same = allDataDF[allDataDF.newBinary==indiv].iloc[0,:].copy() # Get the first individual that matches this one
                same['Generation'] = genCount
                fitDF = fitDF.append(same)
                print("[Using previously calculated fitness = "+ str(same.Fitness)+"]")
                print("[BaselineFitness = "+str(same.baselineFitness)+"]")
                print("[PPI size = "+str(same.PPI_size)+"]")
        print("OldBinary: " + indiv)
        print("NewBinary: "+ newBinary)
    # Convert cols to numeric
    fitDF['Fitness'] = pd.to_numeric(fitDF['Fitness'], errors='ignore')
    fitDF['baselineFitness'] = pd.to_numeric(fitDF['baselineFitness'], errors='ignore')
    fitDF['PPI_size'] = pd.to_numeric(fitDF['PPI_size'], errors='ignore')
    return fitDF

### ****** Run X2K on 1000 random GMTs to get the null rank distribution of each gene ****** ###

# Choose a single binary string to use (e.g. optimized)
#from Python_scripts.Extra_X2K_functions import getFittestIndividual
#selectedBinary = getFittestIndividual(randomGA_df)
selectedBinary='1101100000111110110001010011001101110010111'



def randomX2Kruns(selectedBinary):
    import pandas as pd
    import os
    from shutil import copyfile
    from time import sleep
    # Setup DF
    allDataDF = pd.DataFrame(
        columns=['Generation', 'oldBinary', 'newBinary', 'fitnessMethod', 'Fitness', 'baselineFitness', \
                 'PPI_size', 'CHEA_parameters', 'G2N_parameters', 'KEA_parameters', 'targetKinases', 'predictedKinases'])
    def addRandomGMT():
        # Clear old
        import os
        dir_name = "data/testgmt/"
        files = os.listdir(dir_name)
        for item in files:
            if item.endswith(".txt") or item.endswith(".gmt"):
                os.remove(os.path.join(dir_name, item))
        # Make new
        randomKinaseGMT(300)
        # Copy to testGMT folder
        copyfile("Validation/Perturbation_Data/Random_GMTs/allKinases_randomGenes.gmt", \
                 "data/testgmt/allKinases_randomGenes.gmt")
    ## Startup G2N and KEA
    try:
        os.popen('java -jar x2k_CHEA.jar')
        os.popen('java -jar x2k_G2N.jar')
        os.popen('java -jar x2k_KEA.jar')
    except:
        print("Already loaded")
    # Loop for 1000 random files
    for i in range(1000):
        print('Run: '+str(i))
        # Kill CHEA
        try:
            ind = os.popen('lsof -i tcp:5000').read().split(" ").index("schilder")
            PID = os.popen('lsof -i tcp:5000').read().split(" ")[ind-1]
            os.popen('kill ' + PID)
        except:
            print("No process to kill")
        # Replace testGMT with new file
        addRandomGMT()
        # Startup CHEA again
        os.popen('java -jar x2k_CHEA.jar')
        sleep(5)
        # Calculate fitness (after allowing CHEA some time to load up)
        fitDF = calculateFitness([selectedBinary], allDataDF)
        allDataDF = allDataDF.append(fitDF)
    # Save results
    GA_output_name = 'X2K_NullDist-PKlengthCorrected.npy'
    import os, numpy as np
    results_dir = 'GA_Results/Null_Distribution/'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    np.save(results_dir + GA_output_name, allDataDF)

    return allDataDF, selectedBinary

allDataDF, selectedBinary = randomX2Kruns(selectedBinary)
allDataDF = pd.read_csv('GA_Results/Null_Distribution/X2K_NullDist-PKlengthCorrected.csv')

# Get null distribution s
def nullDistributions(allDataDF):
    import pandas as pd
    KINASES = pd.read_table("data/KEA/kea_ranks.txt", header=None)[0].tolist()
    fullranks = list(range(1, len(KINASES) + 1))
    allPKs = allDataDF.predictedKinases.tolist()
    rankList=[]; predictedKinaseList=[]
    for ind in allPKs:
        print(ind)
        expts = ind.split(";")
        for exp in expts:
            expSplit = exp.split(",")
            misses = set(KINASES) - set(expSplit)
            remainingRanks = fullranks[len(set(expSplit))-1:]
            avgRank = sum(remainingRanks)/len(remainingRanks)
            for k in misses:
                rankList.append(avgRank)
                predictedKinaseList.append(k)
            if expSplit!=['']:
                for i,pkin in enumerate(expSplit):
                    rankList.append(i+1)
                    predictedKinaseList.append(pkin)
    df = pd.DataFrame({"predictedKinase":predictedKinaseList,"ranks":rankList})
    df['Rank'] = df["ranks"].astype(int)
    df['Counts'] = df.groupby(['predictedKinase','Rank']).transform('count')

    DF = df.copy()
    DF = DF.drop_duplicates()

    from scipy.stats.mstats import zscore
    Zdf = pd.DataFrame({'Rank':range(1,473+1)})
    Zdf = Zdf.astype(int) #[newDF['Rank']<=20]
    neverPredicted=[]
    for kinase in KINASES:
        print("Processing: "+str(kinase))
        sub = DF[DF['predictedKinase']==kinase].sort_values(by='Rank')
        #sub = sub[sub['Rank']<=20] # Remove any Ranks over 20
        if len(sub)>0:
            tmpMerge = pd.merge(Zdf, sub, on='Rank', how='left').fillna(0)
            tmpMerge['Zscore'] = zscore(tmpMerge['Counts'].tolist())
            Zdf[str(kinase)] = tmpMerge['Zscore']
        else:
            neverPredicted.append(kinase)
            print("No sig ranks for: "+str(kinase))
    Zdf.to_csv("Validation/Perturbation_Data/Random_GMTs/kinaseNullDistributions_all.csv", index=None)

def plotNullDist():
    # Plot
    ## Mean rank for every kinase
    import numpy as np
    import matplotlib.pyplot as plt
    DF['ranks'] = pd.to_numeric(DF['ranks'])
    avgRank = DF[DF['Rank']<=20].groupby('predictedKinase')['ranks'].mean().reset_index().iloc[1:,:]
    avgRank.plot(x='predictedKinase', y='ranks', kind='bar')
    plt.xticks(rotation=90, ha='center', fontsize=6)

    medianRank = DF[DF['Rank']<=20].groupby('predictedKinase')['ranks'].median().reset_index().iloc[1:,:]
    medianRank.plot(x='predictedKinase', y='ranks', kind='bar', label="", legend=False)
    plt.xticks(rotation=90,  ha='center', fontsize=6)
    plt.title("Median Rank for all Predicted Kinases")

    ## Single kinase dist
    def plotGene(gene):
        geneSub = DF[DF['predictedKinase']==gene][DF['Rank']<=20].copy()
        geneSub = geneSub.sort_values(by='Rank')
        geneSub.plot(x="Rank", y="Counts", kind="bar",color="pink", label="", legend=False)
        plt.title(str(gene)+" : rank frequencies")
        plt.xticks(np.arange(0, 20, 1))
    plotGene('TNIK')
    plotGene('ABL1')
    plotGene('TYK2')
    plotGene('VRK2')

    Zdf = pd.read_csv("Validation/Perturbation_Data/Random_GMTs/kinaseNullDistributions_all.csv", index_col=None)
    def plotZscore(gene):
        Zdf.plot(x='Rank',y=gene, kind='bar', color='c', legend='')
        plt.ylabel('Z-score of Rank Frequencies')
        plt.title(str(gene) + " : Null Distribution")
    plotZscore('TNIK')
    plotZscore('ABL1')
    plotZscore('TYK2')
    plotZscore('MAPK1')






