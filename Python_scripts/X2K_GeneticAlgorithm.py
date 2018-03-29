#############################################################################################
##################### X2K:Parameter Optimization Via Genetic Algorithm #####################
#############################################################################################

###################################
# 1. Create initial population
###################################
from random import choice
def createPopulation(popSize, binaryStringLength):
    populationinit = []
    for i in range(popSize):
        populationinit.append(''.join(choice(('0', '1')) for _ in range(binaryStringLength)) )
        print(populationinit[i])
    return populationinit
# population = createPopulation(20, 43)

###################################
# 2. Calculate fitness
###################################
try:
    del X2K_fitness
except NameError:
    pass
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
# allDataDF = pd.DataFrame(columns=['Generation', 'oldBinary', 'newBinary', 'fitnessMethod', 'Fitness', 'baselineFitness', 'PPI_size', 'CHEA_parameters', 'G2N_parameters', 'KEA_parameters', 'targetKinases', 'predictedKinases'])
# fitDF = calculateFitness(population, allDataDF)



###################################
# 3. Subset only the top fittest individuals
###################################

def selectFittest(topNum, fitDF, selectionMethod='Fitness-proportional'):
    import pandas as pd
    if selectionMethod == 'Fitness-proportional':
        fittestDF = fitDF.sort_values(by=['Fitness'],ascending=False).iloc[:topNum,:]
        print("Top fitnesses:  " + str(fittestDF["Fitness"].values))
    # Tournament selection (less stringent)
    ## Split the population into equal subgroups, and then select the fittest individual from each group
    elif selectionMethod == 'Tournament':
        fittestDF=pd.DataFrame()
        if fitDF.shape[0] % topNum!=0:
            print("Tournament selection requires that populationSize/topNum and childrenPerGeneration/topNum are both whole numbers.")
        subsetSize = int( fitDF.shape[0] / topNum )
        for t in range(topNum):
            subDF = fitDF.sample(n=subsetSize, replace=False)
            fittestDF = fittestDF.append( subDF.sort_values(by=['Fitness'], ascending=False).iloc[0,:].copy())
    elif selectionMethod == 'mixedTournament':
        if fitDF.shape[0]%topNum!=0 or topNum%2!=0:
            print("Tournament selection requires that populationSize/topNum, childrenPerGeneration/topNum, and topNum/2 to be whole numbers.")
        topNumHalf = int(topNum/2)
        sortedDF = fitDF.sort_values(by=['Fitness'], ascending=False).copy()
        # The first half of the new pop are the fittest parents overall
        fittestDF = sortedDF.iloc[:topNumHalf, :].copy()
        # Then run Tournament selection on the rest of the population to get the other half of the new pop
        everybodyElse = sortedDF.iloc[topNumHalf:, :].copy()
        subsetSize = int(everybodyElse.shape[0] / topNumHalf)
        for t in range(topNumHalf):
            subDF = everybodyElse.sample(n=subsetSize, replace=False)
            fittestDF = fittestDF.append(subDF.sort_values(by=['Fitness'], ascending=False).iloc[0, :].copy())
    else:
        print("Use viable 'selectionMethod'")
    return fittestDF
# fittestDF = selectFittest(5, fitDF, selectionMethod="Tournament")



###################################
# 4. Crossover/breed fittest
###################################
###################################
# 5. Introduce random mutations
###################################
# individual1=  population[0]
# individual2=  population[1]
# crossoverPoints=8
def createChild(individual1, individual2, crossoverPoints, crossoverLocations="evenlyDistributed"):
    if crossoverLocations=="evenlyDistributed":
        chunkSize = int(len(individual1) / (crossoverPoints+1))
        ind1Split = [individual1[i:i + chunkSize] for i in range(0, len(individual1), chunkSize)]
        ind2Split = [individual2[i:i + chunkSize] for i in range(0, len(individual2), chunkSize)]
    elif crossoverLocations=='random':
        from random import sample
        cutpoints = sorted(sample(range(1, len(individual1)-1), crossoverPoints)) # randomly generate n non-overlapping numbers
        def splitParent(parent, cutpoints):
            indSplit=[]
            for i,num in enumerate(cutpoints):
                #print("**Cutpoint index= "+str(i))
                if i == 0: # If it's the first cutpoint, take all values up to the first index+1
                    start = 0
                    end = num
                else:
                    start = cutpoints[i-1]
                    end = num
                segment = parent[start:end]
                #print("Cutpoint= " + str(start) + " : " + str(end))
                #print("------- "+segment+" -------")
                indSplit.append(segment)
            # Add the very last segment
            indSplit.append(parent[cutpoints[-1]:])
            return indSplit
        ind1Split = splitParent(individual1, cutpoints)
        ind2Split = splitParent(individual2, cutpoints)
    # Put together the new child
    from random import random
    childFragments=[]
    for fragment in range(len(ind1Split)):
        if int(100*random()) < 50: # Just randomly picks from ParentA or ParentB for each individual parameter
            childFragments.append(ind1Split[fragment])
        else:
            childFragments.append(ind2Split[fragment])
    child = "".join(childFragments)
    return child
# child = createChild( fittest[0] , fittest[1], 3)


def mutateChild(child, mutationRate):
    from random import random
    mutant = ''
    for bit, val in enumerate(child):
        rando = random()
        if rando <= mutationRate and val == '1':
            mutant = str(str(mutant) + '0')
        elif rando <= mutationRate and val == '0':
            mutant = str(str(mutant) + '1')
        else:
            mutant = str(str(mutant) + str(val))
    return mutant

def createChildren(numberOfChildren, fittestDF, mutationRate, breedingVariation, crossoverLocations):
    from random import random
    fittest = fittestDF.newBinary.tolist()
    #breedingChances = []
    # Add noise to fitness score?
    # for b in range(len(Fittest)):
    #     breedingChances.append(np.random.uniform(1 + breedingVariation, 1 - breedingVariation) * int(fittestFitness[b]))
    #     topBreeders = [x for _, x in sorted(zip(breedingChances, Fittest), reverse=True)]
    # Breed n times
    # 'Once you're in, you're in'. After selecting the top fittest individuals, it doesn't matter who is fitter within that group: everyone breeds with everyone else randomly
    children = []
    for i in range(numberOfChildren):
        ind1 = int(random()*len(fittest))
        ind2 = int(random()*len(fittest))
        child = createChild(fittest[ind1], fittest[ind2], 3, crossoverLocations)
        # MUTATE the children!
        child = mutateChild(child, mutationRate)
        children.append(child)
    return children
# createChildren(100, fitDF, .01)



###################################
# Genetic Algorithm
###################################


def GAfunction(initialPopSize, binaryStringLength, numberOfGenerations, topNum, childrenPerGeneration, crossoverPoints, crossoverLocations, breedingVariation, mutationRate, includeFittestParents, fitnessMethod, selectionMethod, setInitialPopulation=False):
    # Store GA settings in in a dictionary
    GAsettings = dict(zip(
        ['initialPopSize', 'binaryStringLength', 'numberOfGenerations', 'topNum', 'childrenPerGeneration',\
         'crossoverPoints', 'crossoverLocations','breedingVariation', 'mutationRate', 'includeFittestParents', 'fitnessMethod','selectionMethod', 'setInitialPopulation'], \
        [initialPopSize, binaryStringLength, numberOfGenerations, topNum, childrenPerGeneration, crossoverPoints,crossoverLocations,\
         breedingVariation, mutationRate, includeFittestParents, fitnessMethod, selectionMethod, setInitialPopulation]))
    # Store all of the relevant results in one dataframe
    import pandas as pd
    allDataDF = pd.DataFrame(columns=['Generation', 'oldBinary', 'newBinary', 'fitnessMethod', 'Fitness', 'baselineFitness', \
                                      'PPI_size', 'CHEA_parameters', 'G2N_parameters', 'KEA_parameters', 'targetKinases','predictedKinases'])
    print("Creating initial population...")
    if setInitialPopulation==False:
        population = createPopulation(initialPopSize, binaryStringLength)
    else:
        population = setInitialPopulation
    genCount = 1
    # Loop through n generations:
    for i in range(numberOfGenerations):
        print()
        print("+++++++ ANALYZING GENERATION: " + str(genCount) + " +++++++")
        # Calculate fitness
        fitDF = calculateFitness(population, allDataDF, genCount, fitnessMethod)
        allDataDF = allDataDF.append(fitDF)
        #allDataDF = allDataDF.append(allDataDF)
        # Select fittest
        fittestDF = selectFittest(topNum, fitDF, selectionMethod)
        # Create new population from the fittest individuals
        population = createChildren(childrenPerGeneration, fittestDF, mutationRate, breedingVariation, crossoverLocations)
        # Include fittest parents?
        if includeFittestParents > 0:
            population.extend( fittestDF.newBinary.values[:includeFittestParents].tolist() )
        genCount += 1
    return allDataDF, GAsettings


GA_df, GAsettings = GAfunction(initialPopSize=10, binaryStringLength=43, numberOfGenerations=2, topNum=2, childrenPerGeneration=8, crossoverPoints=6, crossoverLocations='random',\
                       breedingVariation=0, mutationRate=0.01, includeFittestParents=2, fitnessMethod='rankCorrectedTAO', selectionMethod='mixedTournament', setInitialPopulation=False)
                        #targetAdjustedOverlap


# ***************** SUBSET 1 *****************
# Put Subset1 GMT file into testgmt folder

import os
from shutil import copyfile
## Delete whatever file is in testgmt
dir_name = "data/testgmt/"
files = os.listdir(dir_name)
for item in files:
    if item.endswith(".txt") or item.endswith(".gmt"):
        os.remove(os.path.join(dir_name, item))
    ## Replace it with subset 1
    ### Dataset.A: GEO KINASE PERTURBATION DATA
copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET1.80per.txt")
# Down genes only (performs MUCH better than up genes, which consistently return 0% of all kinase perturbation experiments in GEO)
#copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_up.gmt", "data/testgmt/Kinase_Perturbations_from_GEO_up.gmt")
    ### Dataset.B: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET1.80per.txt", "data/testgmt/Chem_combo_DRH.kinaseInihibitors_SUBSET1.80per.txt")
    ### Dataset.C: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET1.txt", "data/testgmt/LINCS-L1000_KINOMEscan_SUBSET1.txt")
    ### Dataset D: LINCS L1000 + Kinase Perturbations
#copyfile("LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_SUBSET1.80per.txt","data/testgmt/LINCS_L1000_Kinase_Perturbations_SUBSET1.80per.txt")
    ### COMBINED dataset: GEO-KinasePert + L1000-DRH
#copyfile("Validation/Perturbation_Data/Combined/GEO-KinasePert_L1000-DRH_SUBSET1-80per.txt", "data/testgmt/GEO-KinasePert_L1000-DRH_SUBSET1-80per.txt")
    ### Dataset X: Randomly chosen sets of DEGS:
#copyfile("Validation/Perturbation_Data/Random_GMTs/allKinases_randomGenes.gmt", "data/testgmt/allKinases_randomGenes.gmt")

## Run GA with Subset1
Subset1_df, GAsettings = GAfunction(initialPopSize=100, binaryStringLength=43, numberOfGenerations=20, topNum=10, childrenPerGeneration=90,\
                                    crossoverPoints=5, crossoverLocations='random', breedingVariation=0, mutationRate=0.01, includeFittestParents=10, \
                                    fitnessMethod='targetAdjustedOverlap_outputLengthCorrection', selectionMethod='mixedTournament') # modifiedRBO, targetAdjustedOverlap, rankCorrectedTAO



# ***************** SUBSET 2 *****************
# Put Subset2 GMT file into testgmt folder
## Delete whatever file is in testgmt
import os
from shutil import copyfile
dir_name = "data/testgmt/"
files = os.listdir(dir_name)
for item in files:
    if item.endswith(".txt") or item.endswith(".gmt"):
        os.remove(os.path.join(dir_name, item))
    ## Replace it with subset 2
    ### Dataset.A: GEO KINASE PERTURBATION DATA
copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET2.20per.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET2.20per.txt")
    ### Dataset.B: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt", "data/testgmt/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt")
    ### Dataset.C: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET2.txt", "data/testgmt/LINCS-L1000_KINOMEscan_SUBSET2.txt")
    ### Dataset D: LINCS L1000 + Kinase Perturbations
#copyfile("LINCS_L1000_Kinase_pert/LINCS_L1000_Kinase_Perturbations_SUBSET2.20per.txt","data/testgmt/LINCS_L1000_Kinase_Perturbations_SUBSET2.20per.txt")
    ### COMBINED dataset: GEO-KinasePert + L1000-DRH
#copyfile("Validation/Perturbation_Data/Combined/GEO-KinasePert_L1000-DRH_SUBSET2-20per.txt", "data/testgmt/GEO-KinasePert_L1000-DRH_SUBSET2-20per.txt")
    ### Dataset X: Randomly chosen sets of DEGS:
#copyfile("Validation/Perturbation_Data/Random_GMTs/allKinases_randomGenes2.gmt", "data/testgmt/allKinases_randomGenes2.gmt")



## Run GA with Subset2
# xxxxxxxx MAKE SURE YOU SHUT DOWN CHEA FIRST OR ELSE IT WON"T PRE-LOAD THE NEW TESTGMT!!!!!!!! xxxxxxxx
def runSubset2(Subset1_df, GAsettings):
    import pandas as pd
    allDataDF2=pd.DataFrame(columns=Subset1_df.columns)
    for gen in Subset1_df.Generation.unique().tolist():
        genPop = Subset1_df[Subset1_df.Generation==gen].newBinary.values.tolist()
        allDataDF2 = allDataDF2.append( calculateFitness(genPop, genCount=gen, fitnessMethod=GAsettings['fitnessMethod'], allDataDF=allDataDF2))
    return allDataDF2
Subset2_df = runSubset2(Subset1_df, GAsettings)

# Save/load GAresults as file
# Save
GA_output_name = 'GAresults_GEO-PKlengthCorrected.npy'
import os, numpy as np
results_dir = 'GA_Results/GEO/'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
np.save(results_dir+GA_output_name, [Subset1_df, Subset2_df, GAsettings])


# # Load
# import numpy as np
# results_file = 'GA_Results/'+GA_output_name
# np.load(results_file)
#


import os
from shutil import copyfile
dir_name = "data/testgmt/"
files = os.listdir(dir_name)
for item in files:
    if item.endswith(".txt"):
        os.remove(os.path.join(dir_name, item))
copyfile("X2K_Genetic_Algorithm/Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.80per-SCRAMBLED.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET1.80per-SCRAMBLED.txt")
## Run GA with scrambled-DEGs in testGMT:
def runSubset2(Subset1_df, GAsettings):
    import pandas as pd
    allDataDF2=pd.DataFrame(columns=Subset1_df.columns)
    for gen in Subset1_df.Generation.unique().tolist():
        genPop = Subset1_df[Subset1_df.Generation==gen].newBinary.values.tolist()
        allDataDF2 = allDataDF2.append( calculateFitness(genPop, genCount=gen, fitnessMethod=GAsettings['fitnessMethod'], allDataDF=allDataDF2))
    return allDataDF2
Subset2_df = runSubset2(Subset1_df, GAsettings)