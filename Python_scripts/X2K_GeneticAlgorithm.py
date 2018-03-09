#############################################################################################
#############################################################################################
##################### X2K:Parameter Optimization Via Genetic Algorithm #####################
#############################################################################################
#############################################################################################

###################################
# 0. Initialize X2K
###################################
# import os
# # Initialize X2K steps
# os.system('java -jar x2k_CHEA.jar')
# os.system('java -jar x2k_G2N.jar')
# os.system('java -jar x2k_KEA.jar')
# # Kill processes
# os.system('pkill 5000')
# os.system('pkill 5001')
# os.system('pkill 5002')
#



###################################
# 1. Create initial population
###################################
# binary = '0010100101111001111011111110'

from random import choice
def createPopulation(popSize, binaryStringLength):
    populationinit = []
    for i in range(popSize):
        populationinit.append(''.join(choice(('0', '1')) for _ in range(binaryStringLength)) )
        print(populationinit[i])
    return populationinit
# population = createPopulation(10, 43)

###################################
# 2. Calculate fitness
###################################
try:
    del X2K_fitness
except NameError:
    pass
from Python_scripts.X2K_Pipeline import X2K_fitness

fitnessDictionary = {}
ppiSizeDictionary = {}
def calculateFitness(population, genCount='', fitness_method='target-adjusted overlap'):
    testFitness = False
    populationFitness=[]; avg_PPI_size=[]; newPopulation=[]
    for i,indiv in enumerate(population):
        print()
        print("****** Generation "+str(genCount)+"  ::  Individual " + str(i+1) + " ******")
        print(indiv)
        # Delete .DS_Store files
        import os
        print("(Deleting .DS_Store files...)")
        for root, dirs, files in os.walk(os.getcwd()):
            for file in files:
                if file.endswith(".DS_Store"):
                    #print(os.path.join(root, file))
                    os.remove(os.path.join(root, file))
        # Calculate fitness (ONLY if it hasn't been previously calculated)
        if indiv not in fitnessDictionary:
            if testFitness==True:
                new_fitness = sum(map(int,indiv)) # Test fitness
                new_PPIsize=1; newBinary=indiv
                print('Fake fitness= '+str(new_fitness))
            else:
                # REAL FITNESS
                new_fitness, new_PPIsize, newBinary = X2K_fitness(indiv,fitness_method)
            populationFitness.append(new_fitness)
            avg_PPI_size.append(new_PPIsize)
            newPopulation.append(newBinary)
            # Store calculated values in dictionaries
            fitnessDictionary[indiv] = new_fitness
            ppiSizeDictionary[indiv] = new_PPIsize
            print("PPI Size = " +str(new_PPIsize))
        else:
            populationFitness.append( fitnessDictionary[indiv] )
            avg_PPI_size.append(ppiSizeDictionary[indiv])
            newBinary=indiv
            newPopulation.append(newBinary)
            print("{Using previously calculated fitness: " + str( fitnessDictionary[indiv] ) + "}")
            print("PPI Size = " + str(ppiSizeDictionary[indiv]))
    return populationFitness, avg_PPI_size, newPopulation
# populationFitness, avg_PPI_size, newPopulation = calculateFitness(population)


###################################
# 3. Subset only the top fittest individuals
###################################

def selectFittest(topNum, newPopulation, populationFitness, selectionMethod='Fitness-proportional'):
    import numpy as np
    import pandas as pd
    if selectionMethod == 'Fitness-proportional':
        # fitnessIndices = sorted(range(len(populationFitness)), key=lambda i: populationFitness[i])[-topNum:]
        # fittest = list(np.array(newPopulation)[fitnessIndices])
        # # Get Fitness of fittest
        # fittestFitness = list(np.array(populationFitness)[fitnessIndices])

        # DAtaframe method
        df = pd.DataFrame({'newPopulation':newPopulation,'populationFitness':populationFitness})
        df['populationFitness'] = pd.to_numeric(df['populationFitness'],  errors='ignore')
        sortedDF = df.sort_values(by=['populationFitness'],ascending=False)
        fittest =  sortedDF.newPopulation.tolist()[0:topNum]
        fittestFitness = sortedDF.populationFitness.tolist()[0:topNum]
        print("Top fitnesses:  " + str(fittestFitness))
    # elif selectionMethod == 'Tournament':
    #     import random; import operator
        # dataDict = {z[0]:list(z[1:]) for z in zip( list(range(0,len(population)+1)), population, populationFitness,average_PPI_size)}
        # subsetSize = int( len(population) / topNum )
        # keyList = list(dataDict.keys())
        # # Break dict into equal subsets, and get individual with top fitness from each subset
        # for num in topNum:
        #     random.shuffle(keyList)
        #     keySubset = keyList[:subsetSize]
        #     dictSubset = {k: dataDict[k] for k in (keySubset)}
        #     # Get top val
        #     #dicts
        #     max(dictSubset.items(), key=operator.itemgetter(2))[0]
        #     # Remove keys from list to make sure they don't get repeated
        #     keyList = list( set(keyList) - set(keySubset))
    return fittest, fittestFitness

# fittest, fittestFitness = selectFittest(10,newPopulation, populationFitness)



###################################
# 4. Crossover/breed fittest
###################################
###################################
# 5. Introduce random mutations
###################################
# individual1=population[0]
# individual2=population[6]
# crossoverPoints=8
def createChild(individual1, individual2, crossoverPoints, crosspointLocations="evenly distributed"):
    from random import sample, random
    child_fragments = []
    if crosspointLocations=="evenly distributed":
        chunkSize = int(len(individual1) / (crossoverPoints+1))
        ind1Split = [individual1[i:i + chunkSize] for i in range(0, len(individual1), chunkSize)]
        ind2Split = [individual2[i:i + chunkSize] for i in range(0, len(individual2), chunkSize)]
    else:
        ind1Split=[]; ind2Split=[]
        cutpoints = sorted(sample(range(1, len(individual1)-1), crossoverPoints)) # randomly generate n non-overlapping numbers
        ind1Split = [individual1[num:cutpoints[i+1]] for i,num in enumerate(cutpoints)]
        ind2Split =[]
        for i,num in enumerate(cutpoints):
            print(individual1[num:cutpoints[i+1]])
            ind1Split.append(individual1[num:cutpoints[i+1]])

    for fragment in range(len(ind1Split)):
        if int(100 * random()) < 50: # Just randomly picks from ParentA or ParentB for each individual parameter
            child_fragments.append(ind1Split[fragment])
        else:
            child_fragments.append(ind2Split[fragment])
    child = "".join(child_fragments)
    return child

# child = createChild( fittest[0] , fittest[1], 3)


def mutateChild(child, mutationRate):
    from random import random
    mutant = ''
    for bit, val in enumerate(child):
        rando = random()
        if rando <= mutationRate and val == '1':
            mutant = str(str(mutant) + '0')
        #print("1 ==> 0")
        elif rando <= mutationRate and val == '0':
            mutant = str(str(mutant) + '1')
        #print("0 ==> 1")
        else:
            mutant = str(str(mutant) + str(val))
        #print("NO MUTATION")
    #print(str(mutant))
    #print(str(child))
    return mutant

def createChildren(numberOfChildren, fittest, fittestFitness, breedingVariation, mutationRate):
    import numpy as np
    from random import random
    children = []
    #breedingChances = []
    # Add noise to fitness score?
    # for b in range(len(Fittest)):
    #     breedingChances.append(np.random.uniform(1 + breedingVariation, 1 - breedingVariation) * int(fittestFitness[b]))
    #     topBreeders = [x for _, x in sorted(zip(breedingChances, Fittest), reverse=True)]
    # Breed n times
    # 'Once you're in, you're in'. After selecting the top fittest individuals, it doesn't matter who is fitter within that group: everyone breeds with everyone else randomly
    for i in range(numberOfChildren):
        ind1 = int(random()*len(fittest))
        ind2 = int(random()*len(fittest))
        child = createChild(fittest[ind1], fittest[ind2], 3)
        # MUTATE the children!
        child = mutateChild(child, mutationRate)
        children.append(child)
    return children

# createChildren(1000, fittest, fittestFitness, 0, .01)



###################################
# Genetic Algorithm
###################################


def GAfunction(initialPopSize, binaryStringLength, numberOfGenerations, topNum, childrenPerGeneration, crossoverPoints, breedingVariation, mutationRate, includeFittestParents, fitness_method, fitnessDictionary={}):
    print("Creating initial population...")
    population = createPopulation(initialPopSize, binaryStringLength)
    allPopulations=[]; allFitnesses=[]; averageFitness=[]; peakFitness=[]; average_PPI_sizes=[]
    genCount = 0
    from random import random
    # Loop through n generationss
    for i in range(numberOfGenerations):
        genCount += 1

        print()
        print("+++++++ ANALYZING GENERATION: " + str(genCount) + " +++++++")
        # Calculate fitness
        populationFitness, average_PPI_size, newPopulation = calculateFitness(population=population, genCount=genCount, fitness_method=fitness_method)
        allFitnesses.append(populationFitness) # Store all fitness
        average_PPI_sizes.append(average_PPI_size)
        allPopulations.append(newPopulation) # Store ****POST-SHUFFLED**** population!!!!

        # Select fittest
        fittest, fittestFitness = selectFittest(topNum=topNum, newPopulation=newPopulation, populationFitness=populationFitness)
        # Create new population from the fittest individuals
        population = createChildren(numberOfChildren=childrenPerGeneration, fittest=fittest, fittestFitness=fittestFitness, breedingVariation=breedingVariation, mutationRate=mutationRate)
        # Include fittest parents?
        if includeFittestParents > 0:
            population.extend( fittest[:includeFittestParents] )

        # Store average fitness
        averageFitness.append( sum(populationFitness)*1.0 / len(populationFitness) )
        # Store peak fitness
        peakFitness.append( max(populationFitness) )
        # Store GA settings in in a dictionary
        GA_settings = {'initialPopSize':initialPopSize,'parameterLength':binaryStringLength,
                       'numberOfGenerations':numberOfGenerations, 'topNum':topNum,
                       'childrenPerGeneration':childrenPerGeneration, 'crossoverPoints':crossoverPoints,
                       'breedingVariation':breedingVariation, 'mutationRate':mutationRate,
                       'includeFittestParents':includeFittestParents
                       }
    return allPopulations, allFitnesses, averageFitness, peakFitness, GA_settings, average_PPI_sizes, fitnessDictionary


GAresults = GAfunction(initialPopSize=10, binaryStringLength=43, numberOfGenerations=5, topNum=2, childrenPerGeneration=8, crossoverPoints=3, breedingVariation=0, mutationRate=0.01, includeFittestParents=2, \
                       fitness_method='target-adjusted overlap')


allPopulations = GAresults[0]# Get all populations
allFitnesses = GAresults[1] # Get all fitnesses
averageFitness = GAresults[2] # Get averageFitness per generation
peakFitness = GAresults[3] # Get the peakFitness per generation
GA_settings = GAresults[4] # GA settings
average_PPI_sizes = GAresults[5] # Average_PPI_sizes
fitnessDictionary_revived = GAresults[6]

import Python_scripts.Extra_X2K_functions as Ex
Ex.tell_parameters( Ex.getFittestIndividual(GAresults) )
Ex.parameterEvolutionPlot(GAresults)


###################################
# 7. Plot results
###################################
import matplotlib.pyplot as plt
# Plot averageFitness
y1 = averageFitness
y2 = peakFitness
x = range(1,len(y1)+1)

plt.plot(x, y1, 'bo--', markersize=3, label='Average Fitness')
plt.plot(x, y2, 'm^--', markersize=3, label='Peak Fitness')
plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.show()
plt.title('Fitness Over Generations')
plt.legend(loc='lower right')

# axes = plt.gca()
# axes.set_ylim([0,30])



######################################
# Genetic Algorithm: Overfitting Test
######################################
# Calculate fitness for all individuals using only 250 of testgmt experiments, and again in other 250 of experiments
    # Breed & mutate the fittest individuals
# Calculate fitness of the new population using the 1st set of 250 and select top X fittest individuals
    # Calculate fitness of new pop using 2nd set of 250
# Repeat

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
# COMBINED dataset: GEO-KinasePert + L1000-DRH
#copyfile("Validation/Perturbation_Data/Combined/GEO-KinasePert_L1000-DRH_SUBSET1-80per.txt", "data/testgmt/GEO-KinasePert_L1000-DRH_SUBSET1-80per.txt")

## Run GA with Subset1
FITNESS_METHOD='target-adjusted overlap'
GAresults_Subset1 = GAfunction(initialPopSize=100, binaryStringLength=43, numberOfGenerations=20, topNum=10, childrenPerGeneration=90, crossoverPoints=5, breedingVariation=0, mutationRate=0.01, includeFittestParents=10,\
                               fitness_method=FITNESS_METHOD)

# # Recover fitnessDictionary
# def recoverFitnessDictionary(GAresults):
#     all_Fitness = GAresults[1]
#     fitnessDictionary_Subset1 = {}
#     allFitnesses = []
#     for sublist in all_Fitness:
#         for item in sublist:
#             allFitnesses.append(item)
#     allPopulations = GAresults[0]
#     allBinaries = []
#     for sublist in allPopulations:
#         for item in sublist:
#             allBinaries.append(item)
#     FitnessDictionary = dict(zip(allBinaries, allFitnesses))
#     return FitnessDictionary
#
# fitnessDictionary_Subset1 = recoverFitnessDictionary(GAresults_Subset1)
#


# ***************** SUBSET 2 *****************
# Put Subset2 GMT file into testgmt folder
## Delete whatever file is in testgmt
import os
from shutil import copyfile
dir_name = "data/testgmt/"
files = os.listdir(dir_name)
for item in files:
    if item.endswith(".txt"):
        os.remove(os.path.join(dir_name, item))
## Replace it with subset 2
### Dataset.A: GEO KINASE PERTURBATION DATA
copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET2.20per.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET2.20per.txt")
### Dataset.B: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt", "data/testgmt/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt")
### Dataset.C: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET2.txt", "data/testgmt/LINCS-L1000_KINOMEscan_SUBSET2.txt")
# COMBINED dataset: GEO-KinasePert + L1000-DRH
#copyfile("Validation/Perturbation_Data/Combined/GEO-KinasePert_L1000-DRH_SUBSET2-20per.txt", "data/testgmt/GEO-KinasePert_L1000-DRH_SUBSET2-20per.txt")

## Run GA with Subset2
# xxxxxxxx MAKE SURE YOU SHUT DOWN CHEA FIRST OR ELSE IT WON"T PRE-LOAD THE NEW TESTGMT!!!!!!!! xxxxxxxx
fitnessDictionary = {} # Reset fitnessDict
allFitnesses_Subset2 = []
averageFitness_Subset2 = []
genCount=0
average_PPI_sizes_Subset2 = []
for generation in GAresults_Subset1[0]:
    genCount+=1
    calcFitness_output = calculateFitness(generation, genCount=genCount, fitness_method=FITNESS_METHOD)
    new_fitnesses  = calcFitness_output[0]
    average_PPI_sizes_Subset2.append(calcFitness_output[1])
    allFitnesses_Subset2.append(new_fitnesses)
    averageFitness_Subset2.append( sum(new_fitnesses) / len(new_fitnesses) )
peakFitness_Subset2 = []
for generation in allFitnesses_Subset2:
    peakFitness_Subset2.append(sorted(generation, reverse=True)[0])


# Save/load GAresults as file
# Save
GA_output_name = 'GA_results_L1000-DRH.wPPIlimiters_ShuffleCorrected.npy'
import os, numpy as np
results_dir = 'GA_Results/GEO/'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
np.save(results_dir+GA_output_name, [GAresults_Subset1, allFitnesses_Subset2, \
                                          averageFitness_Subset2, peakFitness_Subset2, average_PPI_sizes_Subset2])

# Load
import numpy as np
results_file = 'GA_Results/'+GA_output_name
np.load(results_file)





# Plot trained vs. untrained
## Average fitness
import matplotlib.pyplot as plt
y_s1 = GAresults_Subset1[2] # average fitness for each generation in TRAINED data
y_s2 = averageFitness_Subset2
x = range(len(y_s1))
plt.plot(x, y_s1, 'c1--', markersize=7, label="Trained Results")
plt.plot(x, y_s2, 'gx--', markersize=7, label="Untrained Results")
plt.xlabel('Generation')
plt.ylabel('Average Fitness')
plt.legend(loc='lower right')

# Peak fitness
import matplotlib.pyplot as plt
y_s1 = GAresults_Subset1[3] # average fitness for each generation in TRAINED data
y_s2 = peakFitness_Subset2
x = range(len(y_s1))
plt.plot(x, y_s1, 'c1--', markersize=7, label="Trained Results")
plt.plot(x, y_s2, 'gx--', markersize=7, label="Untrained Results")
plt.xlabel('Generation')
plt.ylabel('Peak Fitness')
plt.legend(loc='lower right')