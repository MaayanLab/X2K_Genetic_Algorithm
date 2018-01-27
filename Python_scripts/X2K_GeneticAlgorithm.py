#############################################################################################
#############################################################################################
##################### X2K:Parameter Optimization Via Genetic Algorithm #####################
#############################################################################################
#############################################################################################

###################################
# 1. Create initial population
###################################
# binary = '0010100101111001111011111110'

from random import choice
def createPopulation(popSize, parameterLength):
    populationinit = []
    for i in range(popSize):
        populationinit.append(''.join(choice(('0', '1')) for _ in range(parameterLength)) )
        print(populationinit[i])
    return populationinit
# population = createPopulation(10, 27)

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
def calculateFitness(population, genCount='', fitness_method='simple'):
    indCount = 1
    fitness = []
    avg_PPI_size = []
    for i in range(len(population)):
        print()
        print("****** Generation "+str(genCount)+"  ::  Individual " + str(indCount) + " ******")
        print(population[i])
        # Delete .DS_Store files
        import os
        print("(Deleting .DS_Store files...)")
        for root, dirs, files in os.walk(os.getcwd()):
            for file in files:
                if file.endswith(".DS_Store"):
                    print(os.path.join(root, file))
                    os.remove(os.path.join(root, file))
        # Calculate fitness (ONLY if it hasn't been previously calculated)

        if population[i] not in fitnessDictionary:
            # FITNESS
            #new_fitness = sum(map(int, population[i])) # Test fitness
            X2K_output = X2K_fitness(population[i],fitness_method)
            new_fitness = X2K_output[0]  # Real fitness
            fitness.append(new_fitness)
            # AVERAGE PPI SIZE
            new_PPIsize = X2K_output[1]
            avg_PPI_size.append(new_PPIsize)
            # Store calculated values in dictionaries
            fitnessDictionary[population[i]] = new_fitness
            ppiSizeDictionary[population[i]] = new_PPIsize
            print("PPI Size = " +str(new_PPIsize))
        else:
            fitness.append( fitnessDictionary[population[i]] )
            avg_PPI_size.append(ppiSizeDictionary[population[i]])
            print("{Using previously calculated fitness: " + str( fitnessDictionary[population[i]] ) + "}")
            print("PPI Size = " + str(ppiSizeDictionary[population[i]]))
        indCount += 1

    return fitness, avg_PPI_size
# popFitness = calculateFitness(population)


###################################
# 3. Subset only the top fittest individuals
###################################
import numpy as np
def selectFittest(topNum, population, genCount='', fitness_method='simple'):
    calcFitness_output = calculateFitness(population, genCount, fitness_method)
    populationFitness = calcFitness_output[0]
    average_PPI_size = calcFitness_output[1]
    # Find indices of fittest
    fitnessIndices  = sorted(range(len(populationFitness)), key=lambda i: populationFitness[i])[-topNum:]
    fittest = list(np.array(population)[fitnessIndices])
    # Get Fitness of fittest
    fittestFitness = list(np.array(populationFitness)[fitnessIndices])
    print("Top fitnesses:  "+str(fittestFitness))

    return fittest, fittestFitness, populationFitness, average_PPI_size

# selectFittest_output = selectFittest(10,populationInit)
# Fittest = selectFittest_output[0]
# fittestFitness = selectFittest_output[1]


###################################
# 4. Crossover/breed fittest
###################################
###################################
# 5. Introduce random mutations
###################################

from random import random
def createChild(individual1, individual2, crossoverPoints):
    child_fragments = []
    chunkSize = int(len(individual1) / (crossoverPoints+1))
    ind1Split = [individual1[i:i + chunkSize] for i in range(0, len(individual1), chunkSize)]
    ind2Split = [individual2[i:i + chunkSize] for i in range(0, len(individual2), chunkSize)]
    for fragment in range(len(ind1Split)):
        if int(100 * random()) < 50: # Just randomly picks from ParentA or ParentB for each individual parameter
            child_fragments.append(ind1Split[fragment])
        else:
            child_fragments.append(ind2Split[fragment])
            ""
    #print(child_fragments)
    child = "".join(child_fragments)
    return child

# child = createChild( Fittest[0] , Fittest[1], 3)


def mutateChild(child, mutationRate):
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


def createChildren(numberOfChildren, Fittest, fittestFitness, breedingVariation, mutationRate):
    children = []
    breedingChances = []
    # Add noise to fitness score?
    for b in range(len(Fittest)):
        breedingChances.append(np.random.uniform(1 + breedingVariation, 1 - breedingVariation) * int(fittestFitness[b]))
        topBreeders = [x for _, x in sorted(zip(breedingChances, Fittest), reverse=True)]
    # Breed n times
    # 'Once you're in, you're in'. After selecting the top fittest individuals, it doesn't matter who is fitter within that group: everyone breeds with everyone else randomly
    for i in range(numberOfChildren):
        ind1 = int(random()*len(topBreeders))
        ind2 = int(random()*len(topBreeders))
        child = createChild(topBreeders[ind1], topBreeders[ind2], 3)
        # MUTATE the children!
        child = mutateChild(child, mutationRate)
        children.append(child)
    return children

# createChildren(1000, Fittest, fittestFitness, 0, .01)



###################################
# Genetic Algorithm
###################################

# def createPopulation(initialPopSize, parameterLength):
#     return GAresults_Subset1[0][0]
# createPopulation(1, 2)

def GAfunction(initialPopSize, parameterLength, numberOfGenerations, topNum, childrenPerGeneration, crossoverPoints, breedingVariation, mutationRate, includeFittestParents, fitness_method='simple', fitnessDictionary={}):
    print("Creating initial population...")
    population = createPopulation(initialPopSize, parameterLength)
    allPopulations = []
    allFitnesses = []
    averageFitness = []
    peakFitness = []
    average_PPI_sizes = []
    genCount = 0
    # Loop through n generationss
    for i in range(numberOfGenerations):
        genCount += 1
        print()
        print("+++++++ ANALYZING GENERATION: " + str(genCount) + " +++++++")

        # Stpre population
        allPopulations.append(population)
        # Calculate fitness & calculate fitness in one step
        fitnessOutput = selectFittest(topNum, population, genCount, fitness_method)
        fittest = fitnessOutput[0]
        fittestFitness = fitnessOutput[1]
        populationFitness = fitnessOutput[2]
        average_PPI_sizes.append(fitnessOutput[3])

        # Create new population from the fittest individuals
        ## Replace population and repeat
        population = createChildren(childrenPerGeneration, fittest, fittestFitness, breedingVariation, mutationRate)
        # Include fittest parents?
        if includeFittestParents > 0:
            population.extend( fittest[-includeFittestParents:] )

        # Store all fitnesses
        allFitnesses.append(populationFitness)
        # Store average fitness
        averageFitness.append( sum(populationFitness)*1.0 / len(populationFitness) ) #*1.0 turns numbers into floats instead of integers
        # Store peak fitness
        peakFitness.append( max(populationFitness) )
        # Store GA settings in in a dictionary
        GA_settings = {'initialPopSize':initialPopSize,'parameterLength':parameterLength,
                       'numberOfGenerations':numberOfGenerations, 'topNum':topNum,
                       'childrenPerGeneration':childrenPerGeneration, 'crossoverPoints':crossoverPoints,
                       'breedingVariation':breedingVariation, 'mutationRate':mutationRate,
                       'includeFittestParents':includeFittestParents
                       }
    return allPopulations, allFitnesses, averageFitness, peakFitness, GA_settings, average_PPI_sizes, fitnessDictionary


GAresults = GAfunction(initialPopSize=5, parameterLength=27, numberOfGenerations=2, topNum=2, childrenPerGeneration=5, crossoverPoints=3, breedingVariation=0, mutationRate=0.01, includeFittestParents=2, fitness_method='LINCS_L1000 + DrugRepurposingHub (adjustedScore)')


allPopulations = GAresults[0]# Get all populations
allFitnesses = GAresults[1] # Get all fitnesses
averageFitness = GAresults[2] # Get averageFitness per generation
peakFitness = GAresults[3] # Get the peakFitness per generation
numberOfGenerations = GAresults[4] # GA settings
average_PPI_sizes = GAresults[5]
fitnessDictionary_revived = GAresults[6]

import Python_scripts.Extra_X2K_functions as Ex
Ex.tell_parameters( Ex.getFittestIndividual(GAresults) )
Ex.parameterEvolutionPlot(GAresults)


###################################
# 7. Plot results
###################################
import matplotlib.pyplot as plt
# Plot averageFitness
x = range(len(averageFitness))
y1 = averageFitness
y2 = peakFitness
plt.plot(x, y1, 'bo--', markersize=3, label='Average Fitness')
plt.plot(x, y2, 'm^--', markersize=3, label='Peak Fitness')
plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.show()
plt.title('Fitness Over Generations')
plt.legend(loc='lower right')

# axes = plt.gca()
# axes.set_ylim([0,30])


# Plot the distribution of ALL fitnesses in 1st, 2nd & final populations
plt.hist( allFitnesses[0], bins=20) # 1st
plt.hist( allFitnesses[1], bins=20) # 2nd
plt.hist( allFitnesses[-1], bins=5) # Last
plt.ylabel('Frequency')
plt.title('Fitness distributions over generations')




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
    if item.endswith(".txt"):
        os.remove(os.path.join(dir_name, item))
## Replace it with subset 1
### Dataset.A: GEO KINASE PERTURBATION DATA
#copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET1.txt")
### Dataset.B: LINCS L1000 + DrugRepurposingHub
copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET1.txt", "data/testgmt/Chem_combo_DRH.kinaseInihibitors_SUBSET1.txt")
### Dataset.C: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET1.txt", "data/testgmt/LINCS-L1000_KINOMEscan_SUBSET1.txt")

## Run GA with Subset1
FITNESS_METHOD='simple'
GAresults_Subset1 = GAfunction(initialPopSize=100, parameterLength=35, numberOfGenerations=10, topNum=10, childrenPerGeneration=90, crossoverPoints=3, breedingVariation=0, mutationRate=0.01, includeFittestParents=10,\
                               fitness_method=FITNESS_METHOD)

#
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
copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET2.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET2.txt")
### Dataset.B: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/DrugRepurposingHub_filtered/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt", "data/testgmt/Chem_combo_DRH.kinaseInihibitors_SUBSET2.txt")
### Dataset.C: LINCS L1000 + DrugRepurposingHub
#copyfile("Validation/Perturbation_Data/LINCS_L1000_Chem/KinomeScan_filtered/LINCS-L1000_KINOMEscan_SUBSET2.txt", "data/testgmt/LINCS-L1000_KINOMEscan_SUBSET2.txt")

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
GA_output_name = 'GA_results.100pop.10gen.GEO.npy'
import os, numpy as np
results_dir = 'GA_Results/GEO/'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
np.save(results_dir+GA_output_name, list([GAresults_Subset1, allFitnesses_Subset2, \
                                          averageFitness_Subset2, peakFitness_Subset2]))

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


