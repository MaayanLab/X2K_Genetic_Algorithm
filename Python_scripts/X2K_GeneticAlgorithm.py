#############################################################################################
#############################################################################################
##################### X2K:Parameter Optimization Via Genetic Algorithm #####################
#############################################################################################
#############################################################################################

###################################
# 1. Create initial population
###################################

from random import choice
def createPopulation(popSize, parameterLength):
    populationinit = []
    for i in range(popSize):
        populationinit.append(''.join(choice(('0', '1')) for _ in range(parameterLength)) )
        print(populationinit[i])
    return populationinit
parameterLength = 27
populationInit = createPopulation(10,parameterLength)

###################################
# 2. Calculate fitness
###################################
try:
    del X2K_fitness
except NameError:
    pass

from Python_scripts.X2K_Pipeline import X2K_fitness


fitnessDictionary = {}
def calculateFitness(population):
    indCount = 1
    fitness = []
    for i in range(len(population)):
        print()
        print("****** Individual " + str(indCount) + " ******")
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
            #new_fitness = sum(map(int, population[i])) # Test fitness
            new_fitness = X2K_fitness(population[i]) # Real fitness
            fitness.append( new_fitness )
            fitnessDictionary[population[i]] = new_fitness
        else:
            fitness.append( fitnessDictionary[population[i]] )
            print("{Using previously calculated fitness: " + str( fitnessDictionary[population[i]] ) + "}")
        indCount += 1
    return fitness
# popFitness = calculateFitness(populationInit)


###################################
# 3. Subset only the top fittest individuals
###################################
import numpy as np
def selectFittest(topNum, population):
    fitness = populationFitness = calculateFitness(population)
    # Find indices of fittest
    fitnessIndices  = sorted(range(len(fitness)), key=lambda i: fitness[i])[-topNum:]
    fittest = list(np.array(population)[fitnessIndices])
    # Get Fitness of fittest
    fittestFitness = list(np.array(fitness)[fitnessIndices])
    print("Top fitnesses:  "+str(fittestFitness))

    return fittest, fittestFitness, populationFitness

# electFittest_output = selectFittest(10,populationInit)
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

def GAfunction(initialPopSize, parameterLength, numberOfGenerations, topNum, childrenPerGeneration, crossoverPoints, breedingVariation, mutationRate, includeFittestParents):
    print("Creating initial population...")
    population = createPopulation(initialPopSize, parameterLength)
    allPopulations = []
    allFitnesses = []
    averageFitness = []
    peakFitness = []
    genCount = 0
    # Loop through n generationss
    for i in range(numberOfGenerations):
        genCount += 1
        print()
        print("+++++++ ANALYZING GENERATION: " + str(genCount) + " +++++++")

        # Stpre population
        allPopulations.append(population)
        # Calculate fitness & calculate fitness in one step
        fitnessOutput = selectFittest(topNum, population)
        fittest = fitnessOutput[0]
        fittestFitness = fitnessOutput[1]
        populationFitness = fitnessOutput[2]

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
    return allPopulations, allFitnesses, averageFitness, peakFitness, GA_settings


fitnessDictionary = {}
GAresults = GAfunction(initialPopSize=5, parameterLength=28, numberOfGenerations=2, topNum=2, childrenPerGeneration=5, crossoverPoints=3, breedingVariation=0, mutationRate=0.01, includeFittestParents=2)


allPopulations = GAresults[0]# Get all populations
allFitnesses = GAresults[1] # Get all fitnesses
averageFitness = GAresults[2] # Get averageFitness per generation
peakFitness = GAresults[3] # Get the peakFitness per generation
numberOfGenerations = GAresults[4]['numberOfGenerations'] # Get recording # of generations

from Python_scripts.Extra_X2K_functions import tell_parameters, getFittestIndividual
tell_parameters( getFittestIndividual(GAresults) )


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
copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET1.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET1.txt")
## Run GA with Subset1
fitnessDictionary = {}
GAresults_Subset1 = GAfunction(initialPopSize=100, parameterLength=28, numberOfGenerations=50, topNum=10, childrenPerGeneration=90, crossoverPoints=4, breedingVariation=0, mutationRate=0.01, includeFittestParents=10)



# Put Subset2 GMT file into testgmt folder
## Delete whatever file is in testgmt
dir_name = "data/testgmt/"
files = os.listdir(dir_name)
for item in files:
    if item.endswith(".txt"):
        os.remove(os.path.join(dir_name, item))
## Replace it with subset 1
copyfile("Validation/Perturbation_Data/GEO/Kinase_Perturbations_from_GEO_SUBSET2.txt", "data/testgmt/Kinase_Perturbations_from_GEO_SUBSET2.txt")
## Run GA with Subset2
fitnessDictionary = {} # Make sure to reset fitnessDictionary to capture new fitnesses
allFitnesses_Subset2 = []
averageFitness_Subset2 = []
for generation in GAresults_Subset1[0]:
    new_fitnesses = calculateFitness(generation)
    allFitnesses_Subset2.append(new_fitnesses)
    averageFitness_Subset2.append( sum(new_fitnesses) / len(new_fitnesses) )
peakFitness_Subset2 = []
for generation in allFitnesses_Subset2:
    peakFitness_Subset2.append(sorted(generation, reverse=True)[0])


# Save/load GAresults as file
# Save

GA_output_name = 'GA_results.npy'
import os, numpy as np
results_dir = 'GA_Results/'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
np.save(results_dir+GA_output_name, list([GAresults_Subset1, allFitnesses_Subset2, \
                                          averageFitness_Subset2,peakFitness_Subset2]))
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

