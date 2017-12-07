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
from X2K_Pipeline import X2K_fitness


fitnessDictionary = {}

def calculateFitness(population):
    indCount = 1
    fitness = []

    for i in range(len(population)):
        # Delete .DS_Store files
        import os
        print()
        print("(Deleting .DS_Store files...)")
        for root, dirs, files in os.walk(os.getcwd()):
            for file in files:
                if file.endswith(".DS_Store"):
                    print(os.path.join(root, file))
                    os.remove(os.path.join(root, file))
        # Calculate fitness (ONLY if it hasn't been previously calculated)
        print()
        print("****** Individual " + str(indCount) + " ******")
        if population[i] not in fitnessDictionary:
            # fitness.append( sum(map(int, population[i])) )
            new_fitness = X2K_fitness(population[i])
            fitness.append( new_fitness )
            fitnessDictionary[population[i]] = new_fitness
        else:
            fitness.append( fitnessDictionary[population[i]] )
            print("Using stored fitness")
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
# Run Genetic Algorithm
###################################

def GAfunction(initialPopSize, parameterLength, numberOfGenerations, topNum, childrenPerGeneration, crossoverPoints, breedingVariation, mutationRate, includeFittestParents):
    print("Creating initial population...")
    population = createPopulation(initialPopSize, parameterLength)
    allPopulations = []
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

        # Store average fitness
        averageFitness.append( sum(populationFitness)*1.0 / len(populationFitness) )
        # Store peak fitness
        peakFitness.append( max(populationFitness) )
    return allPopulations, averageFitness, peakFitness, numberOfGenerations

GAresults = GAfunction(initialPopSize=100, parameterLength=27, numberOfGenerations=100, topNum=4, childrenPerGeneration=8, crossoverPoints=3, breedingVariation=0, mutationRate=0.01, includeFittestParents=2)

finalPopulation = GAresults[0][0] # Get final population
averageFitness = GAresults[1] # Get averageFitness per generation
peakFitness = GAresults[2] # Get the peakFitness per generation
numberOfGenerations = GAresults[3] # Get recording # of generations

###################################
# 7. Plot results
###################################
import matplotlib.pyplot as plt
# Plot averageFitness
y = averageFitness
x = range(len(averageFitness))
plt.plot(x, y, 'bo--')
plt.xlabel('Fitness')
plt.ylabel('Generation')
plt.show()
plt.title('Average Fitness Over Generations')
# axes = plt.gca()
# axes.set_ylim([0,30])

# Plot peakFitness
y = peakFitness
x = range(len(peakFitness))
plt.plot(x, y, 'mo--')
plt.xlabel('Fitness')
plt.ylabel('Generation')
plt.show()
plt.title('Top Fitness Over Generations')
# axes = plt.gca()
# axes.set_ylim([0,30])
