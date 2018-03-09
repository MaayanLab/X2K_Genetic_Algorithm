# binaryString = '0110110100001010110101101101100000001011000'
# GA_Results = GAresults_Subset1
# import numpy as np
# GA_output_name = 'GA_results_GEO.wPPIlimiters.npy'
# results_file = 'X2K_Genetic_Algorithm/GA_Results/GEO/'+GA_output_name
# GAresults = np.load(results_file)[0]
#

# Tells you the parameters for a given binary string
def tell_parameters(binaryString, verbose=True):
    constantBackground = "humanarchs4"

    # CHEA OPTIONS ############
    if constantBackground == False:
        TF_background = {"10": "humanarchs4", "01": "mousearchs4", "11": "hg133", "00": "RESHUFFLE"}
    else:
        TF_background = {"10": constantBackground, "01": constantBackground, "11": constantBackground,
                         "00": constantBackground}
    TF_sort = {"10": "oddsratio", "01": "combined_score", "11": "rank", "00": "pvalue"}
    TF_species = {"10": "human", "01": "mouse", "11": "both", "00": "RESHUFFLE"}
    TF_topTFs = {"10": 5, "01": 10, "11": 20, "00": "RESHUFFLE"}
    # TF_databases = ["ARCHS4", "CHEA", "CREEDS", "ENCODE", "ENRICHR", "HUMAP", "IREF", "JASPAR-TRANSFAC", "TF-PPI","TF-LOF"]
    TF_databases = {"10": "chea", "01": "transfac", "11": "both", "00": "RESHUFFLE"}

    # G2N OPTIONS ############
    PPI_pathLength = {"0": 1, "1": 1}
    PPI_minLength = {"10": 1, "01": 5, "11": 10, "00": 20}
    PPI_maxLength = {"10": 50, "01": 100, "11": 200, "00": 500}
    PPI_finalSize = {"10": 50, "01": 100, "11": 200, "00": 500}
    PPI_databases = ["BIND", "BIOGRID", "BIOCARTA", "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEGG", \
                     "MINT", "MIPS", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL", "HUMAP", "IREF",
                     "BIOPLEX"]  # "KEA", "MURPHY", **** ,
    # KEA OPTIONS ############
    if constantBackground == False:
        KINASE_background = {"10": "humanarchs4", "01": "mousearchs4", "11": "hg133", "00": "RESHUFFLE"}
    else:
        KINASE_background = {"10": constantBackground, "01": constantBackground, "11": constantBackground,
                             "00": constantBackground}
    KINASE_topKinases = 20  # Should be consistent for every individual
    KINASE_sort = {"10": "oddsratio", "01": "combined_score", "11": "rank", "00": "pvalue"}
    # KINASE_databases = ["CREEDS", "iPTMnet", "iREF", "MINT", "HPRD", "PHOSPHOELM", "PHOSPHOPOINT", "PHOSPHOSITE",  "PHOSPHOSITEPLUS", "NETWORKIN"]
    KINASE_databases = {"10": "KP", "01": "P", "11": "RESHUFFLE", "00": "RESHUFFLE"}

    # Flexibly create a dictionary for each bit value and its meaning depending on the particular set of variables
    parameterList = ["TF_sort", "TF_species", "TF_databases", "TF_background", "TF_topTFs", "PPI_databases",
                     "PPI_pathLength", "PPI_minLength", "PPI_maxLength", "PPI_finalSize", "KINASE_sort",
                     "KINASE_databases", "KINASE_background"]
    # def constructBitDict(parameterList, binaryString):
    bit_dict = {}
    prevBits = 0
    for param in parameterList:
        realParam = eval(param)
        # Get the number of necessary bits for the specific parameter
        if type(realParam) == dict:
            numBits = len(next(iter(realParam.keys())))
        elif type(realParam) == list:
            numBits = len(realParam)
        # Assign positions of first and last bits within binaryString
        firstBit = prevBits
        lastBit = firstBit + numBits
        bit_dict[param] = binaryString[firstBit:lastBit]
        prevBits = prevBits + numBits
        # return bitDict

    # bit_dict = constructBitDict(parameterList, binaryString)

    def totalBitLength():
        totalLength = 0
        for string in bit_dict.values():
            totalLength = totalLength + len(string)
        print(totalLength)

    # totalBitLength()


    # General function for converting binary into list of selected databases
    def selectDatabases(databaseType, fullDatabaseList):
        # ******* Turn off specific databases *******
        # ppi_list = list(PPI_string)
        # ppi_list[2] = '0'  # Turn off BIOGRID PPI
        # ppi_list[8] = '0'  # Turn off KEA PPI
        # PPI_string = "".join(ppi_list)
        # bit_dict["PPI_databases"] = PPI_string
        # *******************************************
        ### If no database selected, select at least one
        dbs = []
        binaryString = bit_dict[databaseType + "_databases"]
        from random import randint
        while sum(map(int, binaryString)) == 0:
            binaryList = list(binaryString)
            binaryList[randint(0, (len(binaryList)) - 1)] = '1'
            # ******* Turn off specific databases ******* AGAIN just in case it picked one of these
            # ppi_list[2] = '0'  # Turn off BIOGRID PPI
            # ppi_list[8] = '0'  # Turn off KEA PPI
            binaryString = "".join(binaryList)
            bit_dict[databaseType + "_databases"] = binaryString
        ### Generate list of G2N databases
        for ind, bit in enumerate(bit_dict[databaseType + "_databases"]):
            if bit == "1":
                dbs.append(fullDatabaseList[ind])
        selectedDatabases = ",".join(dbs)
        return selectedDatabases

    # ############ CHEA Databases ############
    # selectedTFdatabases = selectDatabases("TF", TF_databases)
    # ############ G2N Databases ############
    selectedPPIdatabases = selectDatabases("PPI", PPI_databases)
    # ############ KEA Databases ############
    # selectedKINASEdatabases = selectDatabases("KINASE", KINASE_databases)


    # Tell X2K parameters
    TF_params = ';'.join( ["run", TF_sort[bit_dict["TF_sort"]], TF_species[bit_dict["TF_species"]], TF_databases[bit_dict["TF_databases"]], TF_background[bit_dict["TF_background"]], str(TF_topTFs[bit_dict["TF_topTFs"]]) ] )
    PPI_params = ';'.join(["run", selectedPPIdatabases, str(PPI_pathLength[bit_dict["PPI_pathLength"]]), str(PPI_minLength[bit_dict["PPI_minLength"]]), str(PPI_maxLength[bit_dict["PPI_maxLength"]]), str(PPI_finalSize[bit_dict["PPI_finalSize"]]) ])
    KINASE_params = ';'.join(["run", KINASE_sort[bit_dict["KINASE_sort"]], KINASE_background[bit_dict["KINASE_background"]], KINASE_databases[bit_dict["KINASE_databases"]], str(KINASE_topKinases)  ])
    if(verbose==True):
        print()
        print( "___TF (CHEA) Parameters___")
        print(TF_params)
        print()
        print( "___PPI (G2N) Parameters___")
        print(PPI_params)
        print()
        print( "___KINASE (KEA) Parameters___")
        print(KINASE_params)
        print()

    return TF_params, PPI_params, KINASE_params



def getFittestIndividual(GAresults):
    finalPop = GAresults[0][-1] # Final pop binaries
    finalFit = GAresults[1][-1] # Final pop fitnesses
    peak = max(GAresults[3]) # Final peak fitness
    fittestBinary = finalPop[ finalFit.index(peak) ]
    return fittestBinary

# def getPeakParameters(GAresults):
#     #tell_parameters(getFittestIndividual(GAresults))



# Create a dataframe the fitness and parameters of every individual from every generation
def parameterDF(GAresults):
    # Get all individuals in one long list
    superPopulation = []
    for sublist in GAresults[0]:
        for item in sublist:
            superPopulation.append(item)
    # Get all fitness scores in one long list
    superFitnesses = []
    Generation = []
    generation = 1
    for sublist in GAresults[1]:
        for item in sublist:
            superFitnesses.append(item)
            Generation.append(generation)
        generation += 1
    # Get all average_PPI_sizes
    average_PPI_size = []
    for sublist in GAresults[5]:
        for item in sublist:
            average_PPI_size.append(item)
    superParams=[]
    for sublist in GAresults[7]:
        for item in sublist:
            superParams.append(item)

    # Compile parameters into dataframe
    ## Initialize lists
    tf_sort = []
    tf_species = []
    tf_databases = []
    tf_background = []
    tf_topTFs = []

    ppi_databases = []
    ppi_pathLength = []
    ppi_minLength = []
    ppi_maxLength = []
    ppi_finalSize = []

    kinase_sort = []
    kinase_background = []
    kinase_databases = []
    kinase_topKinases = []
    #for i in superPopulation:
        #params = tell_parameters(individual, verbose=False)
    for params in superParams:
        ## CHEA
        tf_params = params[0].split(";")
        tf_sort.append(tf_params[1])
        tf_species.append(tf_params[2])
        tf_databases.append(tf_params[3])
        tf_background.append(tf_params[4])
        tf_topTFs.append(tf_params[5])
        ## G2N
        ppi_params = params[1].split(";")
        ppi_databases.append(ppi_params[1])
        ppi_pathLength.append(ppi_params[2])
        ppi_minLength.append(ppi_params[3])
        ppi_maxLength.append(ppi_params[4])
        ppi_finalSize.append(ppi_params[5])
        ## KEA
        kinase_params = params[2].split(";")
        kinase_sort.append(kinase_params[1])
        kinase_background.append(kinase_params[2])
        kinase_databases.append(kinase_params[3])
        kinase_topKinases.append(kinase_params[4])
    # Turn lists into dictionary
    import pandas as pd
    #ind = list(range(1,len(superPopulation)+1))
    d = {'Binary':superPopulation,
         'Fitness':pd.Series(superFitnesses),
         'Generation':pd.Series(Generation),
         'TF_sort' :pd.Categorical(tf_sort),
         'TF_species':pd.Categorical(tf_species),
         'TF_databases':pd.Categorical(tf_databases),
         'TF_background':pd.Categorical(tf_background),
         'TF_topTFs':pd.Categorical(tf_topTFs),
         'PPI_databases':pd.Categorical(ppi_databases),
         'PPI_pathLength':pd.Categorical(ppi_pathLength),
         'PPI_minLength':pd.Categorical(ppi_minLength),
         'PPI_maxLength':pd.Categorical(ppi_maxLength),
         'PPI_finalSize':pd.Categorical(ppi_finalSize),
         'KINASE_sort':pd.Categorical(kinase_sort),
         'KINASE_background':pd.Categorical(kinase_background),
         'KINASE_databases':pd.Categorical(kinase_databases),
         'KINASE_topKinases':pd.Categorical(kinase_topKinases),
         'Average_PPI_size':pd.Series(average_PPI_size)
         }
    # Turn dictionary into dataframe
    df = pd.DataFrame(d)
    #print(df[:5]) # preview
    return df
# data = parameterDF(GAresults)



# Create frequency table of parameters for each generation
def freqTable(data, parameter):
    import pandas as pd
    tab = pd.crosstab(index=data['Generation'], columns=data[parameter])
    #t_tab = np.transpose(tab[gen-1:gen])
    return tab


# import numpy as np
# GA = np.load("GA_Results/GA_results.100pop.50gen.3sortOptions.npy")
# GAresults = GA[0]

# Make super awesome plot showing evolution of parameter distribution over generations/time
def parameterEvolutionPlot(GAresults, figsize=(24,8), chance=4.22):
    data = parameterDF(GAresults)
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    # Setup dimension vars
    parameters = sorted(list(data.columns[4:]), reverse=True)  # Flip around to make plot same order as pipeline E
    # Remove specific parameters
    parameters = [x for x in parameters if x != 'TF_background'] # Remove specific parameter
    parameters = [x for x in parameters if x != 'KINASE_background']  # Remove specific parameter
    parameters = [x for x in parameters if x != 'KINASE_topKinases']  # Remove specific parameter
    param_num = len(parameters)
    fitness_rows = 5
    PPI_size_rows = 1
    nrows = fitness_rows+PPI_size_rows+param_num+2
    padRight=.6; padLeft=.2
    # Setup Fitness data
    Fitness_avg = data[['Fitness', 'Generation']].groupby('Generation').mean()
    Fitness_peak = data[['Fitness', 'Generation']].groupby('Generation').max()
    #Fitness_peak = pd.DataFrame({'Fitness_peak': pd.Series(GAresults[3], list(range(1, len(pd.Series(GAresults[3])) + 1)))})
    # Fitness plot
    # Calculate study for average Fitness
    yerr1 = np.std(GAresults[1], axis=1)
    x = Fitness_avg.index
    ax1 = plt.subplot2grid((nrows, 1), (0, 0), facecolor='whitesmoke', rowspan=fitness_rows)
    plt.gcf().set_facecolor('white')  # Change plot border background color
    #ax1.plot(Fitness_both.index, Fitness_both['Fitness'], color='blue', linestyle='-', marker='.', markersize=2)
    plt.errorbar(x, Fitness_avg['Fitness'], yerr=yerr1, color='blue', marker='o',markersize=5, capsize=2, label=" Average Fitness")
    ax1.plot(x, Fitness_peak['Fitness'], 'purple', linestyle='-', marker='^', markersize=5, label="Peak Fitness")
    #ax1.axhline(y=chance, linestyle="--", color='r', label="Chance levels")
    #plt.title('Fitness Over Generations', fontsize=20)
    plt.xlabel('Generation', fontsize=12)
    plt.ylabel('Fitness', fontsize=12)
    plt.tick_params(axis='x', labelsize=12)
    plt.legend(loc='lower left', borderaxespad=2)
    #plt.ylim([0, max(data['Fitness'])+5])
    plt.subplots_adjust(right=padRight, left=padLeft)
    plt.legend(loc='center left', bbox_to_anchor=(1, .5), ncol=1)  # bbox_to_anchor=(horizontal, vertical)
    plt.xticks(np.arange(0, max(x) + 1, 5))

    # Setup PPI size data:
    ## Calculate the population average PPI of individual average PPIs
    PPI_std = np.std(GAresults[5],axis=1,dtype=float)
    PPI_dat = data[['Average_PPI_size', 'Generation']].groupby('Generation').mean()
    # PPI Size plot
    ax2 = plt.subplot2grid((nrows, 1), ((fitness_rows+PPI_size_rows), 0), facecolor='whitesmoke', rowspan=PPI_size_rows+1)
    ax2.plot(PPI_dat.index, PPI_dat['Average_PPI_size'], 'firebrick', marker='d')
    ax2.errorbar(x=PPI_dat.index, y=PPI_dat['Average_PPI_size'], yerr=PPI_std, capsize=3, color='mediumpurple')
    plt.ylabel('PPI Size', rotation=0, labelpad=50, fontsize=12)
    plt.tick_params(axis='x', labelbottom='off')
    plt.tick_params(axis='y', labelsize=7)
    #plt.yticks(np.arange(0, max(data['Average_PPI_size']), 3500))
    plt.xticks(np.arange(0, max(x) + 1, 5))

    # Parameter plots
    import seaborn as sns
    for i,parameter in enumerate(parameters):
        row_num = i+(fitness_rows+PPI_size_rows)+2 # Skip the first 3 rows since the other plot is occupying that space
        tab1 = freqTable(data, parameter)
        # Set palette
        sns.set_palette(sns.color_palette("BuPu", len(tab1.columns)))  # Run this BEFORE you create your subplot
        # Turn ax string into variable
        new_ax = 'ax' + str(row_num)
        exec(new_ax + " = plt.subplot2grid((nrows, 1), (row_num, 0))")
        # Run ax variable from string
        tab1.plot(kind='area', ax=eval(new_ax), stacked=True, figsize=figsize)  # , colors=['hotpink','slateblue','navy']
        plt.legend(loc='center left', bbox_to_anchor=(1, .5), ncol=70)  # bbox_to_anchor=(horizontal, vertical)
        if parameter == 'PPI_databases':
            bestPPIdbs = tell_parameters(getFittestIndividual(GAresults))[1].split(";")[1]
            plt.legend(["Optimized database combination (of "+str(len(tab1.columns))+"):\n"+bestPPIdbs],loc='center left', bbox_to_anchor=(1, .5), ncol=70)
        plt.xticks(rotation=0)
        plt.tick_params(axis='y', labelsize=7)
        plt.tick_params(axis='x', labelsize=12)
        plt.xlabel('Generation',fontsize=12)
        plt.ylabel(parameter, rotation=0, labelpad=50, fontsize=12)
        plt.xticks(np.arange(0, max(x) + 1, 5))
        plt.subplots_adjust(right=padRight, left=padLeft)  # Expand plot area to get more of legends in view
        if i != param_num-1: # Turn of xtick labels for all but bottom plot
            plt.tick_params(axis='x', labelbottom='off')
            plt.xlabel('')

# parameterEvolutionPlot(GAresults)


def ParameterViolinPlots(GAresults, numRows=2, numCols=5, figSize=(10,6)):
    import seaborn as sns
    import matplotlib.pyplot as plt
    df = parameterDF(GAresults)
    tested_parameters = [x for x in df.columns if x not in ['Average_PPI_size','Binary','KINASE_topKinases','KINASE_background','PPI_databases','TF_background']]  # Remove specific parameter
    data =  df[[col for col in df.columns if col in tested_parameters]] # Filter database by revelvant paramters
    plt.figure(figsize=figSize)
    for i,parameter in enumerate(reversed(data.columns[1:])):
        # Turn ax string into variable
        #sns.set_palette(sns.color_palette("BuPu", len(data[parameter].unique()) ))  # Run this BEFORE you create your subplot
        ax = plt.subplot(numRows, numCols, i+1)
        #pd.DataFrame.boxplot(data, column='Fitness', by=parameter, ax=ax, grid=False, notch=True, patch_artist=True, rot=45)
        sns.violinplot(x=parameter, y="Fitness", data=data, palette="BuPu", ax=ax)
        plt.xticks(rotation=45, ha='right')
        plt.title(parameter)
        plt.ylabel(''); plt.xlabel('')
        if i == 0 or i == 3:
            plt.ylabel('Fitness')
    plt.subplots_adjust(hspace=.75, wspace=.5, bottom=.2)
    plt.gcf().set_facecolor('white')





def fitnessHistogramCurves(allFitnesses, genSpacing=2, figSize=(10,6)):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    import seaborn as sns
    import numpy as np
    #whichGens = [0] + list(range(-1, len(allFitnesses)+1))[0::gen_spacing][1:]
    count = 0
    whichGens = [0]
    while count <= len(allFitnesses)-genSpacing-1:
        whichGens.append(count+genSpacing)
        count += genSpacing

    color = cm.rainbow(np.linspace(0, 1, len(whichGens)))

    plt.figure()
    for i,gen in enumerate(whichGens):
        # LINE METHOD
        datos=allFitnesses[i]
        # Use kernal density plot function from seaborn
        sns.kdeplot(datos, shade=True, color=color[i], bw=1, label="Gen: "+str(gen+1), gridsize=100)
    plt.xlabel('Fitness')
    plt.ylabel('Frequency')
    #plt.title('Distribution of All Fitnesses: Every %.0f Generation(s)' % (gen_spacing))
    plt.gcf().set_facecolor('white')
    return plt
#fitnessHistogramCurves(allFitnesses, gen_spacing=1)





def parameterStats(GAresults, writeExcel='No'):
    import pandas as pd
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    import numpy as np
    data = parameterDF(GAresults)
    newTable = pd.DataFrame()
    tested_parameters = list(data.columns[4:][:-1]) # Test which parameters to run through ANOVA
    tested_parameters = [x for x in tested_parameters if x not in ['KINASE_topKinases','KINASE_background','TF_background','PPI_pathLength','PPI_minLength','PPI_maxLength', 'PPI_finalSize']]

    for parameter in reversed(tested_parameters):
        # Rearrange data
        #grps = pd.unique(data[parameter].values)
        #d_data = {grp: data['Fitness'][data[parameter] == grp] for grp in grps}

        # One way ANOVA
        mod = ols('Fitness ~ '+parameter,data=data).fit()
        aov_table = sm.stats.anova_lm(mod, typ=2)
        esq_sm = aov_table['sum_sq'][0] / (aov_table['sum_sq'][0] + aov_table['sum_sq'][1])
        p = aov_table['PR(>F)'][0]
        # Add P-val summary
        if p >= 0.05:
            aov_table['Sig'] = "non-sig"
            P = "≥ 0.05"
        elif p < 0.05:
            aov_table['Sig'] = '*'
            P = "< 0.05"
        elif p < 0.01:
            aov_table['Sig'] = '**'
            P = "< 0.01"
        elif p < 0.001:
            aov_table['Sig'] ='***'
            P = "< 0.001"
        elif p < 0.0001:
            aov_table['Sig'] = '****'
            P = "< 0.0001"

        # mod.summary() # For full summary
        SS = int(round(float(aov_table.iloc[0].sum_sq), 0))
        Sig= aov_table['Sig'][0]
        F = int(round(float(aov_table.iloc[0].F),0))
        df = int(aov_table.df[0])
        newTable = newTable.append( pd.DataFrame(np.column_stack([parameter, SS, df, F, P, Sig]), columns=["X2K Parameter","SS","DF","F-value","P-value","Sig."]))

        #newTable = newTable.append(aov_table.iloc[0])
        ## Using individual numbers
        #frames.append(pd.Series([parameter, esq_sm, p]))
    # Turn rownames into first col
    #.index.name = 'Parameter'
    # df.reset_index(inplace=True)
    # parameter_AOV_results.columns = ['Parameter','SS','df','F','p-value','Sig.']
    # parameter_AOV_results[['SS', 'F', 'p-value','Sig.']].apply(pd.to_numeric)
    # parameter_AOV_results = pd.concat(frames).round(3)
    # parameter_AOV_results = parameter_AOV_results.fillna("")
    # parameter_AOV_results['df'] = pd.to_parameter_AOV_results['df'].to_numeric()
    # print(parameter_AOV_results)
    # round(parameter_AOV_results)

    if writeExcel != 'No':
        # Write AOV results to excel file
        print("Writing AOV results to excel file...")
        writer = pd.ExcelWriter(writeExcel)
        newTable.to_excel(writer, 'Sheet1')
        writer.save()

    return newTable

# parameter_AOV_results = parameterStats(GAresults, True)


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