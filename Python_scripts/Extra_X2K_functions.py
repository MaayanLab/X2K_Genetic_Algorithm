# Tells you the parameters for a given binary string
def tell_parameters(binary, verbose=True):
    # Readjust bit parameters if length of binary string changes
    bit_dict = {"CHEA_sort":binary[0:2],
                "CHEA_species":binary[2:4],
                "CHEA_databases":binary[4:6],
                "CHEA_background":binary[6:8],
                "CHEA_topTFs":binary[8:10],
                "G2N_databases":binary[10:20],
                "G2N_pathLength":binary[20:21],
                "KEA_sort":binary[21:23],
                "KEA_interactome":binary[23:25],
                "KEA_background":binary[25:27]
                }
    KEA_topKinases = 10 # Should be consistent for every individual

    # CHEA OPTIONS ############
    ## sort
    CHEA_sort = {"10":"oddsratio", "01":"combined_score", "11":"rank", "00":"RESHUFFLE"}
    CHEA_species = {"10":"human", "01":"mouse", "11":"both", "00":"RESHUFFLE"}
    CHEA_databases = {"10":"chea", "01":"transfac", "11":"both", "00":"RESHUFFLE"}
    CHEA_background = {"10":"humanarchs4", "01":"mousearchs4", "11":"hg133", "00":"RESHUFFLE"}
    CHEA_topTFs = {"10":5, "01":10, "11":20, "00":"RESHUFFLE"}

    # G2N OPTIONS ############
    ## Databases
    all_G2N_databases = ["BIND","BIOCARTA","BIOGRID","DIP","FIGEYS","HPRD","INNATEDB","INTACT","KEA","KEGG","MINT","MIPS","MURPHY","PDZBASE","PPID","PREDICTEDPPI","SNAVI","STELZL","VIDAL","HUMAP"]
    G2N_dbs = []
    G2N_string = bit_dict["G2N_databases"]
    ### If no database selected, select at least one
    from random import randint
    if sum(map(int, G2N_string)) == 0:
        g2n_list = list(G2N_string)
        g2n_list[randint(0, (len(g2n_list)) - 1)] = '1'
        bit_dict["G2N_databases"] = "".join(g2n_list)
    ### Generate list of G2N databases
    for ind,bit in enumerate(G2N_string):
        if bit == "1":
            G2N_dbs.append( all_G2N_databases[ind] )
    G2N_databases = ",".join(G2N_dbs)
    ## Path length
    G2N_pathLength = {"0":1, "1":2}


    ############ KEA OPTIONS ############
    KEA_sort = {"10":"oddsratio", "01":"combined_score", "11":"rank", "00":"RESHUFFLE"}
    KEA_interactome = {"10":"KP", "01":"P", "11":"RESHUFFLE", "00":"RESHUFFLE"}
    KEA_background = {"10":"humanarchs4", "01":"mousearchs4", "11":"hg133", "00":"RESHUFFLE"}

    parametersList = ["CHEA_sort","CHEA_species","CHEA_databases","CHEA_background", "CHEA_topTFs", "G2N_pathLength", "KEA_sort", "KEA_interactome", "KEA_background"] #g2n_databases
    from random import choice
    # RESHUFFLE any bits that aren't legitimate options
    for param in parametersList:
        bad_bits = [k for k, v in eval(param).items() if v == "RESHUFFLE"]
        good_bits = [k for k, v in eval(param).items() if v != "RESHUFFLE"]
        current_bits = bit_dict[param]
        if current_bits in bad_bits:
            bit_dict[param] = choice( good_bits )
            #print("Reshuffling...")
    # Tell X2K parameters
    chea_params = ';'.join( ["run", CHEA_sort[bit_dict["CHEA_sort"]], CHEA_species[bit_dict["CHEA_species"]], CHEA_databases[bit_dict["CHEA_databases"]], CHEA_background[bit_dict["CHEA_background"]], str(CHEA_topTFs[bit_dict["CHEA_topTFs"]]) ] )
    g2n_params = ';'.join(["run",G2N_databases, str(G2N_pathLength[bit_dict["G2N_pathLength"]]) ])
    kea_params = ';'.join(["run", KEA_sort[bit_dict["KEA_sort"]], KEA_background[bit_dict["KEA_background"]], KEA_interactome[bit_dict["KEA_interactome"]], str(KEA_topKinases)  ])
    if(verbose==True):
        print()
        print( "___CHEA Parameters___")
        print(chea_params)
        print()
        print( "___G2N Parameters___")
        print(g2n_params)
        print()
        print( "___KEA Parameters___")
        print(kea_params)
        print()

    return chea_params, g2n_params, kea_params



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
    generation = 0
    for sublist in GAresults[1]:
        for item in sublist:
            superFitnesses.append(item)
            Generation.append(generation)
        generation += 1

    # Compile parameters into dataframe
    ## Initialize lists
    chea_sort = []
    chea_species = []
    chea_databases = []
    chea_background = []
    chea_topTFs = []
    g2n_databases = []
    g2n_pathLength = []
    kea_sort = []
    kea_background = []
    kea_interactome = []
    kea_topKinases = []
    for individual in superPopulation:
        params = tell_parameters(individual, verbose=False)
        ## CHEA
        chea_params = params[0].split(";")
        chea_sort.append(chea_params[1])
        chea_species.append(chea_params[2])
        chea_databases.append(chea_params[3])
        chea_background.append(chea_params[4])
        chea_topTFs.append(chea_params[5])
        ## G2N
        g2n_params = params[1].split(";")
        g2n_databases.append(g2n_params[1])
        g2n_pathLength.append(g2n_params[2])
        ## KEA
        kea_params = params[2].split(";")
        kea_sort.append(kea_params[1])
        kea_background.append(kea_params[2])
        kea_interactome.append(kea_params[3])
        kea_topKinases.append(kea_params[4])
    # Turn lists into dictionary
    import pandas as pd
    #ind = list(range(1,len(superPopulation)+1))
    d = {'Binary':superPopulation,
         'Fitness':pd.Series(superFitnesses),
         'Generation':pd.Series(Generation),
         'chea_sort' :pd.Categorical(chea_sort),
         'chea_species':pd.Categorical(chea_species),
         'chea_databases':pd.Categorical(chea_databases),
         'chea_background':pd.Categorical(chea_background),
         'chea_topTFs':pd.Categorical(chea_topTFs),
         'g2n_databases':pd.Categorical(g2n_databases),
         'g2n_pathLength':pd.Categorical(g2n_pathLength),
         'kea_sort':pd.Categorical(kea_sort),
         'kea_background':pd.Categorical(kea_background),
         'kea_interactome':pd.Categorical(kea_interactome),
         'kea_topKinases':pd.Categorical(kea_topKinases)}
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


import numpy as np
GA = np.load("GA_Results/GA_results.npy")
GAresults = GA[0]

# Make super awesome plot showing evolution of parameter distribution over generations/time
def parameterEvolutionPlot(GAresults, figsize=(24,8)):
    data = parameterDF(GAresults)
    import matplotlib.pyplot as plt
    import pandas as pd
    # Setup Fitness data
    Fitness_avg = data[['Fitness', 'Generation']].groupby('Generation').mean()
    PeakFits = pd.DataFrame({'Fitness_peak': pd.Series(GAresults[3], list(range(1, len(pd.Series(GAresults[3])) + 1)))})
    Fitness_both = pd.concat([Fitness_avg, PeakFits], axis=1)
    # Setup dimension vars
    parameters = list(data.columns[3:])
    param_num = len(parameters)
    fitness_rows = 4
    nrows = param_num + fitness_rows

    # Fitness plot
    ax1 = plt.subplot2grid((nrows, 1), (0, fitness_rows), facecolor='whitesmoke')
    plt.gcf().set_facecolor('white')  # Change plot border background color
    ax1.plot(Fitness_both.index, Fitness_both['Fitness'], 'b|-', Fitness_both.index, Fitness_both['Fitness_peak'], 'g.--')

    plt.title('Fitness Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend(('Average Fitness', 'Peak Fitness'), loc='lower right')
    plt.ylim([0, 50])
    plt.xlim([1, 50])
    plt.subplots_adjust(right=0.5, left=0.1)
    # Parameter plots
    import seaborn as sns
    for i,parameter in enumerate(parameters):
        row_num = i+3 # Skip the first 3 rows since the other plot is occupying that space
        tab1 = freqTable(data, parameter)
        # Set palette
        sns.set_palette(sns.color_palette("BuPu", len(tab1.columns)))  # Run this BEFORE you create your subplot
        # Turn ax string into variable
        new_ax = 'ax' + str(row_num)
        exec(new_ax + " = plt.subplot2grid((nrows, 1), (row_num, 0))")
        # Run ax variable from string
        tab1.plot(kind='area', ax=eval(new_ax), stacked=True, figsize=figsize)  # , colors=['hotpink','slateblue','navy']
        plt.legend(loc='center left', bbox_to_anchor=(1, .5), ncol=70)  # bbox_to_anchor=(horizontal, vertical)
        plt.xticks(rotation=0)
        plt.xlabel('')
        plt.ylabel(parameter, rotation=0, labelpad=50)
        plt.xlim([1, 50])
        plt.subplots_adjust(right=0.5, left=0.1)  # Expand plot area to get more of legends in view
        if i != 10: # Turn of xtick labels for all but bottom plot
            plt.tick_params(axis='x', labelbottom='off')

# parameterEvolutionPlot(GAresults)

def ParameterBoxplots(GAresults):
    import pandas as pd
    data = parameterDF(GAresults)
    for parameter in data.columns[3:]:
        pd.DataFrame.boxplot(data, 'Fitness', parameter)


def parameterStats(GAresults, write_Excel=False):
    import pandas as pd
    data = parameterDF(GAresults)
    frames = []
    tested_parameters = list(data.columns[3:][:-1]) # Test which parameters to run through ANOVA
    for parameter in tested_parameters:
        # Rearrange data
        grps = pd.unique(data[parameter].values)
        d_data = {grp: data['Fitness'][data[parameter] == grp] for grp in grps}
        # One way ANOVA
        import statsmodels.api as sm
        from statsmodels.formula.api import ols

        mod = ols('Fitness ~ '+parameter,data=data).fit()
        aov_table = sm.stats.anova_lm(mod, typ=2)
        p = aov_table['PR(>F)'][0]
        # Add P-val summary
        if p >= 0.05:
            aov_table['Sig'] = "non-sig"
        if p < 0.05:
            aov_table['Sig'] = '*'
        if p < 0.01:
            aov_table['Sig'] = '**'
        if p < 0.001:
            aov_table['Sig'] ='***'
        if p < 0.0001:
            aov_table['Sig'] = '****'
        # mod.summary() # For full summary
        frames.append(aov_table)
        parameter_AOV_results = pd.concat(frames)
    print(parameter_AOV_results)
    if write_Excel == True:
        # Write AOV results to excel file
        print("Writing AOV results to excel file...")
        writer = pd.ExcelWriter('GA_Results/Parameter_AOV_Results.xlsx')
        parameter_AOV_results.to_excel(writer, 'Sheet1')
        writer.save()
    return parameter_AOV_results

# parameter_AOV_results = parameterStats(GAresults, True)

