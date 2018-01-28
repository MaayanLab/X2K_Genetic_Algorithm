
# Tells you the parameters for a given binary string
def tell_parameters(binary, verbose=True):
    # Readjust bit parameters if length of binary string changes
    bit_dict = {"TF_sort": binary[0:2],
                "TF_species": binary[2:4],
                "TF_databases": binary[4:6],
                "TF_background": binary[6:8],
                "TF_topTFs": binary[8:10],
                "PPI_databases": binary[10:28],
                "PPI_pathLength": binary[28:29],
                "KINASE_sort": binary[29:31],
                "KINASE_interactome": binary[31:33],
                "KINASE_background": binary[33:35]
                }
    KINASE_topKinases = 20 # Should be consistent for every individual

    # CHEA OPTIONS ############
    TF_sort = {"10":"oddsratio", "01":"combined_score", "11":"rank", "00":"pvalue"}
    TF_species = {"10":"human", "01":"mouse", "11":"both", "00":"RESHUFFLE"}
    TF_databases = {"10":"chea", "01":"transfac", "11":"both", "00":"RESHUFFLE"}
    TF_background = {"10":"humanarchs4", "01":"humanarchs4", "11":"humanarchs4", "00":"RESHUFFLE"}
    TF_topTFs = {"10":5, "01":10, "11":20, "00":"RESHUFFLE"}

    # G2N OPTIONS ############
    ## Databases
    all_PPI_databases = ["BIND", "BIOCARTA", "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEGG",
                         "MINT", "MIPS", "MURPHY", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL",
                         "HUMAP"]  # "BIOGRID", "KEA"
    PPI_dbs = []
    PPI_string = bit_dict["PPI_databases"]
    bit_dict["PPI_databases"] = PPI_string
    ### Generate list of G2N databases
    for ind,bit in enumerate(PPI_string):
        if bit == "1":
            PPI_dbs.append( all_PPI_databases[ind] )
    PPI_databases = ",".join(PPI_dbs)
    ## Path length
    PPI_pathLength = {"0":1, "1":2}


    ############ KEA OPTIONS ############
    KINASE_sort = {"10":"oddsratio", "01":"combined_score", "11":"rank", "00":"pvalue"}
    KINASE_interactome = {"10":"KP", "01":"P", "11":"RESHUFFLE", "00":"RESHUFFLE"}
    KINASE_background = {"10":"humanarchs4", "01":"humanarchs4", "11":"humanarchs4", "00":"RESHUFFLE"}

    parametersList = ["TF_sort","TF_species","TF_databases","TF_background", "TF_topTFs", "PPI_pathLength", "KINASE_sort", "KINASE_interactome", "KINASE_background"]
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
    TF_params = ';'.join( ["run", TF_sort[bit_dict["TF_sort"]], TF_species[bit_dict["TF_species"]], TF_databases[bit_dict["TF_databases"]], TF_background[bit_dict["TF_background"]], str(TF_topTFs[bit_dict["TF_topTFs"]]) ] )
    PPI_params = ';'.join(["run",PPI_databases, str(PPI_pathLength[bit_dict["PPI_pathLength"]]) ])
    KINASE_params = ';'.join(["run", KINASE_sort[bit_dict["KINASE_sort"]], KINASE_background[bit_dict["KINASE_background"]], KINASE_interactome[bit_dict["KINASE_interactome"]], str(KINASE_topKinases)  ])
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


    # Compile parameters into dataframe
    ## Initialize lists
    tf_sort = []
    tf_species = []
    tf_databases = []
    tf_background = []
    tf_topTFs = []
    ppi_databases = []
    ppi_pathLength = []
    kinase_sort = []
    kinase_background = []
    kinase_interactome = []
    kinase_topKinases = []
    for individual in superPopulation:
        params = tell_parameters(individual, verbose=False)
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
        ## KEA
        kinase_params = params[2].split(";")
        kinase_sort.append(kinase_params[1])
        kinase_background.append(kinase_params[2])
        kinase_interactome.append(kinase_params[3])
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
         'KINASE_sort':pd.Categorical(kinase_sort),
         'KINASE_background':pd.Categorical(kinase_background),
         'KINASE_interactome':pd.Categorical(kinase_interactome),
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
def parameterEvolutionPlot(GAresults, figsize=(24,8), chance=4.22,saveFig='No'):
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
    Fitness_peak = pd.DataFrame({'Fitness_peak': pd.Series(GAresults[3], list(range(1, len(pd.Series(GAresults[3])) + 1)))})
    # Fitness plot
    # Calculate stdv for average Fitness
    yerr1 = np.std(GAresults[1], axis=1)
    x = Fitness_avg.index
    ax1 = plt.subplot2grid((nrows, 1), (0, 0), facecolor='whitesmoke', rowspan=fitness_rows)
    plt.gcf().set_facecolor('white')  # Change plot border background color
    #ax1.plot(Fitness_both.index, Fitness_both['Fitness'], color='blue', linestyle='-', marker='.', markersize=2)
    plt.errorbar(x, Fitness_avg['Fitness'], yerr=yerr1, color='blue', marker='o',markersize=5, capsize=2, label=" Average Fitness")
    ax1.plot(x, Fitness_peak['Fitness_peak'], 'purple', linestyle='-', marker='^', markersize=5, label="Peak Fitness")
    ax1.axhline(y=chance, linestyle="--", color='r', label="Chance levels")
    #plt.title('Fitness Over Generations', fontsize=20)
    plt.xlabel('Generation', fontsize=12)
    plt.ylabel('Fitness', fontsize=12)
    plt.tick_params(axis='x', labelsize=12)
    plt.legend(loc='lower right', borderaxespad=2)
    plt.ylim([0, max(data['Fitness'])+5])
    plt.subplots_adjust(right=padRight, left=padLeft)
    plt.legend(loc='center left', bbox_to_anchor=(1, .5), ncol=70)  # bbox_to_anchor=(horizontal, vertical)
    plt.xticks(np.arange(1, max(x) + 1, 1))


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
    plt.yticks(np.arange(0, max(data['Average_PPI_size']), 3500))
    plt.xticks(np.arange(1, max(x) + 1, 1))

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
        plt.xticks(np.arange(1, max(x) + 1, 1))
        plt.subplots_adjust(right=padRight, left=padLeft)  # Expand plot area to get more of legends in view
        if i != param_num-1: # Turn of xtick labels for all but bottom plot
            plt.tick_params(axis='x', labelbottom='off')
            plt.xlabel('')
    if saveFig!='No':
        import os
        if not os.path.exists('Figures/'):
            os.makedirs('Figures/')
        plt.savefig('Figures/'+saveFig+'.eps', format='eps', dpi=1000)

# parameterEvolutionPlot(GAresults_Subset1)
#
# def ParameterBoxplots(GAresults):
#     import pandas as pd
#     data = parameterDF(GAresults)
#     for parameter in data.columns[3:]:
#         pd.DataFrame.boxplot(data, 'Fitness', parameter)
#
# def fitnessHistogramCurves(allFitnesses, gen_spacing=10):
#     from scipy.stats import norm
#     import matplotlib.mlab as mlab
#     import matplotlib.pyplot as plt
#     from matplotlib.pyplot import cm
#     which_gens = [0] + list(range(-1, len(allFitnesses)+1))[0::gen_spacing][1:]
#
#     color = cm.rainbow(np.linspace(0, 1, len(which_gens)))
#
#     for i,gen in enumerate(which_gens):
#         # LINE METHOD
#         datos=allFitnesses[gen-1]
#         # best fit of data
#         (mu, sigma) = norm.fit(datos)
#         # the histogram of the data
#         n, bins, patches = plt.hist(datos, 5, normed=1, facecolor=color[i], alpha=0.75) # Make the histogram bar invisible by coloring white
#         # add a 'best fit' line
#         y = mlab.normpdf(bins, mu, sigma)
#         l = plt.plot(bins, y, color[i], linewidth=1)
#         # plot
#         plt.xlabel('Generation')
#         plt.ylabel('Frequency')
#         plt.title('Distribution of All Fitnesses: Every %.0f Generations' % (gen_spacing))
#         plt.grid(True)
#         plt.show()
#
#     for i gen in enumerate(which_gens):
#         #MORIGINAL METHOD
#         plt.hist(allFitnesses[0], bins=20, color=color[i])  # 1st
#
#




def parameterStats(GAresults, write_Excel=False):
    import pandas as pd
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    data = parameterDF(GAresults)
    frames = []
    tested_parameters = list(data.columns[4:][:-1]) # Test which parameters to run through ANOVA
    tested_parameters = [x for x in tested_parameters if x != 'KINASE_topKinases'] # Remove specific parameter
    tested_parameters = [x for x in tested_parameters if x != 'KINASE_background']  # Remove specific parameter
    tested_parameters = [x for x in tested_parameters if x != 'TF_background']  # Remove specific parameter
    for parameter in tested_parameters:
        # Rearrange data
        #grps = pd.unique(data[parameter].values)
        #d_data = {grp: data['Fitness'][data[parameter] == grp] for grp in grps}

        # One way ANOVA
        mod = ols('Fitness ~ '+parameter,data=data).fit()
        aov_table = sm.stats.anova_lm(mod, typ=1)
        p = aov_table['PR(>F)'][0]
        # Add P-val summary
        if p > 0.05:
            aov_table['Sig'] = "non-sig"
        if p <= 0.05:
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

