

def parameterDF(GA_df):
    import pandas as pd
    df = GA_df.copy()
    df["TF_sort"] = df["CHEA_parameters"].str.split(";").apply(lambda x: x[1])
    df["TF_species"] = df["CHEA_parameters"].str.split(";").apply(lambda x: x[2])
    df["TF_databases"] = df["CHEA_parameters"].str.split(";").apply(lambda x: x[3])
    df["TF_background"] = df["CHEA_parameters"].str.split(";").apply(lambda x: x[4])
    df["TF_topTFs"] = df["CHEA_parameters"].str.split(";").apply(lambda x: x[5])

    df["PPI_databases"] = df["G2N_parameters"].str.split(";").apply(lambda x: x[1])
    df["PPI_pathLength"] = df["G2N_parameters"].str.split(";").apply(lambda x: x[2])
    df["PPI_minLength"] = df["G2N_parameters"].str.split(";").apply(lambda x: x[3])
    df["PPI_maxLength"] = df["G2N_parameters"].str.split(";").apply(lambda x: x[4])
    df["PPI_finalSize"] = df["G2N_parameters"].str.split(";").apply(lambda x: x[5])

    df["KINASE_sort"] = df["KEA_parameters"].str.split(";").apply(lambda x: x[1])
    df["KINASE_background"] = df["KEA_parameters"].str.split(";").apply(lambda x: x[2])
    df["KINASE_databases"] = df["KEA_parameters"].str.split(";").apply(lambda x: x[3])
    df["KINASE_topKinases"] = df["KEA_parameters"].str.split(";").apply(lambda x: x[4])

    df['Generation'] = pd.to_numeric(df['Generation'])
    df['uniqueID'] = range(0,len(df))
    return df
# data = parameterDF(GA_df)

# Create frequency table of parameters for each generation
def freqTable(data, parameter):
    import pandas as pd
    tab = pd.crosstab(index=data['Generation'], columns=data[parameter])
    #t_tab = np.transpose(tab[gen-1:gen])
    return tab


def getDatabaseFrequencies(parameter, data):
    import pandas as pd; import numpy as np
    dbProportions=pd.DataFrame()
    def makeCountsDF(parameter, databaseList):
        cols=["Generation","Type","Name","Count","WithinGenPercent"]
        allMeasuresDF=pd.DataFrame(columns=cols)
        subDF = data.loc[:, ['Generation', parameter]].copy()
        for i,gen in enumerate(subDF.Generation.unique().tolist()):
            # DB Counts
            dbCount = pd.Series({w: subDF[subDF.Generation==gen][parameter].str.contains(w, case=False).sum() for w in databaseList})
            # Percent: Within a given generation
            genPercent = dbCount/ dbCount.sum()*100 # Distribution of db percents WITHIN one generation
            Type = parameter.split("_")[0]
            df = pd.DataFrame(np.column_stack([[gen]*len(dbCount),[Type]*len(dbCount), genPercent.index, dbCount, genPercent]),columns=cols)
            allMeasuresDF['WithinGenPercent'] = pd.to_numeric(allMeasuresDF['WithinGenPercent'])
            # Percent: Across all previous generation and the current one
            if i==0: # On the first iteration, allMeasuresDF doesn't have anything in it yet
                df['cumulPercent'] = df['WithinGenPercent']
            else:
                df['cumulPercent'] = allMeasuresDF.groupby('Name')["WithinGenPercent"].mean().values
            # Append to overall DF
            allMeasuresDF = allMeasuresDF.append(df)
        allMeasuresDF["cumulPercent"] = pd.to_numeric(allMeasuresDF['cumulPercent'])
        allMeasuresDF["WithinGenPercent"] = pd.to_numeric(allMeasuresDF['WithinGenPercent'])
        return allMeasuresDF
    if parameter=='PPI_databases':
        PPI_databases = ["BIND", "BIOGRID", "BIOCARTA", "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEGG", \
                         "MINT", "MIPS", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL", "HUMAP", "IREF", "BIOPLEX"]
        dbProportions = dbProportions.append(makeCountsDF(parameter, PPI_databases))
    if parameter=='TF_databases':
        TF_databases = ['chea','transfac','both']
        dbProportions = dbProportions.append(makeCountsDF(parameter, TF_databases))
    if parameter=='KINASE_databases':
        KINASE_databases = ['KP', 'P']
        dbProportions = dbProportions.append(makeCountsDF(parameter, KINASE_databases))
    return dbProportions




# Make super awesome plot showing evolution of parameter distribution over generations/time
def parameterEvolutionPlot(GA_df, figsize=(24,20), padRight=.8, padLeft=.2):
    data = parameterDF(GA_df)
    import matplotlib.pyplot as plt
    import numpy as np
    # Setup dimension vars
    parameters = ['TF_sort','TF_species','TF_databases','TF_topTFs','PPI_databases','PPI_pathLength','PPI_minLength','PPI_maxLength','PPI_finalSize','KINASE_sort','KINASE_databases']
    param_num = len(parameters)
    fitness_rows = 5; PPI_size_rows = 1
    nrows = fitness_rows+PPI_size_rows+param_num+2
    # Setup Fitness data
    Fitness_avg = data.groupby('Generation')["Fitness"].mean()
    Fitness_peak = data.groupby('Generation')["Fitness"].max()

    # Fitness plot
    yerr1 = data.groupby('Generation')["Fitness"].std()
    x = Fitness_avg.index.tolist()
    ax1 = plt.subplot2grid((nrows, 1), (0, 0), facecolor='whitesmoke', rowspan=fitness_rows)
    plt.gcf().set_facecolor('white')  # Change plot border background color
    plt.errorbar(x, Fitness_avg, yerr=yerr1, color='blue', marker='o',markersize=5, capsize=2, label=" Average Fitness")
    ax1.plot(x, Fitness_peak, 'purple', linestyle='-', marker='^', markersize=5, label="Peak Fitness")
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
    PPI_std = data.groupby('Generation')["PPI_size"].std()
    PPI_avg = data.groupby('Generation')["PPI_size"].mean()
    # PPI Size plot
    ax2 = plt.subplot2grid((nrows, 1), ((fitness_rows+PPI_size_rows), 0), facecolor='whitesmoke', rowspan=PPI_size_rows+1)
    ax2.plot(PPI_avg.index, PPI_avg, 'firebrick', marker='d')
    ax2.errorbar(x=PPI_avg.index, y=PPI_avg, yerr=PPI_std, capsize=3, color='mediumpurple')
    plt.ylabel('PPI Size', rotation=0, labelpad=50, fontsize=12)
    plt.tick_params(axis='x', labelbottom='off')
    plt.tick_params(axis='y', labelsize=7)
    #plt.yticks(np.arange(0, max(data['Average_PPI_size']), 3500))
    plt.xticks( np.arange(0, max(x) + 1, 5) )

    # Parameter plots
    import seaborn as sns
    for i,parameter in enumerate(parameters):
        row_num = i+(fitness_rows+PPI_size_rows)+2 # Skip the first 3 rows since the other plot is occupying that space
        new_ax = 'ax' + str(row_num)
        if parameter in ["PPI_databases"]:
            dbProportions = getDatabaseFrequencies(parameter, data).loc[:,['Generation','Name','WithinGenPercent']] # cumulPercent
            dbNames = dbProportions.Name.unique().tolist()
            unstackedData = dbProportions.groupby(['Generation','Name']).mean().unstack()
            sns.set_palette(sns.color_palette("hls", len(dbNames) ))  # husl (for equal color intensity as humans percieve them)
            exec(new_ax + " = plt.subplot2grid((nrows, 1), (row_num, 0))")
            unstackedData.plot(kind='area', ax=eval(new_ax), stacked=True)
            leg = plt.legend(loc='center left', bbox_to_anchor=(1, .5), ncol=7, labels=dbNames, fontsize='small', markerscale=1, columnspacing=None)
            # set the linewidth of each legend object
            for legobj in leg.legendHandles:
                legobj.set_linewidth(5.0)
            #plt.legend(["Optimized database combination (of "+str(len(tab1.columns))+"):\n"+bestPPIdbs],loc='center left', bbox_to_anchor=(1, .5), ncol=70)
            plt.xlabel('')
            plt.ylabel(parameter, rotation=0, labelpad=50, fontsize=12)

        else:
            tab1 = freqTable(data, parameter)
            sns.set_palette(sns.color_palette("BuPu", len(tab1.columns)))  # Run this BEFORE you create your subplot
            exec(new_ax + " = plt.subplot2grid((nrows, 1), (row_num, 0))")
            # Run ax variable from string
            tab1.plot(kind='area', ax=eval(new_ax), stacked=True, figsize=figsize)  # , colors=['hotpink','slateblue','navy']
            plt.legend(loc='center left', bbox_to_anchor=(1, .5), ncol=70)  # bbox_to_anchor=(horizontal, vertical)
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


def ParameterViolinPlots(GA_df, numRows=2, numCols=5, figSize=(10,6)):
    import seaborn as sns
    import matplotlib.pyplot as plt
    df = parameterDF(GA_df)
    tested_parameters = ["TF_sort","TF_species","TF_databases","TF_topTFs","PPI_minLength","PPI_maxLength","PPI_finalSize","KINASE_sort","KINASE_databases"]
    data =  df[[col for col in df.columns if col in ["Fitness"]+tested_parameters]] # Filter database by revelvant paramters
    plt.figure(figsize=figSize)
    for i,parameter in enumerate(tested_parameters):
        # Turn ax string into variable
        #sns.set_palette(sns.color_palette("BuPu", len(data[parameter].unique()) ))  # Run this BEFORE you create your subplot
        ax = plt.subplot(numRows, numCols, i+1)
        #pd.DataFrame.boxplot(data, column='Fitness', by=parameter, ax=ax, grid=False, notch=True, patch_artist=True, rot=45)
        sns.violinplot(x=parameter, y="Fitness", data=data, palette="BuPu", ax=ax, boxprops=dict(alpha=.7))
        plt.xticks(rotation=45, ha='right')
        plt.title(parameter)
        plt.ylabel(''); plt.xlabel('')
        if i == 0 or i == 6:
            plt.ylabel('Fitness')
    plt.subplots_adjust(hspace=.75, wspace=.5, bottom=.2)
    plt.gcf().set_facecolor('white')




def fitnessHistogramCurves(GA_df, genSpacing=2, figSize=(10,6)):
    df = GA_df.copy()
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    import seaborn as sns
    import numpy as np
    import pandas as pd
    df["Generation"] = pd.to_numeric(df.Generation)
    if genSpacing==0:
        whichGens=df.Generation.unique()
    else:
        count = 1
        whichGens = [1]
        while count <= len(df.Generation.unique())-genSpacing-1:
            whichGens.append(count+genSpacing)
            count += genSpacing
    color = cm.rainbow(np.linspace(0, 1, len(whichGens)))

    plt.figure()
    for i,gen in enumerate(whichGens):
        datos=df[df.Generation==gen].Fitness
        # Use kernal density plot function from seaborn
        sns.kdeplot(datos, shade=True, color=color[i], bw=1, label="Gen: "+str(gen), gridsize=100)
    plt.xlabel('Fitness')
    plt.ylabel('Frequency')
    #plt.title('Distribution of All Fitnesses: Every %.0f Generation(s)' % (gen_spacing))
    plt.gcf().set_facecolor('white')
    return plt
#fitnessHistogramCurves(allFitnesses, gen_spacing=1)





def parameterStats(GA_df, writeExcel='No'):
    import pandas as pd
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    import numpy as np
    def anovaTable(Data, independentVar, dependentVar, dependentVarName, dependentVarLabel):
        # One way ANOVA
        mod = ols(independentVar+' ~ '+dependentVar,data=Data).fit()
        aov_table = sm.stats.anova_lm(mod, typ=2)
        #esq_sm = aov_table['sum_sq'][0] / (aov_table['sum_sq'][0] + aov_table['sum_sq'][1])
        p = aov_table['PR(>F)'][0]
        # Add P-val summary
        if p >= 0.05:
            aov_table['Sig'] = "non-sig"
            P = "≥ 0.05"
        if p < 0.05:
            aov_table['Sig'] = '*'
            P = "< 0.05"
        if p < 0.01:
            aov_table['Sig'] = '**'
            P = "< 0.01"
        if p < 0.001:
            aov_table['Sig'] ='***'
            P = "< 0.001"
        if p < 0.0001:
            aov_table['Sig'] = '****'
            P = "< 0.0001"
        # mod.summary() # For full summary
        SS = int(round(float(aov_table.iloc[0].sum_sq), 3))
        Sig= aov_table['Sig'][0]
        F = int(round(float(aov_table.iloc[0].F),0))
        df = int(aov_table.df[0])
        newTable = pd.DataFrame(np.column_stack([dependentVarName, SS, df, F, P, Sig]), columns=[dependentVarLabel,"SS","DF","F-value","P-value","Sig."])
        return newTable
    # Fitness vs. Parameter selection
    data = parameterDF(GA_df)
    #tested_parameters = [x for x in data.columns if x not in ['KINASE_topKinases','KINASE_background','TF_background','PPI_pathLength', 'newBinary','Generation', 'Fitness']]
    tested_parameters=['PPI_size','TF_sort','TF_databases','TF_topTFs','PPI_databases','PPI_minLength','PPI_maxLength','PPI_finalSize',\
                       'KINASE_sort','KINASE_databases'] #TF_background, PPI_pathLength, KINASE_topKinases, KINASE_background
    parameterResults=pd.DataFrame()
    for parameter in tested_parameters:
        parameterResults = parameterResults.append(anovaTable(data, "Fitness", parameter, parameter, "X2K Parameter"))
    # For each generation test Actual Fitness vs. Baselines Fitness
    baselineResults = pd.DataFrame()
    for gen in pd.to_numeric(GA_df.Generation).unique():
        genSub = data[data.Generation==gen]
        baselineResults = baselineResults.append(anovaTable(genSub, "Fitness", "baselineFitness", gen, "Generation"))

    if writeExcel != 'No':
        # Write AOV results to excel file
        print("Writing AOV results to excel file...")
        writer = pd.ExcelWriter(writeExcel)
        print("***** Fitness.Vs.Parameters *****")
        parameterResults.to_excel(writer, 'Fitness.Vs.Parameters')
        print()
        print("***** Fitness.Vs.baselineFitness *****")
        baselineResults.to_excel(writer, 'Fitness.Vs.baselineFitness')
        writer.save()
    return parameterResults, baselineResults
# parameter_AOV_results = parameterStats(GAresults, True)





# Get all individuals with the peak fitness in the final generation
def topParametersReport(GA_df):
    fittestDF = GA_df.sort_values(by=['Fitness'], ascending=False).iloc[:10, :]
    data = parameterDF(fittestDF)
    uniqueOptimizations = data.iloc[:,4:].drop_duplicates()

    # if len(uniqueOptimizations)>1:
    #     """
    #     Get the frequency of each database used amongst only the fittest individuals
    #     to determine how much variation there was for the "best" parameter choice
    #     """
    #     ppiDbFreq  = getDatabaseFrequencies("PPI_databases",data).groupby('Name').mean()
    #     freqTable("TF_sort",data)
    return uniqueOptimizations

# Plot different aspects of fitness for Training Data, Test Data, and baselineFitness
def plotAverageFitness(GA_df1, GA_df2, barsOrFill):
    import matplotlib.pyplot as plt
    import numpy as np
    data1 = parameterDF(GA_df1)
    data2 = parameterDF(GA_df2)

    fig, ax1 = plt.subplots(nrows=1, ncols=1, sharey=True)
    def tsplot(ax, data, measure, linestyle='', color='', marker='',barsOrFill="fill", **kw):
        est =  data.groupby('Generation')[measure].mean()
        x = range(1, len(est)+1)
        sd = data.groupby('Generation')[measure].std()
        cis = (est - sd, est + sd)
        if barsOrFill=="fill":
            ax.fill_between(x, cis[0], cis[1], alpha=0.2, color=color, **kw)
        elif barsOrFill=="bars":
            ax.errorbar(x, est, yerr=sd, color=color, marker=marker, capsize=2)
        ax.plot(x, est, linestyle=linestyle, color=color, marker=marker)
    errorSelection="fill"
    tsplot(ax1, data1, measure="Fitness", linestyle="-", color='c', marker='o', barsOrFill=errorSelection)
    tsplot(ax1, data2, measure="Fitness", linestyle="-", color='g', marker='s', barsOrFill=errorSelection)
    tsplot(ax1, data1, measure="baselineFitness", linestyle="--", color='r', marker='x', barsOrFill=errorSelection)
    ax1.legend(loc='lower center', borderaxespad=0.5, labels=["Training","Test Data","Randomized Baseline"], ncol=70)
    plt.xlabel('Generation')
    plt.ylabel('Average Fitness')
    plt.xticks(np.arange(0, max(data1.Generation) + 1, 5))
    plt.ylim([0,30])
    plt.gcf().set_facecolor('white')  # Change plot border background color


def plotPeakFitness(GA_df1, GA_df2):
    import matplotlib.pyplot as plt
    import numpy as np
    data1 = parameterDF(GA_df1)
    data2 = parameterDF(GA_df2)
    # Peak Fitness plot
    fig, ax2 = plt.subplots(nrows=1, ncols=1, sharey=True)
    ## Subset 1
    y1 = data1.groupby('Generation')["Fitness"].max()
    x = range(1, len(y1)+1)
    ## Subset 2
    y2 = data2.groupby('Generation')["Fitness"].max()
    ## Randomized Baseline
    y3 = data1.groupby('Generation')["baselineFitness"].max()
    # Plot
    ax2.errorbar(x, y1, color='c', marker='o', capsize=2, label="Training Data")
    ax2.errorbar(x, y2, color='g', marker='s', capsize=2, label="Test Data")
    ax2.errorbar(x, y3, linestyle="--", color='r', marker='x',capsize=2, label="Randomized Baseline")
    ax2.legend(loc='lower right', borderaxespad=2)
    plt.xlabel('Generation')
    plt.ylabel('Peak Fitness')
    plt.ylim([0,30])
    plt.xticks(np.arange(0, max(x) + 1, 5))
    plt.gcf().set_facecolor('white')  # Change plot border background color

