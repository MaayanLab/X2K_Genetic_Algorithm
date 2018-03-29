

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
    # Get the average number of kinases produced by X2K for each experiment
    numKinases = GA_df['predictedKinases'].str.split(';').apply(lambda parts: [len(x.split(",")) for x in parts]).copy()
    df['avgNumKinases'] = numKinases.apply(lambda x: (sum(x) / len(x)))

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
def parameterEvolutionPlot(GA_df, figsize=(24,30), padRight=.8, padLeft=.2):
    data = parameterDF(GA_df)
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools
    # Setup dimension vars
    parameters = ['TF_sort','TF_species','TF_databases','TF_topTFs','PPI_databases','PPI_pathLength','PPI_minLength','PPI_maxLength','PPI_finalSize','KINASE_sort','KINASE_databases']
    param_num = len(parameters)
    fitness_rows = 5; PPI_size_rows = 1
    ppi_dbs = len(getDatabaseFrequencies("PPI_databases", data).Name.unique())-1
    nrows = fitness_rows+PPI_size_rows+param_num+2+ppi_dbs
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
            dbProportions = getDatabaseFrequencies(parameter, data).loc[:,['Generation','Name','cumulPercent']] # withinGenPercent
            dbNames = dbProportions.Name.unique().tolist()
            """
            #---- Plot cumulative or withinGen percent (all dbs in same subplot)
            unstackedData = dbProportions.groupby(['Generation','Name']).mean().unstack()
            sns.set_palette(sns.color_palette("hls", len(dbNames) ))  # husl (for equal color intensity as humans percieve them)
            exec(new_ax + " = plt.subplot2grid((nrows, 1), (row_num, 0))")
            unstackedData.plot(kind='area', ax=eval(new_ax), stacked=True)
            leg = plt.legend(loc='center left', bbox_to_anchor=(1, .5), ncol=7, labels=dbNames, fontsize='small', markerscale=1, columnspacing=None)
            # plt.legend(["Optimized database combination (of "+str(len(tab1.columns))+"):\n"+bestPPIdbs],loc='center left', bbox_to_anchor=(1, .5), ncol=70)
            # set the linewidth of each legend object
            for legobj in leg.legendHandles:
                legobj.set_linewidth(5.0)
            plt.xlabel('')
            plt.ylabel(parameter, rotation=0, labelpad=50, fontsize=12)
            """
            #---- Plot cumulative or withinGen percent (each db in different subplot)
            palette = itertools.cycle(sns.color_palette("hls", len(dbNames)))
            for db in dbNames:
                plotCount=1
                row_num+=plotCount
                exec(new_ax + " = plt.subplot2grid((nrows, 1), (row_num, 0))")
                dbSub = dbProportions[dbProportions.Name==db]
                dbSub.plot(x="Generation",y="cumulPercent", kind="area", color=next(palette), label=db, ax=eval(new_ax))
                plt.ylabel(parameter, rotation=0, labelpad=50, fontsize=12)
                plt.legend(loc='center left', bbox_to_anchor=(1, .5))
                plt.tick_params(axis='x', labelbottom='off')
                plt.xlabel('')
                plt.ylabel(parameter, rotation=0, labelpad=50, fontsize=12)
                plotCount+=1
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
            plt.tick_params(axis='x', labelbottom='off')
            plt.xlabel('')
    plt.xticks(np.arange(0, max(x) + 1, 5))
    plt.xlabel("Generation")
    eval(new_ax).yaxis.set_major_locator(plt.NullLocator())
    eval(new_ax).xaxis.set_major_formatter(plt.NullFormatter())
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
def plotFitness(GA_df1, GA_df2, barsOrFill):
    import matplotlib.pyplot as plt
    import numpy as np
    data1 = parameterDF(GA_df1)
    data2 = parameterDF(GA_df2)

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False)
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
    # Average # PKs
    tsplot(ax0, data1, measure="avgNumKinases", linestyle="-", color='m', marker='o', barsOrFill=errorSelection)
    ax0.legend(loc='center left', bbox_to_anchor=(1, .5), labels=['Average #\nPredicted Kinases'])
    ax0.set_ylabel('Average # \nPredicted Kinases')

    # Average fitness
    tsplot(ax1, data1, measure="Fitness", linestyle="-", color='c', marker='o', barsOrFill=errorSelection)
    tsplot(ax1, data2, measure="Fitness", linestyle="-", color='g', marker='s', barsOrFill=errorSelection)
    tsplot(ax1, data1, measure="baselineFitness", linestyle="--", color='r', marker='x', barsOrFill=errorSelection)
    ax1.legend(loc='center left', bbox_to_anchor=(1, .5), labels=["Training Data","Test Data","Baseline Fitness"])
    ax1.set_ylabel('Average Fitness')

    # Peak fitness
    import pandas as pd
    x = pd.to_numeric(data1['Generation'].unique()).tolist()
    ax2.errorbar(x, data1.groupby('Generation')["Fitness"].max(), color='c', marker='o', label="Training Data")
    ax2.errorbar(x, data2.groupby('Generation')["Fitness"].max(), color='g', marker='s', label="Test Data")
    ax2.errorbar(x, data1.groupby('Generation')["baselineFitness"].max(), linestyle="--", color='r', marker='x', label="Randomized Baseline")
    ax2.legend(loc='center left', bbox_to_anchor=(1, .5), labels=["Training Data","Test Data","Baseline Fitness"])
    ax2.set_xlabel('Generation')
    ax2.set_ylabel('Peak Fitness')

    plt.xticks(np.arange(0, max(data1.Generation) + 1, 5))
    plt.gcf().set_facecolor('white')  # Change plot border background color


def estimateParameterCombinations():
    TF_sort = 4; TF_species = 3; Top_TFs = 4
    PPI_pathLength = 2; Min_len = 4; Max_len = 4; Max_PPIsize = 4
    KIN_sort = 4
    PPI_dbs = ["BIND", "BIOGRID", "BIOCARTA",  "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEGG", \
                     "MINT", "MIPS", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL", "HUMAP", "IREF", "BIOPLEX"]
    PPI_combinations = 2**len(PPI_dbs)-1
    TF_combinations = 2**10
    KIN_combinations = 2**10
    totalCombinations = TF_sort * TF_species * Top_TFs * PPI_pathLength * Min_len * Max_len * Max_PPIsize * KIN_sort * PPI_combinations * TF_combinations * KIN_combinations
    len(str(totalCombinations))
    from decimal import Decimal
    print( '%.2E' % Decimal(totalCombinations) )

def checkPredictedKinaseFreq(GA_df):
    from collections import Counter
    import pandas as pd
    # Target Kinases
    TKs = GA_df.iloc[0,]['targetKinases']
    TKdf = pd.DataFrame(Counter(TKs.split(",")).most_common(),columns=["Gene","TKcount"])
    TKdf["TKpercent"] = TKdf['TKcount']/TKdf['TKcount'].sum()*100
    # Predicted Kinases
    allPKs=[]
    for i,row in GA_df.iterrows():
        PKs = row['predictedKinases'].split(",")
        for pk in PKs:
            allPKs.append(pk)
    PKdf = pd.DataFrame( Counter(allPKs).most_common(), columns=["Gene", "PKcount"])
    PKdf["PKpercent"] =  PKdf['PKcount'] / PKdf['PKcount'].sum() * 100
    # Merge the DFs
    mergedDF = pd.merge(TKdf, PKdf, on='Gene').fillna(0)[['Gene','TKpercent','PKpercent']]
    # Plot
    import seaborn as sns
    melted_df = pd.melt(mergedDF, id_vars='Gene', var_name='kinaseCategory', value_name='Percent')
    ## Line plot
    lp = sns.pointplot(data=melted_df, x="Gene", y="Percent", hue="kinaseCategory",
              palette={"TKpercent": "c", "PKpercent": "m"}, markers=["^", "o"], linestyles=["-", "--"], labels=["targetKinase", "predictedKinases"])
    lp.set_xticklabels(g.get_xticklabels(), rotation=45, fontsize=8)
    ## Barplot
    bp = sns.factorplot(x="Gene", y="Percent", hue="kinaseCategory", data=melted_df, kind="bar", palette={"TKpercent": "c", "PKpercent": "m"})
    bp.set_xticklabels(rotation=30, fontsize=8)
    ## Stacked barplot
    melted_df.plot(kind='bar', stacked=True)



def predictedKinaseRanks(GA_df):
    import pandas as pd
    geneList=[]; rankList=[]
    #entryCount=0
    for I,row in GA_df.iterrows():
        eachExperiment = row['predictedKinases'].split(';')
        for experiment in eachExperiment:
            experimentSplit = experiment.split(",")
            for i,gene in enumerate(experimentSplit):
                geneList.append(gene)
                rankList.append(i+1)
                #rankDF.loc[entryCount] = [gene, i+1]
                #entryCount+=1
    import numpy as np
    rankDF = pd.DataFrame(np.column_stack([geneList,rankList]), columns=['Gene','Rank'])
    rankDF['Rank'] = pd.to_numeric(rankDF.Rank)
    meanRankDF = rankDF.groupby('Gene').mean()
    meanRankDF['std'] = rankDF.groupby('Gene')['Rank'].std().fillna(0)
    meanRankDF['count'] = rankDF.groupby('Gene').count()
    ## Calculate correction factor for each gene
    mean_std = meanRankDF['Rank']-meanRankDF['std']
    normalized = (mean_std- min(mean_std)) / (max(mean_std) - min(mean_std))
    meanRankDF['correctionFactor'] = 1-normalized
    meanRankDF = meanRankDF[meanRankDF.index!='']
    meanRankDF.to_csv("Validation/predictedKinaseCorrectionFactors.csv")
    return meanRankDF


