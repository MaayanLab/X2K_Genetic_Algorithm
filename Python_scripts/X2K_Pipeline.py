#############################################################################################
##################### X2K: Pipeline & Fitness Function #####################
#############################################################################################


# For checking whether the bit_dict is recording the right parameters
# def makeString(length):
#     import string
#     import random
#     letters = []
#     for i in range(length):
#         letters.append(random.choice(string.ascii_lowercase))
#     return ''.join(letters)
# letterString = makeString(43)


directory = "/Users/schilder/Desktop/X2K_Genetic_Algorithm/data/" # Change to your directory
#directory='/Users/maayanlab/PycharmProjects/X2K_Genetic_Algorithm/data/'
#directory='/Users/maayanlab/Desktop/X2K_Genetic_Algorithm/data/'


# binaryString = '0001011000100101100010011100011101100101011'

# X2K Fitness Function for GA
def X2K_fitness(binaryString, fitnessMethod):
    ######################################
    #### Convert Binary to Parameters ####
    ######################################
    constantBackground = "humanarchs4"

    # CHEA OPTIONS ############
    if constantBackground == False:
        TF_background = {"10": "humanarchs4", "01": "mousearchs4", "11": "hg133", "00": "RESHUFFLE"}
    else:
        TF_background = {"10": constantBackground, "01": constantBackground, "11": constantBackground,"00": constantBackground}
    TF_sort = {"10": "oddsratio", "01": "combined_score", "11": "rank", "00": "pvalue"}
    TF_species = {"10": "human", "01": "mouse", "11": "both", "00": "RESHUFFLE"}
    TF_topTFs = {"10": 5, "01": 10, "11": 20, "00": 30}
    #TF_databases = ["ARCHS4", "CHEA", "CREEDS", "ENCODE", "ENRICHR", "HUMAP", "IREF", "JASPAR-TRANSFAC", "TF-PPI","TF-LOF"]
    TF_databases = {"10": "chea", "01": "transfac", "11": "both", "00": "RESHUFFLE"}

    # G2N OPTIONS ############
    PPI_pathLength = {"0":1, "1":1}
    PPI_minLength = {"10":1, "01":5, "11":10, "00":20}
    PPI_maxLength = {"10":50, "01":100, "11":200, "00":500}
    PPI_finalSize = {"10":50, "01":100, "11":200, "00":500}
    PPI_databases = ["BIND", "BIOGRID", "BIOCARTA",  "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEGG", \
                     "MINT", "MIPS", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL", "HUMAP", "IREF", "BIOPLEX"]  # "KEA", "MURPHY", **** ,
    # KEA OPTIONS ############
    if constantBackground==False:
        KINASE_background = {"10": "humanarchs4", "01": "mousearchs4", "11": "hg133", "00": "RESHUFFLE"}
    else:
        KINASE_background = {"10": constantBackground, "01": constantBackground, "11": constantBackground, "00": constantBackground}
    KINASE_topKinases = 20  # Should be consistent for every individual
    KINASE_sort = {"10": "oddsratio", "01": "combined_score", "11": "rank", "00": "pvalue"}
    #KINASE_databases = ["CREEDS", "iPTMnet", "iREF", "MINT", "HPRD", "PHOSPHOELM", "PHOSPHOPOINT", "PHOSPHOSITE",  "PHOSPHOSITEPLUS", "NETWORKIN"]
    KINASE_databases = {"10": "KP", "01": "P", "11": "RESHUFFLE", "00": "RESHUFFLE"}

    # Flexibly create a dictionary for each bit value and its meaning depending on the particular set of variables
    parameterList = ["TF_sort", "TF_species","TF_databases","TF_background","TF_topTFs", "PPI_databases","PPI_pathLength","PPI_minLength","PPI_maxLength","PPI_finalSize","KINASE_sort","KINASE_databases","KINASE_background"]
    #def constructBitDict(parameterList, binaryString):
    # import pandas as pd
    # import numpy as np
    # bit_dict={}; keyDF=pd.DataFrame()
    # prevBits=0
    # for param in parameterList:
    #     realParam = eval(param)
    #     # Get the number of necessary bits for the specific parameter
    #     if type(realParam)==dict:
    #         numBits = len(next(iter(realParam.keys())))
    #     elif type(realParam)==list:
    #         numBits = len(realParam)
    #     # Assign positions of first and last bits within binaryString
    #     firstBit = prevBits
    #     lastBit = firstBit + numBits
    #     bit_dict[param] = binaryString[firstBit:lastBit]
    #     keyDF = keyDF.append(pd.DataFrame(np.column_stack([param, binaryString[firstBit:lastBit], str(firstBit)+":"+str(lastBit)]), columns=["Parameter", "binary", "indices"]) )
    #     prevBits = prevBits+numBits
        #return bitDict
    #bit_dict = constructBitDict(parameterList, binaryString)
    import pandas as pd
    import numpy as np
    keyDF=pd.DataFrame()
    prevBits=0
    for param in parameterList:
        realParam = eval(param)
        # Get the number of necessary bits for the specific parameter
        if type(realParam)==dict:
            numBits = len(next(iter(realParam.keys())))
        elif type(realParam)==list:
            numBits = len(realParam)
        # Assign positions of first and last bits within binaryString
        firstBit = prevBits
        lastBit = firstBit + numBits
        keyDF = keyDF.append(pd.DataFrame(np.column_stack([param, binaryString[firstBit:lastBit], str(firstBit)+":"+str(lastBit)]), columns=["Parameter", "binary", "indices"]) )
        prevBits = prevBits+numBits



    def totalBitLength(keyDF):
        totalLength=0
        for string in keyDF.binary:
            totalLength = totalLength+len(string)
        print(totalLength)
    # totalBitLength(keyDF)


    # General function for converting binary into list of selected databases
    # def selectDatabases(databaseType, fullDatabaseList):
    #     dbs = []; newDict ={}
    #     dbString = bit_dict[databaseType+"_databases"]
    #     from random import randint
    #     while sum(map(int, dbString)) == 0:
    #         binaryList = list(dbString)
    #         binaryList[randint(0, (len(binaryList)) - 1)] = '1'
    #         dbString = "".join(binaryList)
    #         bit_dict[databaseType+"_databases"] = dbString
    #     # Add to new dictionary
    #     newDict[databaseType+"_databases"] = dbString
    #     ### Generate list of G2N databases
    #     for ind, bit in enumerate(bit_dict[databaseType+"_databases"]):
    #         if bit == "1":
    #             dbs.append(fullDatabaseList[ind])
    #     selectedDatabases = ",".join(dbs)
    #     return selectedDatabases, newDict

    ## Create database list and update keyDF (in case DB list was all zeros)
    def selectDatabases(databaseType, fullDatabaseList, keyDF):
        dbs = []
        dbString  = keyDF[keyDF.Parameter==databaseType+"_databases"].binary[0]
        from random import randint
        while sum(map(int, dbString)) == 0:
            binaryList = list(dbString)
            binaryList[randint(0, (len(binaryList)) - 1)] = '1'
            dbString = "".join(binaryList)
            keyDF.loc[keyDF.Parameter==databaseType+"_databases",'binary'] = dbString
        ### Generate list of G2N databases
        for ind, bit in enumerate(dbString):
            if bit == "1":
                dbs.append(fullDatabaseList[ind])
        selectedDatabases = ",".join(dbs)
        return selectedDatabases, keyDF


    # ############ CHEA Databases ############
    # selectedTFdatabases = selectDatabases("TF", TF_databases)
    # ############ G2N Databases ############
    selectedPPIdatabases, keyDF = selectDatabases("PPI",PPI_databases, keyDF)
    # ############ KEA Databases ############
    # selectedKINASEdatabases = selectDatabases("KINASE", KINASE_databases)

    #RESHUFFLE any bits that aren't legitimate options
    #from random import choice
    # for param in parameterList:
    #     if type(eval(param))==dict:
    #         bad_bits = [k for k, v in eval(param).items() if v == "RESHUFFLE"]
    #         good_bits = [k for k, v in eval(param).items() if v != "RESHUFFLE"]
    #         current_bits = bit_dict[param]
    #         if current_bits in bad_bits:
    #             bit_dict[param] = choice(good_bits)
    #             newDict[param] = choice(good_bits)
    #         else:
    #             newDict[param] = bit_dict[param]
    #RESHUFFLE any bits that aren't legitimate options
    from random import choice
    for param in parameterList:
        if type(eval(param)) == dict:
            bad_bits = [k for k, v in eval(param).items() if v == "RESHUFFLE"]
            good_bits = [k for k, v in eval(param).items() if v != "RESHUFFLE"]
            current_bits = keyDF.loc[keyDF.Parameter==param,'binary'].copy()[0]
            if current_bits in bad_bits:
                keyDF.loc[keyDF.Parameter==param,'binary'] = choice(good_bits)

    # Reconstruct corrected binary string
    newBinary=''
    for i,row in keyDF.iterrows():
        param = row.Parameter
        newBinary = keyDF.loc[keyDF.Parameter==param,'binary'].copy()[0] + newBinary


    ############################################################
    ##################### RUN X2K PIPELINE #####################
    ############################################################
    def bits(parameter):
        bits = keyDF.loc[keyDF.Parameter==parameter,'binary'][0]
        return bits
    # def getParam(parameter):
    #     bits = keyDF.loc[keyDF.Parameter==parameter,'binary'][0]
    #     return str(eval(parameter,globals())[bits])

    import time
    import socket
    HOST = "localhost"
    PORT = 5000
    PORT2 = 5001
    PORT3 = 5002
    start_time = time.time()

    print("Initiating X2K . . . ")
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect((HOST, PORT))
    buffer_size = 1024
    data = ""
    allData = ""

    print("Running ChEA...")
    chea_parameters = ';'.join( ["run", TF_sort[bits("TF_sort")], TF_species[bits("TF_species")], TF_databases[bits("TF_databases")], TF_background[bits("TF_background")], str(TF_topTFs[bits("TF_topTFs")]) ] ) #"run,oddsratio,both,chea,humanarchs4,20")
    print("CHEA Parameters:   "+chea_parameters)
    sock.sendall(bytes(chea_parameters+"\n", 'utf-8'))

    while 1:
        data = sock.recv(buffer_size).decode("utf-8")
        allData = allData + data
        if allData.endswith("messageComplete\n"):
            break

    sock.send(bytes("kill\n", 'utf-8'))
    sock.close()

    allData = allData.replace(chea_parameters+"\n", "")
    text_file = open(directory+"output/chea_out.txt", "w")
    text_file.write(allData)
    text_file.close()

    sock2 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock2.connect((HOST, PORT2))

    buffer_size = 1024
    data = ""
    allData2 = ""

    allData.replace("messageComplete\n", "")

    print("Running G2N...")
    g2n_string = ';'.join(["run", selectedPPIdatabases, str(PPI_pathLength[bits("PPI_pathLength")]), str(PPI_minLength[bits("PPI_minLength")]), str(PPI_maxLength[bits("PPI_maxLength")]), str(PPI_finalSize[bits("PPI_finalSize")]) ])
    print("G2N Parameters:   "+g2n_string)
    g2n_parameters = g2n_string + "\n"+allData+"messageComplete\n"
    sock2.sendall(bytes(g2n_parameters+"\n", 'utf-8'))

    while 1:
        #print("d: "+data)
        data = sock2.recv(buffer_size).decode("utf-8")
        allData2 = allData2 + data
        if allData2.endswith("messageComplete\n"):
            break

    allData2.replace("messageComplete\n", "")

    sock2.send(bytes("kill\n", 'utf-8'))
    sock2.close()

    # ******************************* #
    # *** Record average PPI size ***
    PPIs = allData2.split('\n')[:-2]# Get rid of messageComplete and '' lines at end
    ppi_sizes = []
    for i,ppi in enumerate(PPIs):
        split_outcome = ppi.split(',')
        if len(split_outcome) == 1: # Will get just the row name if there's no genes in the G2N output
            ppi_sizes.append(0)
        else:
            ppi_sizes.append( len(ppi.split(','))-1 )
    average_PPI_size = sum(ppi_sizes) / float(len(ppi_sizes))
   # ***                          ***
   # ******************************* #


    allData2 = allData2.replace("messageComplete\n", "").replace(g2n_parameters, "")
    text_file = open(directory+"output/g2n_out.txt", "w")
    text_file.write(allData2)
    text_file.close()

    sock2 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock2.connect((HOST, PORT3))

    buffer_size = 1024
    data = ""
    allData3 = ""

    print("Running KEA...")
    kea_string = ';'.join(["run", KINASE_sort[bits("KINASE_sort")], KINASE_background[bits("KINASE_background")], KINASE_databases[bits("KINASE_databases")], str(KINASE_topKinases)  ])
    print("KEA Parameters:   "+kea_string)
    kea_parameters = kea_string + "\n"+allData2+"messageComplete\n"
    kea_parameters.replace("messageComplete\n", "")
    sock2.sendall(bytes(kea_parameters+"\n", 'utf-8'))

    while 1:
        #print("d: "+data)
        data = sock2.recv(buffer_size).decode("utf-8")
        allData3 = allData3 + data
        if allData3.endswith("messageComplete\n"):
            break

    allData3.replace("messageComplete\n", "")

    sock2.send(bytes("kill\n", 'utf-8'))
    sock2.close()

    allData3 = allData3.replace("messageComplete\n", "").replace(kea_parameters, "")
    text_file = open(directory+"output/kea_out.txt", "w")
    text_file.write(allData3)
    text_file.close()

    time.time() - start_time





    ##################### FITNESS CALCULATION/VALIDATION #####################

    def detectValidationData():
        import os
        gmtName = os.listdir("data/testgmt")
        if gmtName[0].startswith("GEO-KinasePert_L1000-DRH"):
            dataType = 'GEO-L1000'
        if gmtName[0].startswith('Kinase_Perturbations_from_GEO'):
            dataType = 'GEO'
        if gmtName[0].startswith('Chem_combo_DRH'):
            dataType = 'L1000'
        else:
            dataType='GEO'
            print("Couldn't detect Validation Data type. Defaulting to GEO")
        return(dataType)
    dataType = detectValidationData()

    def parseTargetsPredicted(line, dataType):
        if dataType == 'GEO':
            kinaseTargets = line.split(",")[0].split("_")[0]
            predictedKinases = line.split(",")[2:]
        elif dataType == 'L1000':
            kinaseTargets = line.split(",")[0].split("_")[3].strip("[").strip("]").split("|")
            predictedKinases = line.split(",")[2:]
        elif dataType == 'GEO-L1000':
            kinaseTargets = line.split(",")[0].split("_")[-1].strip("[").strip("]").split("|")
            predictedKinases = line.split(",")[2:]
        else:
            print("Unrecognized Validation Data Structure")
        # If it's just a string, convert to list
        if type(kinaseTargets)==str:
            kinaseTargets = [kinaseTargets]
        if type(predictedKinases)==str:
            predictedKinases = [predictedKinases]
        # Strip away any lingering \n
        if kinaseTargets != []:
            kinaseTargets[-1] = kinaseTargets[-1].strip()
        if predictedKinases != []:
            predictedKinases[-1] = predictedKinases[-1].strip()
        return( [kinaseTargets, predictedKinases] )

    # Get synonyms
    import pickle
    with open('../X2K_Summaries/General_Resources/synDict.pkl', 'rb') as f:
        synDict = pickle.load(f)

    def getSyn(gene):
        if gene in synDict.keys():
            syns = synDict[gene]
        else:
            syns =[gene]
        return syns

    def randomizedBaselineFitness(KEAlines, fitnessMETHOD, dataType):
        from random import choice
        shuffled_KEAlines = []
        """
        Names = []; KinasePredictions = []; 
        for line in KEAlines:
            lineSp = line.split(',')
            Names.append(lineSp[0])
            KinasePredictions.append(lineSp[1:])
        shuffle(Names)
        """
        # Get list of ALL kinases
        import pandas as pd
        KINASES = pd.read_table("../X2K_Summaries/General_Resources/Kinase_Families_Maxime/kinases_fam.tsv", \
                                header=None, names=["Name", "Group", "Family"], index_col=False)['Name'].tolist()
        for line in KEAlines:
            info = line.split(',')[0].split("_")[1:]
            newKinase = choice(KINASES) # Randomly select any kinase
            KinasePredictions = line.split(",")[1:]
            shuffled_KEAlines.append(newKinase+"_"+"_".join(info)+ "," + ",".join(KinasePredictions))
        # Calculate randomized baseline fitness
        baselineFitness, xTKs, xPKs = fitnessMETHOD(shuffled_KEAlines, dataType)
        return baselineFitness


    #@@@@@ FITNESS METHODS #@@@@@

    # [1] Simple Presence/absence fitness:
    ## Percentage of experiments whose target kinase(s) were in the predicted kinases.
    # Makes the most sense when there's only one target, so it's consistent across experiments.
    def targetAdjustedOverlap(KEAlines, dataType):
        scaledOverlapScores=[]; TKs=[]; PKs=[]
        for line in KEAlines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            TKs.append(targetKinases); PKs.append(predictedKinases)
            # Create lists of recovered and missed kinases
            overlappingKinases = list( set(targetKinases) & set(predictedKinases) )
            scaledOverlapScores.append( len(overlappingKinases) / len(targetKinases) *100 )
        # Calculate individual's  fitness score
        if scaledOverlapScores!=[]:
            fitness_score = sum(scaledOverlapScores) / len(scaledOverlapScores)
        else:
            fitness_score = 0.0
        return fitness_score, TKs, PKs

    def targetAdjustedOverlap_outputLengthCorrection(KEAlines, dataType):
        scaledOverlapScores = [];
        TKs = [];
        PKs = []
        for line in KEAlines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            TKs.append(targetKinases); PKs.append(predictedKinases)
            # Create lists of recovered and missed kinases
            overlappingKinases = list(set(targetKinases) & set(predictedKinases))
            if len(overlappingKinases)>0:
                scaledOverlapScores.append(len(overlappingKinases) / len(targetKinases) / len(predictedKinases) * 100)
            else:
                scaledOverlapScores.append(0.0)
        # Calculate individual's  fitness score
        if scaledOverlapScores != []:
            fitness_score = sum(scaledOverlapScores) / len(scaledOverlapScores)
        else:
            fitness_score = 0.0
        return fitness_score, TKs, PKs

    def WilcoxonRankSum(KEAlines, dataType):
        # Assumes BOTH lists are ordered (would work for KINOMEscan)
        import scipy.stats as ss
        import numpy as np
        stats=[]; TKs=[]; PKs=[]
        for line in KEAlines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            TKs.append(targetKinases); PKs.append(predictedKinases)
            wrm = ss.ranksums(targetKinases, predictedKinases)
            if np.isnan(wrm.statistic):
                stats.append(0)
            else:
                stats.append(wrm.statistic)
        fitness_score = np.nanmean(stats)
        return fitness_score, TKs, PKs

    def rankBiasedOverlap(KEALines, dataType):
        from Python_scripts.rbo import rbo
        rboScores = []; TKs=[]; PKs=[]
        for line in KEALines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            TKs.append(targetKinases); PKs.append(predictedKinases)
            # As long as the kinase list is not blank, get rid of \n in last item
            if predictedKinases != []:
                predictedKinases[-1] = predictedKinases[-1].strip()
            # Conduct RankedBasedOrder statistic (0-1 score)
            if len(predictedKinases) > 0:
                rboResults = rbo(targetKinases, predictedKinases, .9)
                rboScores.append(rboResults['res'])
            else:
                rboScores.append(0)
            # Get the average RBO score
        fitness_score = sum(rboScores) / len(rboScores)
        return fitness_score, TKs, PKs

    # [4] Modified RankBiasedOverlap
    ## RBO normally assumes both lists are ranked. But in some cases, they're not.
    # So to account for this I'm repeating the RBO for each kinaseTarget individually and then taking the average RBO score.
    ## Developed by Brian M. Schilder & Moshe Silverstein
    def modifiedRBO(KEAlines, dataType):
        from Python_scripts.rboScore import getRBOScore
        rboScores=[]; TKs=[]; PKs=[]
        for line in KEAlines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            TKs.append(targetKinases); PKs.append(predictedKinases)
            # If there's any target hits in the predicted output, calculate score. Otherwise, save time by just appending a 0.0.
            if len(set(targetKinases).intersection(set(predictedKinases))) > 0:
                rboScores.append(getRBOScore(predicted=predictedKinases, target=targetKinases))
            else:
                rboScores.append(0.0)
        fitness_score = sum(rboScores) / len(rboScores)
        return fitness_score, TKs, PKs

    def rankCorrectedTAO(KEAlines, dataType):
        import pandas as pd
        correctionDF = pd.read_csv("Validation/predictedKinaseCorrectionFactors.csv")
        correctedScores = []; TKs = []; PKs = []
        for line in KEAlines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            TKs.append(targetKinases); PKs.append(predictedKinases)
            # Create lists of recovered and missed kinases
            overlappingKinases = list(set(targetKinases) & set(predictedKinases))
            if len(overlappingKinases)>0:
                for hit in overlappingKinases:
                    if hit in correctionDF['Gene'].tolist():
                        correctedScores.append(100 * correctionDF[correctionDF['Gene'] == hit]['correctionFactor'].values[0] / len(targetKinases))
                    else:
                        correctedScores.append(100 / len(targetKinases)) # If there's no associated correction factor for a hit, just assume it has full weight (1.0)
            else:
                correctedScores.append(0)
        fitnessScore = sum(correctedScores) / len(correctedScores)
        return fitnessScore, TKs, PKs

    def simpleRankFitness(KEAlines, dataType):
        from random import randint
        ranks = [];TKs = []; PKs = []
        for line in KEAlines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            TKs.append(targetKinases); PKs.append(predictedKinases)
            # Add target kinase synonyms
            targetKinases_syns=[]
            for tk in targetKinases:
                targetKinases_syns = list(set(targetKinases_syns + getSyn(tk)))
            # Get the ranks of each hit
            kinaseHits  = set(targetKinases_syns).intersection(set(predictedKinases))
            for i in range(len(targetKinases)): # iterate over each tk for L1000-DRH
                if len(kinaseHits)>0:
                    for hit in kinaseHits:
                        ranks.append( predictedKinases.index(hit) )
                else:
                    # If not a hit, randomly select a rank from all the non-hits
                    ranks.append(randint(0,473-len(predictedKinases)) )
        fitnessScore = sum(ranks)/len(ranks)
        converrtedFitness = (472 - fitnessScore)/472
        return converrtedFitness, TKs, PKs



    # ----------------------------- Calculate fitness -----------------------------#
    def calculateFitness(fitnessMETHOD, dataType):
        print("Calculating '" + str(fitnessMethod) + "' fitness...")
        with open(directory + "output/kea_out.txt") as KEA_output:
            KEAlines = KEA_output.readlines()
        fitnessScore, TKs, PKs = fitnessMETHOD(KEAlines, dataType)
        baselineFitness = randomizedBaselineFitness(KEAlines, fitnessMETHOD, dataType)
        print("Individual's fitness = " + str(fitnessScore))
        print("Baseline fitness = " + str(baselineFitness))
        return fitnessScore, baselineFitness, TKs, PKs
    fitnessScore, baselineFitness, TKs, PKs = calculateFitness(fitnessMETHOD=eval(fitnessMethod), dataType=dataType)

    def nestedListToString(nestedList):
        tmpList=[]
        for sublist in nestedList:
            tmpList.append(",".join(sublist))
        return ";".join(tmpList)

    return fitnessScore, average_PPI_size, newBinary, chea_parameters, g2n_string, kea_string, baselineFitness, nestedListToString(TKs), nestedListToString(PKs)








