#############################################################################################
##################### X2K: Pipeline & Fitness Function #####################
#############################################################################################
# lsof -i tcp:5001
# COMMAND   PID     USER   FD   TYPE             DEVICE SIZE/OFF NODE NAME
# java    70871 schilder    5u  IPv6 0x7d91bfdfeea593bd      0t0  TCP *:commplex-link (LISTEN)
# kill 70871

# # Create test individual
# from random import choice
# def createPopulation(popSize, parameterLength):
#     populationinit = []
#     for i in range(popSize):
#         populationinit.append(''.join(choice(('0', '1')) for _ in range(parameterLength)) )
#         print(populationinit[i])
#     return populationinit
# parameterLength = 43
# populationInit = createPopulation(10,parameterLength)
# binaryString = populationInit[0]

# For checking whether the bit_dict is recording the right parameters
# def makeString(length):
#     import string
#     import random
#     letters = []
#     for i in range(length):
#         letters.append(random.choice(string.ascii_lowercase))
#     return ''.join(letters)
# binary = makeString(100)


directory = "/Users/schilder/Desktop/X2K_Genetic_Algorithm/data/" # Change to your directory
#directory='/Users/maayanlab/PycharmProjects/X2K_Genetic_Algorithm/data/'

# binaryString = '0110110100001010110101101101100000001011000'

# X2K Fitness Function for GA
def X2K_fitness(binaryString, fitness_method='target-adjusted overlap'):
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
    TF_topTFs = {"10": 5, "01": 10, "11": 20, "00": "RESHUFFLE"}
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
    bit_dict={}
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
        bit_dict[param] = binaryString[firstBit:lastBit]
        prevBits = prevBits+numBits
        #return bitDict
    #bit_dict = constructBitDict(parameterList, binaryString)

    def totalBitLength():
        totalLength=0
        for string in bit_dict.values():
            totalLength = totalLength+len(string)
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
        binaryString = bit_dict[databaseType+"_databases"]
        from random import randint
        while sum(map(int, binaryString)) == 0:
            binaryList = list(binaryString)
            binaryList[randint(0, (len(binaryList)) - 1)] = '1'
            # ******* Turn off specific databases ******* AGAIN just in case it picked one of these
            # ppi_list[2] = '0'  # Turn off BIOGRID PPI
            # ppi_list[8] = '0'  # Turn off KEA PPI
            binaryString = "".join(binaryList)
            bit_dict[databaseType+"_databases"] = binaryString
        ### Generate list of G2N databases
        for ind, bit in enumerate(bit_dict[databaseType+"_databases"]):
            if bit == "1":
                dbs.append(fullDatabaseList[ind])
        selectedDatabases = ",".join(dbs)
        return selectedDatabases

    # ############ CHEA Databases ############
    # selectedTFdatabases = selectDatabases("TF", TF_databases)
    # ############ G2N Databases ############
    selectedPPIdatabases = selectDatabases("PPI",PPI_databases)
    # ############ KEA Databases ############
    # selectedKINASEdatabases = selectDatabases("KINASE", KINASE_databases)

    # Reshuffle
    from random import choice
    # RESHUFFLE any bits that aren't legitimate options
    for param in parameterList:
        if type(eval(param))==dict:
            bad_bits = [k for k, v in eval(param).items() if v == "RESHUFFLE"]
            good_bits = [k for k, v in eval(param).items() if v != "RESHUFFLE"]
            current_bits = bit_dict[param]
            if current_bits in bad_bits:
                bit_dict[param] = choice(good_bits)
            #print("Reshuffling...")
    newBinary = "".join(bit_dict.values())

    ############################################################
    ##################### RUN X2K PIPELINE #####################
    ############################################################
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
    chea_parameters = ';'.join( ["run", TF_sort[bit_dict["TF_sort"]], TF_species[bit_dict["TF_species"]], TF_databases[bit_dict["TF_databases"]], TF_background[bit_dict["TF_background"]], str(TF_topTFs[bit_dict["TF_topTFs"]]) ] ) #"run,oddsratio,both,chea,humanarchs4,20")
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
    g2n_string = ';'.join(["run", selectedPPIdatabases, str(PPI_pathLength[bit_dict["PPI_pathLength"]]), str(PPI_minLength[bit_dict["PPI_minLength"]]), str(PPI_maxLength[bit_dict["PPI_maxLength"]]), str(PPI_finalSize[bit_dict["PPI_finalSize"]]) ])
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
    kea_string = ';'.join(["run", KINASE_sort[bit_dict["KINASE_sort"]], KINASE_background[bit_dict["KINASE_background"]], KINASE_databases[bit_dict["KINASE_databases"]], str(KINASE_topKinases)  ])
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
        elif gmtName[0].startswith('Chem_combo_DRH'):
            dataType = 'L1000'
        else:
            print("Couldn't detect Validation Data type.")
        return(dataType)
    dataType = detectValidationData()

    def parseTargetsPredicted(line, dataType):
        if dataType == 'GEO':
            kinaseTargets = line.split(",")[0].split("_")[0]
            predictedKinases = line.split(",")[1:]
        elif dataType == 'L1000':
            kinaseTargets = line.split(",")[0].split("_")[3].strip("[").strip("]").split("|")
            predictedKinases = line.split(",")[1:]
        elif dataType == 'GEO-L1000':
            kinaseTargets = line.split(",")[0].split("_")[-1].strip("[").strip("]").split("|")
            predictedKinases = line.split(",")[1:]
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

    # [1] Simple Presence/absence fitness:
    # [1] Simple Presence/absence fitness:
    ## Percentage of experiments whose target kinase(s) were in the predicted kinases.
    # Makes the most sense when there's only one target, so it's consisent across experiments.
    if fitness_method == 'target-adjusted overlap':
        print("Calculating 'presence/absence' fitness...")
        with open (directory+"output/kea_out.txt") as KEA_output:
            KEA_lines = KEA_output.readlines()
        scaledOverlapScores=[]

        for line in KEA_lines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            # Create lists of recovered and missed kinases
            overlappingKinases = list( set(targetKinases) & set(predictedKinases) )
            scaledOverlapScores.append( len(overlappingKinases) / len(targetKinases)*100 )
        # Calculate individual's  fitness score
        if scaledOverlapScores!=[]:
            fitness_score = sum(scaledOverlapScores) / len(scaledOverlapScores)
        else:
            fitness_score = 0.0
            print("Empty KEA file...")
        print("Individual's fitness = " + str(fitness_score))


    # [2] Rank-weighted fitness
    ## Weights fitness according to 1) the number of (unranked) targetKinases recovered
    #### and 2) the rank of each Kinase after removing all other targetKinases from the predictedKinases list (to prevent never reaching 1)
    ## Scale from 0 (no targetKinases recovered) to 1 (all targetKinases recovered at the top of the list in any order)
    ## Developed by Brian M. Schilder
    if fitness_method == 'Rank-Weighted Mean':
        print("Calculating <"+fitness_method+"> fitness...")
        with open(directory + "output/kea_out.txt") as KEA_output:
            KEA_lines = KEA_output.readlines()
        for line in KEA_lines:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            # Repeat for each targetKinase
            tkRankRatios=[]
            for tk in targetKinases:
                # Remove any other targetKinases form predictedKinases list
                if len(targetKinases)>1:
                    predictedSubset = list( set(predictedKinases)- (set(targetKinases) - set([tk])) )
                else:
                    predictedSubset = predictedKinases
                # If targetKinase is in predictedKinases find 1/rank
                if tk in predictedSubset:
                    rankRatio = (len(predictedSubset)-predictedSubset.index(tk)) / len(predictedSubset)
                    # Add a power function to make top-weighted (lower ranks count exponentially less)
                    topWeighted_rankRatio = pow(rankRatio, 2)
                    tkRankRatios.append(topWeighted_rankRatio)
                else:
                    tkRankRatios.append(0)
            # if sum(tkRankRatios) !=0:
            #     # Divide by the number of targetKinases to adjust for the number of possible hits
            #     tkRankRatios_avgs.append( sum(tkRankRatios) / len(targetKinases) )
            # else:
            #     tkRankRatios_avgs.append(0)
        fitness_score = sum(tkRankRatios) / len(tkRankRatios)

        print("Individual's fitness = " + str(fitness_score))


    if fitness_method == "Wilcoxon rank-sum":
        print("Calculating <" + fitness_method + "> fitness...")
        with open(directory + "output/kea_out.txt") as KEA_output:
            KEA_lines = KEA_output.readlines()
            stats=[]
            import scipy.stats as ss
            for line in KEA_lines:
                targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
                predictedBits=[]
                targetBits = list('1'*len(targetKinases))
                for gene in predictedKinases:
                    if gene in targetKinases:
                        predictedBits.append(1)
                    else:
                        predictedBits.append(0)
                wrm = ss.ranksums(targetBits, predictedBits)
                stats.append(wrm.statistic)
        fitness_score = sum(stats) / len(stats)
        print("Individual's fitness = " + str(fitness_score))

    # [3] Rank Biased Overlap
    ## From the following blog post/package: https://ragrawal.wordpress.com/2013/01/18/comparing-ranked-list/
    if fitness_method == 'Rank-Biased Overlap':
        print("Calculating "+fitness_method+" fitness...")
        from Python_scripts.rbo import rbo
        with open(directory+"output/kea_out.txt") as KEA_output:
            KEA_out = KEA_output.readlines()
        rboScores = []
        for line in KEA_out:
            targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
            # As long as the kinase list is not blank, get rid of \n in last item
            if predictedKinases != []:
                predictedKinases[-1] = predictedKinases[-1].strip()
            # Conduct RankedBasedOrder statistic (0-1 score)
            if len(predictedKinases)>0:
                rboResults = rbo(targetKinases, predictedKinases, .9)
                rboScores.append(rboResults['res'])
            else: rboScores.append(0)
            # Get the average RBO score
            fitness_score = sum(rboScores) / len(rboScores)

            print("Individual's fitness = " + str(fitness_score))

        # [4] Modified RankBiasedOverlap
        ## RBO normally assumes both lists are ranked. But in some cases, they're not.
        # So to account for this I'm repeating the RBO for each kinaseTarget individually and then taking the average RBO score.
        ## Developed by Brian M. Schilder
        if fitness_method == 'modified RankBiasedOverlap':
            print("Calculating <" + fitness_method + "> fitness...")

            from Python_scripts.rboScore import getRBOScore

            with open(directory + "output/kea_out.txt") as KEA_output:
                KEA_lines = KEA_output.readlines()
            rboScores=[]
            for line in KEA_lines:
                targetKinases, predictedKinases = parseTargetsPredicted(line, dataType)
                print(line)
                # If there's any target hits in the predicted output, calculate score. Otherwise, save time by just appending a 0.0.
                if len(set(targetKinases).intersection(set(predictedKinases))) > 0:
                    rboScores.append( getRBOScore(targetKinases, predictedKinases) )
                else: rboScores.append(0.0)
            # Calculate average rank-weighted fitness
            fitness_score = sum(rboScores) / len(rboScores)
            print("Individual's fitness = " + str(fitness_score))

    return fitness_score, average_PPI_size, newBinary



