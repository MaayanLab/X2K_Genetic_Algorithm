#############################################################################################
##################### X2K: Pipeline & Fitness Function #####################
#############################################################################################

# # Create test individual
# from random import choice
# def createPopulation(popSize, parameterLength):
#     populationinit = []
#     for i in range(popSize):
#         populationinit.append(''.join(choice(('0', '1')) for _ in range(parameterLength)) )
#         print(populationinit[i])
#     return populationinit
# parameterLength = 28
# populationInit = createPopulation(10,parameterLength)
# binary = populationInit[0]

# For checking whether the bit_dict is recording the right parameters
# def makeString(length):
#     import string
#     import random
#     letters = []
#     for i in range(length):
#         letters.append(random.choice(string.ascii_lowercase))
#     return ''.join(letters)
# binary = makeString(35)

directory = "/Users/schilder/Desktop/X2K_Genetic_Algorithm/data/" # Change to your directory
#directory='/Users/maayanlab/PycharmProjects/X2K_Genetic_Algorithm/data/'

# binary = '11010101111110001100010001110001111'

# X2K Fitness Function for GA
def X2K_fitness(binary, fitness_method='simple'):
    ######################################
    #### Convert Binary to Parameters ####
    ######################################

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
    KINASE_topKinases = 20  # Should be consistent for every individual

    # CHEA OPTIONS ############
    TF_sort = {"10": "oddsratio", "01": "combined_score", "11": "rank", "00": "pvalue"}
    TF_species = {"10": "human", "01": "mouse", "11": "both", "00": "RESHUFFLE"}
    TF_databases = {"10": "chea", "01": "transfac", "11": "both", "00": "RESHUFFLE"}
    TF_background = {"10": "humanarchs4", "01": "humanarchs4", "11": "humanarchs4", "00": "RESHUFFLE"}
    TF_topTFs = {"10": 5, "01": 10, "11": 20, "00": "RESHUFFLE"}

    # G2N OPTIONS ############
    ## Databases
    all_PPI_databases = ["BIND", "BIOCARTA", "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEGG",
                         "MINT", "MIPS", "MURPHY", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL", "HUMAP"] # "BIOGRID", "KEA"
    PPI_dbs = []
    PPI_string = bit_dict["PPI_databases"]
    # ******* Turn off specific databases *******
    #ppi_list = list(PPI_string)
    #ppi_list[2] = '0'  # Turn off BIOGRID PPI
    #ppi_list[8] = '0'  # Turn off KEA PPI
    #PPI_string = "".join(ppi_list)
    #bit_dict["PPI_databases"] = PPI_string

    # *******************************************
    ### If no database selected, select at least one
    from random import randint
    while sum(map(int, PPI_string)) == 0:
        ppi_list = list(PPI_string)
        ppi_list[randint(0, (len(ppi_list)) - 1)] = '1'
        # ******* Turn off specific databases ******* AGAIN just in case it picked one of these
        #ppi_list[2] = '0'  # Turn off BIOGRID PPI
        #ppi_list[8] = '0'  # Turn off KEA PPI
        PPI_string = "".join(ppi_list)
        bit_dict["PPI_databases"] = PPI_string
    ### Generate list of G2N databases
    for ind, bit in enumerate(bit_dict["PPI_databases"]):
        if bit == "1":
            PPI_dbs.append(all_PPI_databases[ind])
    PPI_databases = ",".join(PPI_dbs)
    ## Path length
    PPI_pathLength = {"0": 1, "1": 2}

    ############ KEA OPTIONS ############
    KINASE_sort = {"10": "oddsratio", "01": "combined_score", "11": "rank", "00": "pvalue"}
    KINASE_interactome = {"10": "KP", "01": "P", "11": "RESHUFFLE", "00": "RESHUFFLE"}
    KINASE_background = {"10": "humanarchs4", "01": "humanarchs4", "11": "humanarchs4", "00": "RESHUFFLE"}

    parametersList = ["TF_sort", "TF_species", "TF_databases", "TF_background", "TF_topTFs", "PPI_pathLength",
                      "KINASE_sort", "KINASE_interactome", "KINASE_background"]
    from random import choice
    # RESHUFFLE any bits that aren't legitimate options
    for param in parametersList:
        bad_bits = [k for k, v in eval(param).items() if v == "RESHUFFLE"]
        good_bits = [k for k, v in eval(param).items() if v != "RESHUFFLE"]
        current_bits = bit_dict[param]
        if current_bits in bad_bits:
            bit_dict[param] = choice(good_bits)
            # print("Reshuffling...")



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
    g2n_parameters = ';'.join(["run",PPI_databases, str(PPI_pathLength[bit_dict["PPI_pathLength"]]) ]) + "\n"+allData+"messageComplete\n"
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
    kea_parameters = ';'.join(["run", KINASE_sort[bit_dict["KINASE_sort"]], KINASE_background[bit_dict["KINASE_background"]], KINASE_interactome[bit_dict["KINASE_interactome"]], str(KINASE_topKinases)  ]) + "\n"+allData2+"messageComplete\n"
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

    # LINCS_L1000 Kinase perturbations
    # kinase_perturbations_up = directory+"Validation/Perturbation_Data/L1000/LINCS_L1000_Kinase_Perturbations_up.txt"
    # kinase_perturbations_dn = directory+"Validation/Perturbation_Data/L1000/LINCS_L1000_Kinase_Perturbations_down.txt"


    # [1] CREEDS_ManualDrugs_DRH validation
    # Get list of predicted genes
    # import re
    # with open(directory+"output/kea_out.txt") as KEA_output,\
    #         open(directory+"Validation/Drug-Target_Data/Drug Repurposing Hub/Repurposing_Hub_export.txt") as DRH_validation, \
    #         open(directory+"Validation/Perturbation_Data/CREEDS/CREEDS_metadata_single_drug_perturbations-v1.0.csv") as CREEDS_manualDrugs_meta:
    #     validation_report = 0  # list of strings, to be summed together and divided by length of list to get fitness scores
    #     kea_predicted_kinases = []
    #
    #     for lineMeta in CREEDS_manualDrugs_meta:
    #         print(lineMeta)
    #
    #
    #     for KEA_line in KEA_output:
    #         KEA_rowName = KEA_line.split(",")[0].split(";")[1][:-7].upper() # Get row name
    #         KEA_genes = KEA_line.split(",")[1:]
    #         # Find the matching line in
    #         for DRH_line in DRH_validation:
    #             if re.split(r'\t+', DRH_line)[0].upper() == KEA_rowName:
    #                 print("Found it!")
    #                 validation_report = validation_report + 1
    #             else: print("Nope...")
    #     KEA_output.close()
    #     DRH_validation.close()


    # [2] LINCS1000 KINASE PERTURBATIONS: Presence/absence
    if fitness_method == 'simple':
        print("Calculating 'presence/absence' fitness...")
        with open (directory+"output/kea_out.txt") as KEA_output:
            lineCount = 0
            kinases_recovered = []
            kinases_missed = []

            for line in KEA_output:
                kinaseName = line.split(",")[0].split("_")[0]
                predictedKinases = line.split(",")[1:]
                if predictedKinases != []:
                    predictedKinases[-1] = predictedKinases[-1].strip()
                # Create lists of recovered and missed kinases
                if kinaseName in predictedKinases:
                    kinases_recovered.append(kinaseName)
                    #print("YAAAsss Queen!")
                else:
                    kinases_missed.append(kinaseName)
                    #print("Not so much...")
                lineCount += 1
            # Calculate individual's  fitness score
            if len(kinases_recovered) > 0:
                fitness_score = len(kinases_recovered) / lineCount * 100 # Number of kinases recovered over the total number of instances x 100
            else:
                fitness_score = 0
            KEA_output.close()
            print("Individual's fitness = " + str(fitness_score))




    # [3] LINCS1000 KINASE PERTURBATIONS: Rank-weighted fitness
    if fitness_method == 'rank-weighted':
        print("Calculating 'rank-weighted' fitness...")
        with open(directory + "output/kea_out.txt") as KEA_output:
            lineCount = 0
            kinases_recovered = []
            kinases_missed = []
            rankedScores = []

            for line in KEA_output:
                kinaseName = line.split(",")[0].split("_")[0]
                predictedKinases = line.split(",")[1:]
                # As long as the kinase list is not blank, get rid of \n in last item
                if predictedKinases != []:
                    predictedKinases[-1] = predictedKinases[-1].strip()
                # Create lists of recovered and missed kinases
                if kinaseName in predictedKinases:
                    rankedScores.append( (len(predictedKinases) - predictedKinases.index(kinaseName)) / len(predictedKinases))
                    kinases_recovered.append(kinaseName)
                else:
                    rankedScores.append(0)
                    kinases_missed.append(kinaseName)
                lineCount += 1
            # Calculate individual's  fitness score
            if len(kinases_recovered) > 0:
                temp_score = sum(rankedScores) / len(rankedScores)
                fitness_score = temp_score / KEA_topKinases

            else:
                fitness_score = 0
            KEA_output.close()
            print("Individual's fitness = " + str(fitness_score))

    if fitness_method == 'LINCS_L1000 + DrugRepurposingHub (adjustedScore)':
        print("Calculating"+fitness_method+"fitness...")
        with open(directory + "output/kea_out.txt") as KEA_output:
            lineCount = 0; kinases_recovered = []; kinases_missed = []; adjustedScores = []
            KEA_out = KEA_output.readlines()
            #line = "CPC005_PC3_24H-indirubin-10.0_[CDK1|CDK5|GSK3A]_up		USP15	DUSP6	YIPF	IDUA	CHRNA5	SYNRG	AURKA	IFT122	PHKA	GAL	PTMA	IDS	REEP4	HMOX	NPHP4	RREB	ATP5SL"
            for line in KEA_out:
                kinaseTargets = line.split("_")[3].strip("[").strip("]").split("|")
                predictedKinases = line.split(",")[1:]
                # As long as the kinase list is not blank, get rid of \n in last item
                if predictedKinases != []:
                    predictedKinases[-1] = predictedKinases[-1].strip()
                # Create lists of recovered and missed kinases
                kinases_recov = []; kinase_miss = []
                for target in kinaseTargets:
                    if target in predictedKinases:
                        kinases_recov.append(target)
                    else:
                        kinase_miss.append(target)
                # The adjustedScore is the percentage of targets that were captured by KEA for a given experiment (since some drugs have multiple targets and this creates more chances for a target to be discovered)
                if len(kinases_recov)>0:
                    adjustedScores.append( len(kinases_recov)/ len(kinaseTargets))
                    kinases_recovered.append(",".join(kinases_recov))
                else:
                    adjustedScores.append(0) # Make sure you append a 0 if it missed all kinases
                    kinases_missed.append(",".join(kinase_miss))
                lineCount += 1
            # Calculate individual's  fitness score
            if len(kinases_recovered) > 0:
                fitness_score = sum(adjustedScores) / len(adjustedScores)*100  # Number of kinases recovered over the total number of instances x 100
            ## Workaround to limit PPI size:
            # elif average_PPI_size > 100:
            #     fitness_score = 0
            else:
                fitness_score = 0
            KEA_output.close()
            print("Individual's fitness = " + str(fitness_score))



    if fitness_method == 'LINCS_L1000 + KINOMEscan (adjustedScore)':
        print("Calculating"+fitness_method+"fitness...")
        import Python_scripts.rbo as rbo

        with open(directory + "output/kea_out.txt") as KEA_output:
            lineCount = 0; kinases_recovered = []; kinases_missed = []; adjustedScores = []
            KEA_out = KEA_output.readlines()
            #line = "LJP009_PC3_24H-sunitinib-1.11_[MAPK7|TNK2|MATK|CDK5|TSSK1B|MKNK2|PIP5K1A|MAP2K6|NTRK3|MAP3K4|PLK2|PTK6|MYO3B|MAPK10|PRKCQ|NEK7|CIT|MKNK1|TIE1|RIOK3|MARK4|BRSK1|PDPK1|CDK4|MAP3K9|MET|SIK1|EPHB4|MYO3A|DDR2|AKT2|DCLK2|IGF1R|BMPR1B|EPHA7|MAPK9|PAK6|PIM3|RPS6KA6|DYRK1B|EGFR|PAK4|BTK|EPHA3|FGFR4|RIPK4|SRC|DDR1|SgK110|LTK|AURKA|MAP2K3|CDK18|CAMKK2|CAMK2B|NEK2|MAP3K15|MAP3K11|CDPK1|PKN2|DSTYK|STK35|EPHA5|JAK3|MARK1|CDK17|BRSK2|CDKL2|FER|WEE1|ABL2|EPHB6|CAMK1|STK38L|EPHA6|FES|IRAK3|CSNK1G1|CSNK2A1|CAMK4|TAOK1|HCK|MGC42105|TLK1|PKN1|MAP2K4|CAMK2G|DYRK2|TNK1|TYK2|RPS6KA1|EIF2AK2|PAK7|SNRK|LATS1|IKBKE|ERN1|FGFR3|NTRK2|RPS6KA3|SIK2|BMPR2|FGFR2|CSNK1A1L|FRK|OXSR1|FGFR1|FYN|JAK1|CHUK|CAMK1D|HUNK|INSR|MUSK|EPHB1|ICK|LATS2|ROCK1|CAMK1G|PTK2|INSRR|RPS6KA2|CAMK2D|RPS6KA4|CAMKK1|JAK2|MARK3|STK38|PRPF4B|AURKB|PRKD2|DCLK1|RIPK1|MELK|MST4|CDK7|TLK2|ANKK1|MARK2|PRKD1|CHEK1|GRK1|RPS6KA5|STK25|MYLK|PRKD3|FGR|LYN|CDK14|SRPK1|STK16|ABL1|CSNK1G3|LCK|AURKC|MAP3K3|SGK3|TAOK3|DYRK1A|MAST1|SBK1|PLK4|SRPK2|EIF2AK4|GRK7|MAP4K3|ALK|CSNK2A2|HIPK4|DAPK2|CDK16|NUAK2|GRK4|MAP4K4|ROCK2|STK39|MAP2K1|DAPK1|TBK1|YES1|CSNK1G2|DCLK3|STK17B|MAP2K2|MAP3K12|NTRK1|CSNK1A1|MAP3K13|LRRK2|MAP3K7|PRKAA2|pknB|PTK2B|CAMK2A|IRAK4|BLK|STK24|TTK|SRPK3|MAP3K2|STK3|SLK|HIPK1|FLT4|MYLK2|RIOK2|TYRO3|NUAK1|RPS6KB1|MAP2K5|ULK3|HIPK3|MAP4K5|PIP4K2B|STK11|RIOK1|MAP4K2|HIPK2|CLK4|MINK1|MERTK|TNIK|MYLK3|ULK1|CLK1|DAPK3|CLK2|KIT|GAK|PRKAA1|STK10|STK4|STK33|YSK4|MAP4K1|PAK3|CSNK1D|MYLK4|IRAK1|CSNK1E|ITK|RET|ULK2|AAK1|CHEK2|AXL|PHKG2|BMP2K|PHKG1|CSF1R|FLT1|FLT3|KDR|STK17A|PDGFRA|PDGFRB]_up		HOMER	DKK	PAWR	TM4SF	QKI	CXCL13	NUP98	ID	UBE2S	SPINT2	IL8	MT1H	IL6ST	GAS6	NNMT	CHIC2	SERPINE	AURKA	CXCL2	NEDD9	RELB	LAPTM4B	SOCS2	DNAJB6	TBL1X	MAGEA6	PDLIM5	ARHGEF12	POSTN	PLK	DAPK	FAM46A	RFTN	PLK2	HDGFRP3	ADAM	COBL	LEPREL	KCTD5	EIF2S3	SERPINB9	NR4A3	NR4A2	RPS26	RAB15	DLC	HIST1H2BG	CRYAB	FOLR	GPM6B	TYRP	HIST1H2BK	IRS2	THBS	NR3C	NTS	PTPLA	PPM1H	MICAL2	BAG2	MFAP3L"

            import re
            for line in KEA_out:
                kinaseTargets = re.split(r'\t+', line)[0].split("_")[3].strip("[").strip("]").split("|")
                predictedKinases = re.split(r'\t+', line)[1:]
                # As long as the kinase list is not blank, get rid of \n in last item
                if predictedKinases != []:
                    predictedKinases[-1] = predictedKinases[-1].strip()
                # Create lists of recovered and missed kinases
                kinases_recov = []; kinase_miss = []
                ## Need to take into the account:
                ### 1. The ranking of the KINOMEscan target
                ### 2. The ranking of the KEA predicted kinase
                for i,target in enumerate(kinaseTargets):
                    if target in predictedKinases:
                        kinases_recov.append(target)
                    else:
                        kinase_miss.append(target)
                # The adjustedScore is the percentage of targets that were captured by KEA for a given experiment (since some drugs have multiple targets and this creates more chances for a target to be discovered)
                if len(kinases_recov)>0:
                    adjustedScores.append( len(kinases_recov)/ len(kinaseTargets))
                    kinases_recovered.append(",".join(kinases_recov))
                else:
                    adjustedScores.append(0) # Make sure you append a 0 if it missed all kinases
                    kinases_missed.append(",".join(kinase_miss))
                lineCount += 1
            # Calculate individual's  fitness score
            if len(kinases_recovered) > 0:
                fitness_score = sum(adjustedScores) / len(adjustedScores)*100  # Number of kinases recovered over the total number of instances x 100
            else:
                fitness_score = 0
            KEA_output.close()
            print("Individual's fitness = " + str(fitness_score))

    if fitness_method == 'L1000_DRH - RankBasedOrder':
        print("Calculating"+fitness_method+"fitness...")
        from Python_scripts.rbo import rbo

        with open(directory + "output/kea_out.txt") as KEA_output:
            rboScores=[]
            KEA_out = KEA_output.readlines()
            #line = "CPC005_PC3_24H-indirubin-10.0_[CDK1|CDK5|GSK3A]_up		USP15	DUSP6	YIPF	IDUA	CHRNA5	SYNRG	AURKA	IFT122	PHKA	GAL	PTMA	IDS	REEP4	HMOX	NPHP4	RREB	ATP5SL"
            for line in KEA_out:
                # Set up your two lists (targets and predicted)
                kinaseTargets = line.split("_")[3].strip("[").strip("]").split("|")
                predictedKinases = line.split(",")[1:]
                # As long as the kinase list is not blank, get rid of \n in last item
                if predictedKinases != []:
                    predictedKinases[-1] = predictedKinases[-1].strip()

                # Get the simple intersection between the two lists
                ## intersection = [x for x in kinaseTargets if x in predictedKinases]
                ## percentOverlap = len(intersection)/len(kinaseTargets)*100
                if len(predictedKinases)>0:
                    # Conduct RankedBasedOrder statistic
                    rboResults = rbo(kinaseTargets, predictedKinases, .9)
                    rboScores.append(rboResults['res'])
                else: rboScores.append(0)
                # Get the average RBO score
                fitness_score = sum(rboScores) / len(rboScores)

            print("Individual's fitness = " + str(fitness_score))

    return fitness_score, average_PPI_size



