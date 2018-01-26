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
# binary = makeString(28)

#directory = "/Users/schilder/Desktop/x2k/data/" # Change to your directory
directory='/Users/maayanlab/PycharmProjects/X2K_Genetic_Algorithm/data/'

# binary = '111011001000000000000100100'

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
                "PPI_databases": binary[10:20],
                "PPI_pathLength": binary[20:21],
                "KINASE_sort": binary[21:23],
                "KINASE_interactome": binary[23:25],
                "KINASE_background": binary[25:27]
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
    all_PPI_databases = ["BIND", "BIOCARTA", "BIOGRID", "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEA", "KEGG",
                         "MINT", "MIPS", "MURPHY", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL", "HUMAP"]
    PPI_dbs = []
    PPI_string = bit_dict["PPI_databases"]
    # ******* Turn off specific databases *******
    ppi_list = list(PPI_string)
    ppi_list[2] = '0'  # Turn off BIOGRID PPI
    ppi_list[8] = '0'  # Turn off KEA PPI
    PPI_string = "".join(ppi_list)
    bit_dict["PPI_databases"] = PPI_string

    # *******************************************
    ### If no database selected, select at least one
    from random import randint
    while sum(map(int, PPI_string)) == 0:
        ppi_list = list(PPI_string)
        ppi_list[randint(0, (len(ppi_list)) - 1)] = '1'
        # ******* Turn off specific databases ******* AGAIN just in case it picked one of these
        ppi_list[2] = '0'  # Turn off BIOGRID PPI
        ppi_list[8] = '0'  # Turn off KEA PPI
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
            else:
                fitness_score = 0
            KEA_output.close()
            print("Individual's fitness = " + str(fitness_score))


    if fitness_method == 'LINCS_L1000 + KINOMEscan (adjustedScore)':
        print("Calculating"+fitness_method+"fitness...")
        import Python_scripts.rbo_function
        with oppen(directory + "output/kea_out.txt") as KEA_output:
            lineCount = 0; kinases_recovered = []; kinases_missed = []; adjustedScores = []
            KEA_out = KEA_output.readlines()
            #line = "LJP009_MCF7_24H-JNK-IN-5A-10_[ACVR1B|ACVRL1|CABC1|ADCK4|AKT2|ALK|AURKA|BMPR1A|BMX|BRSK1|BRSK2|CAMK4|CDC2L2|CDK5|CDK8|CIT|CSNK1A1|CSNK1D|MATK|DCLK1|DCLK2|DCLK3|DDR2|STK17A|STK17B|DYRK1A|DYRK2|EIF2AK1|EPHA2|EPHA3|EPHA4|EPHB4|ERBB3|ERBB4|FES|FGFR2|FGR|GRK1|GRK4|HIPK2|HIPK4|CHUK|IKBKE|INSR|INSRR|IRAK1|IRAK3|IRAK4|JAK2|LIMK2|MAP3K15|MAP3K2|MAP4K2|MAPKAPK5|MARK1|MARK2|MAST1|MAP2K2|MAP2K5|MAP2K6|MERTK|MAP2K7|MKNK1|MYLK3|MAP3K10|CDC42BPA|MTOR|MYLK2|STK38|NEK1|NEK11|NEK3|MGC42105|OXSR1|MAPK14|MAPK11|PAK1|PAK4|PAK7|CDK16|MAL13P1.279|PHKG1|PIK3C2B|PIK3CB|PIK3CD|PI4KB|PIM1|PIM3|PIP5K1A|PIP5K1C|PIP4K2B|PIP4K2C|pknB|PLK2|PLK3|PRKCD|PRKCE|PRKD1|EIF2AK2|BRAF|PRKX|ABL1|PRPF4B|RAF1|RIOK1|RIOK3|DSTYK|ROCK1|ROCK2|SBK1|SNRK|SRPK1|SRPK2|SRPK3|STK39|TAOK1|TAOK2|TGFBR1|TGFBR2|TEK|TLK2|RPS6KA3|TNIK|NTRK1|NTRK2|TRPM6|TSSK1B|TXK|TYK2|ULK1|VRK2|YSK4|BUB1|GSG2|SGK1|WNK3|FGFR1|MYLK|PIK3CG|SYK|EPHB3|JAK3|MAP2K1|MAP3K11|MYLK4|CDK4|STK38L|CDK14|ULK2|CDK3|CDKL2|CSF1R|MARK4|MUSK|PAK2|STK32C|BMPR1B|MAPK7|MAP3K3|RIPK4|RPS6KA4|KDR|FGFR3|CDK11B|CDK2|CSK|HUNK|MINK1|MST1R|CDK18|NUAK2|TIE1|AKT1|MAP3K5|AURKB|CAMK2G|CAMKK2|CSNK2A2|EPHA6|LRRK2|MAP2K4|STK4|PAK3|PKN1|ULK3|MAP3K6|CAMK2D|LATS1|MAPKAPK2|CDK17|RET|STK25|PIK3CA|CAMK2A|FER|HCK|MKNK2|NEK6|SGK3|SIK1|AAK1|BLK|CAMKK1|DAPK2|EPHB1|PTK2|FLT4|HIPK1|STK3|RIOK2|SgK110|TEC|TTK|BTK|RPS6KA1|CAMK1G|CAMK2B|MAPK4|ERN1|RPS6KA5|LTK|ROS1|RPS6KB1|CSNK2A1|DYRK1B|EPHA5|MAPK1|MAP4K1|ITK|MARK3|FLT3|PDGFRA|STK16|WNK1|AXL|CDK19|DAPK3|GRK7|MAP2K3|MYO3B|PAK6|PRKD2|PRKG1|TBK1|AURKC|CDK13|CLK2|ERBB2|FLT1|KIT|PIK3C2G|EGFR|PRKCH|RIPK1|TAOK3|ZAK|ZAP70|ABL2|RPS6KA2|BMP2K|MET|MAPK6|MAP3K4|PRKCQ|PRKG2|STK33|NTRK3|CAMK1D|DAPK1|EPHA1|EPHB6|LCK|JAK1|MELK|PRKAA2|BMPR2|CHEK2|CLK4|EPHA7|NEK7|PRKD3|SRMS|CDK9|DMPK|MAP3K13|PHKG2|TESK1|STK32B|AKT3|CDK7|CDKL3|CHEK1|EIF2AK4|MAPK12|PLK4|SLK|MAP3K7|CDC42BPG|ICK|STK11|PDPK1|CDPK1|PKN2|PLK1|SRC|STK32A|PRKAA1|CDC42BPB|STK24|PRKACB|MAP3K12|EPHA8|MAP4K3|MST4|NEK5|PTK2B|ANKK1|LYN|PRKACA|SIK2|YES1|NUAK1|PKMYT1|TNK2|IGF1R|LIMK1|NEK4|TLK1|TYRO3|EPHB2|LATS2|MAPK3|TNK1|CAMK1|RPS6KA6|CDKL1|CLK3|GSK3A|IKBKB|MAP3K9|CASK|DDR1|CDKL5|CSNK1E|NEK9|WEE2|ACVR2A|CLK1|CDK15|PIM2|CSNK1G3|ACVR2B|GAK|STK10|CSNK1A1L|GSK3B|MAP4K4|STK36|MAP4K5|NEK2|NLK|SIK3|ACVR1|MAPK15|STK35|WEE1|FYN|PRKCI|MAPK8|MAP3K1|FRK|RIPK2|CSNK1G2|MAK|TNNI3K|CSNK1G1|FGFR4|PTK6|MAPK13|PDGFRB|HIPK3|MYO3A|MAPK9|MAPK10]		MREG	GINS	C5ORF13	UBE2C	DUSP2	CADM	UBE2A	C7ORF68	CGRRF	PDGFRL	IFI44L	GNPDA	KRT19	C7ORF23	MGP	YME1L	CDK	EAPP	SRSF7	PCNA	TMED	BTG3	SERPINA	IFITM	TOP2A	UCP2	CCNF	GLS	SLC5A3	SCUBE2	NRIP	CTSK	PACSIN3	DNAJB9	CYP1B	CAPN	HES	FIS	EED	RRM2	DDX17	IGFBP5	SMAD	QPCT	AKR1C3	RNASE4	MAFB	ESR	TPD52L	RSL1D	ASPN	RYBP	STAG2	ARMCX2	PIK3R	CCNE2	RNF128	NET	SCGB2A	DTL	KRAS	SELT	SLC24A3	TIMM9	FEN	CPB	CDCA4	DPT	LGALS3BP	WFDC2	COPS7A	PLAT	KIAA	HSD17B12	GATA3	IFIT	FAM134B	G3BP	PCBP2	GPNMB	LRIG	TSPAN6	NUSAP	DEPTOR	HADH"
            import re
            for line in KEA_out:
                kinaseTargets = re.split(r'\t+', line)[0].split("_")[-1].strip("[").strip("]").split("|")
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

    return fitness_score, average_PPI_size



