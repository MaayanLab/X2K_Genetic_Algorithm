#############################################################################################
##################### X2K: Pipeline & Fitness Function #####################
#############################################################################################
directory = "/Users/schilder/Desktop/x2k/data/" # Change to your directory

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
# binary = '1001100010110111000111010000'

# X2K Fitness Function for GA
def X2K_fitness(binary):
    ######################################
    #### Convert Binary to Parameters ####
    ######################################

    # Readjust bit parameters if length of binary string changes
    bit_dict = {"CHEA_sort": binary[0:2],
                "CHEA_species": binary[2:4],
                "CHEA_databases": binary[4:6],
                "CHEA_background": binary[6:8],
                "CHEA_topTFs": binary[8:10],
                "G2N_databases": binary[10:20],
                "G2N_pathLength": binary[20:21],
                "KEA_sort": binary[21:23],
                "KEA_interactome": binary[23:25],
                "KEA_background": binary[25:27]
                }
    KEA_topKinases = 50  # Should be consistent for every individual

    # CHEA OPTIONS ############
    ## sort
    CHEA_sort = {"10": "oddsratio", "01": "RESHUFFLE", "11": "rank", "00": "RESHUFFLE"}
    CHEA_species = {"10": "human", "01": "mouse", "11": "both", "00": "RESHUFFLE"}
    CHEA_databases = {"10": "chea", "01": "transfac", "11": "both", "00": "RESHUFFLE"}
    CHEA_background = {"10": "humanarchs4", "01": "mousearchs4", "11": "hg133", "00": "RESHUFFLE"}
    CHEA_topTFs = {"10": 5, "01": 10, "11": "RESHUFFLE", "00": "RESHUFFLE"}

    # G2N OPTIONS ############
    ## Databases
    all_G2N_databases = ["BIND", "BIOCARTA", "BIOGRID", "DIP", "FIGEYS", "HPRD", "INNATEDB", "INTACT", "KEA", "KEGG",
                         "MINT", "MIPS", "MURPHY", "PDZBASE", "PPID", "PREDICTEDPPI", "SNAVI", "STELZL", "VIDAL",
                         "HUMAP"]
    G2N_dbs = []
    G2N_string = bit_dict["G2N_databases"]
    ### If no database selected, select at least one
    from random import randint
    if sum(map(int, G2N_string)) == 0:
        g2n_list = list(G2N_string)
        g2n_list[randint(0, (len(g2n_list)) - 1)] = '1'
        bit_dict["G2N_databases"] = "".join(g2n_list)
    ### Generate list of G2N databases
    for ind, bit in enumerate(G2N_string):
        if bit == "1":
            G2N_dbs.append(all_G2N_databases[ind])
    G2N_databases = ",".join(G2N_dbs)
    ## Path length
    G2N_pathLength = {"0": 1, "1": 2}

    ############ KEA OPTIONS ############
    KEA_sort = {"10": "oddsratio", "01": "combined_score", "11": "rank", "00": "RESHUFFLE"}
    KEA_interactome = {"10": "KP", "01": "P", "11": "RESHUFFLE", "00": "RESHUFFLE"}
    KEA_background = {"10": "humanarchs4", "01": "mousearchs4", "11": "hg133", "00": "RESHUFFLE"}

    parametersList = ["CHEA_sort", "CHEA_species", "CHEA_databases", "CHEA_background", "CHEA_topTFs", "G2N_pathLength",
                      "KEA_sort", "KEA_interactome", "KEA_background"]  # g2n_databases
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
    chea_parameters = ';'.join( ["run", CHEA_sort[bit_dict["CHEA_sort"]], CHEA_species[bit_dict["CHEA_species"]], CHEA_databases[bit_dict["CHEA_databases"]], CHEA_background[bit_dict["CHEA_background"]], str(CHEA_topTFs[bit_dict["CHEA_topTFs"]]) ] ) #"run,oddsratio,both,chea,humanarchs4,20")
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
    g2n_parameters = ';'.join(["run",G2N_databases, str(G2N_pathLength[bit_dict["G2N_pathLength"]]) ]) + "\n"+allData+"messageComplete\n"
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
    kea_parameters = ';'.join(["run", KEA_sort[bit_dict["KEA_sort"]], KEA_background[bit_dict["KEA_background"]], KEA_interactome[bit_dict["KEA_interactome"]], str(KEA_topKinases)  ]) + "\n"+allData2+"messageComplete\n"
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
    # print("Calculating 'rank-weighted' fitness...")
    # with open(directory + "output/kea_out.txt") as KEA_output:
    #     lineCount = 0
    #     kinases_recovered = []
    #     kinases_missed = []
    #     rankedScores = []
    #
    #     for line in KEA_output:
    #         kinaseName = line.split(",")[0].split("_")[0]
    #         predictedKinases = line.split(",")[1:]
    #         # As long as the kinase list is not blank, get rid of \n in last item
    #         if predictedKinases != []:
    #             predictedKinases[-1] = predictedKinases[-1].strip()
    #         # Create lists of recovered and missed kinases
    #         if kinaseName in predictedKinases:
    #             rankedScores.append( (len(predictedKinases) - predictedKinases.index(kinaseName)) / len(predictedKinases))
    #             kinases_recovered.append(kinaseName)
    #         else:
    #             rankedScores.append(0)
    #             kinases_missed.append(kinaseName)
    #         lineCount += 1
    #     # Calculate individual's  fitness score
    #     if len(kinases_recovered) > 0:
    #         temp_score = sum(rankedScores) / len(rankedScores)
    #         fitness_score = temp_score / KEA_topKinases
    #
    #     else:
    #         fitness_score = 0
    #     KEA_output.close()
    #
    #     print("Individual's fitness = " + str(fitness_score))

    return fitness_score


