## ********** FORMATTING DATABASE GMT FILES FOR X2K ********** ##

import pandas as pd
import os


os.getcwd()
os.chdir("../")
import pandas as pd
# Progress bar for for loop
def percentComplete(loopLength, iteration):
    percent = round((iteration + 1) / loopLength * 100, 2)
    if percent % 5 == 0:
        print(str(percent) + "% complete")


# Add all synonyms (both within and between mouse/human) to geneList
def addSynonyms(geneList):
    geneSyn = pd.read_table("X2K_Databases/General_Resources/Moshe_mapping/mappingFile_2017.txt", header=None)
    intersect = list(set(geneSyn[0]).intersection(geneList))
    for gene in intersect:
        syn = geneSyn[geneSyn[0] == gene][1].values.tolist()
        for s in syn:
            if s not in geneList:
                print("Adding synonym: " + s)
                geneList.append(s)
    print("New geneList length = " + str(len(geneList)))
    return geneList


# Get all TFs
## Human TFs from review by Lambert et.al. 2018:
Lambert_TF = pd.read_excel("X2K_Databases/General_Resources/Lambert-et-al-2018_Human-TFs/mmc2.xlsx",
                           sheet_name="Table S1. Related to Figure 1B", skiprows=1)
Lambert_TF.head()
TFs = list(Lambert_TF[Lambert_TF.iloc[:, 3] == "Yes"].Name.unique())
##  Mouse TFs from mTFkb: # http://sunlab.cpy.cuhk.edu.hk/mTFkb/download.php
mTF = pd.read_excel("X2K_Databases/General_Resources/mTFkb/mTFkb database.xlsx", sheet_name="TFdb")
for tf in mTF.Symbol.unique():
    if tf not in TFs:
        TFs.append(tf)
    if tf.upper() not in TFs:
        TFs.append(tf.upper())
# Add synonyms
TFs = addSynonyms(TFs)
pd.Series(TFs[1]).to_csv("X2K_Databases/General_Resources/compiled-TFs_Mouse-Human.csv", index=False)

# Get all Kinase Names, Groups, and Families
# ## KEA
# with open("X2K_Genetic_Algorithm/data/KEA/kea_ranks.txt") as KEA_kinases:
#     KINASES = []
#     for line in KEA_kinases:
#         KINASES.append(line.split("\t")[0])
# ## Kinase.com
# kinaseCom_HsapUpdated = pd.read_excel("X2K_Databases/General_Resources/Kinase.com/Kinome_Hsap_updated.xls")
# kinaseCom_Hsap = pd.read_excel("X2K_Databases/General_Resources/Kinase.com/Kinome_Hsap.xls")
# kinaseCom_Mmus = pd.read_excel("X2K_Databases/General_Resources/Kinase.com/Kinome_Mmus.xls")
# kinaseCom_Mmus = kinaseCom_Mmus.rename(index=str, columns={"Gene Name": "Name"})
# kinaseCom_Mmus.columns
#
# kinaseGroups = pd.concat([ kinaseCom_HsapUpdated[["Name","Group","Family"]],\
#           kinaseCom_Hsap[["Name","Group","Family"]],\
#           kinaseCom_Mmus[["Name","Group","Family"]] ])
# kinaseGroups.to_csv("X2K_Databases/General_Resources/Kinase.com/Kinase_Groups_&_Families.csv", header=True)
# kinaseCom = list(set(kinaseCom_HsapUpdated.Name) | set(kinaseCom_Hsap.Name) | set(kinaseCom_Mmus.Name))
# ## Merge lists
# for k in kinaseCom:
#     if k not in KINASES:
#         KINASES.append(k)

## Maxime's Kinase Groups & Families (from Kinome.com)
kinaseGroups = pd.read_table("X2K_Databases/General_Resources/Kinase_Families_Maxime/kinases_fam.tsv", header=None, names=["Name","Group","Family"], index_col=False)
### Make kinase group/family dictionary
kgroupDict = dict(zip(kinaseGroups.Name, zip(kinaseGroups.Group, kinaseGroups.Family)))

def addSynonymsToKinaseGroups(kgroupDict, kinaseList):
    geneSyn = pd.read_table("X2K_Databases/General_Resources/Moshe_mapping/mappingFile_2017.txt", header=None)
    for i, kinase in enumerate(kinaseList):
        if kinase in geneSyn[0] and kinase in kgroupDict.keys():
            syns = list(geneSyn[geneSyn[0] == kinase][1].values)
            for s in syns:
                print("Adding to kgroupDict: "+s)
                kgroupDict[s] = kgroupDict[kinase]
    return kgroupDict
KEA = pd.read_table("X2K_Genetic_Algorithm/data/KEA/kea_ranks.txt", header=None)
kgroupDict = addSynonymsToKinaseGroups(kgroupDict, KEA[0].tolist())
KINASES = kinaseGroups.Name.tolist()
## Add synonyms
KINASES = addSynonyms(KINASES)
pd.Series(KINASES).to_csv("X2K_Databases/General_Resources/compiled-Kinases_Mouse-Human.csv", index=False)
### Add Group and Family info to name
def getKinaseGroupAndFamily(kinaseList):
    geneSyn = pd.read_table("X2K_Databases/General_Resources/Moshe_mapping/mappingFile_2017.txt", header=None, names=["Gene","GeneSyn"])
    kinaseGroups = pd.read_csv("X2K_Databases/General_Resources/Kinase.com/Kinase_Groups_&_Families.csv")
    kgroupDict = dict(zip(kinaseGroups.Name, zip(kinaseGroups.Group, kinaseGroups.Family)))
    for k in kinaseList:
        syns=[k]
        if k in list(geneSyn.Gene):
            newSyns = list( geneSyn[geneSyn.Gene == k].GeneSyn.unique() )
            for s in newSyns:
                syns.append(s)
        overlap = list( set(syns).intersection(set(kgroupDict.keys())) )
        if len(overlap) >0:
            results = kgroupDict[overlap[0]]
        else:
            results = ["NA","NA"]
    return results



# Summary table
def getDatasetsSummary(writeExcel=False):

    def extract_TF_KINASE_stats(dataset_path):
        with open(dataset_path) as DATA:
            data = DATA.readlines()
        unique_TF_Ki = [];
        targets_substrates = [];
        targets_substrates_len = []
        if dataset_path.endswith(tuple([".gmt", ".txt"])):
            for line in data:
                # Get TF/Kinase
                TF_Ki = line.split("\t")[0].split("_")[0]
                if TF_Ki not in unique_TF_Ki:
                    unique_TF_Ki.append(TF_Ki)
                # Get DEGs
                DEGs = [gene.replace(',1.0', '') for gene in line.split("\t")[2:]]
                targets_substrates_len.append(len(DEGs))
                for gene in DEGs:
                    if gene not in targets_substrates:
                        targets_substrates.append(gene)
            num_TF_Kinase_Sets = len(data)
            num_unique_TF_Ki = len(unique_TF_Ki)
            num_targets_substrates = len(targets_substrates)
            avg_set_size = round( sum(targets_substrates_len) / len(targets_substrates_len), 2)
        return num_TF_Kinase_Sets, num_unique_TF_Ki, num_targets_substrates, avg_set_size

    def extract_PPI_stats(dataset_path):
        if dataset_path.endswith("sig"):
            data = pd.read_table(dataset_path, header=None, delim_whitespace=True)
            totalInteractions = data.shape[0]
            uniqueInteractions = [];
            interactionsPerProtein = []
            for i, row in data.iterrows():
                gene = row.iloc[0]
                interactor = row.iloc[5]
                if gene + "@" + interactor not in uniqueInteractions and interactor + "@" + gene not in uniqueInteractions:
                    uniqueInteractions.append(gene + "@" + interactor)
            uniqueProteins = list(data[0].unique())
            for prot in uniqueProteins:
                prot_sub = data[data[0] == prot]
                interactionsPerProtein.append(prot_sub.__len__())
            avgInteractionsPerProtein = round( sum(interactionsPerProtein) / len(interactionsPerProtein),2 )
            return totalInteractions, len(uniqueInteractions), len(uniqueProteins), avgInteractionsPerProtein

    # Make summary table
    import os
    import pandas as pd
    import numpy as np
    root = "/Users/schilder/Desktop/X2K_Databases"
    pathList = [root + "/TF", root + "/PPI", root + "/KINASE"]
    TF_summaryTable = pd.DataFrame(); KINASE_summaryTable = pd.DataFrame(); PPI_summaryTable = pd.DataFrame()
    for path in pathList:
        dataType = path.split("/")[-1]
        # Get list of only folders
        foldersOnly = [x for x in os.listdir(path) if os.path.isdir(os.path.join(path, x))]
        for folder in foldersOnly:
            dir_contents = os.listdir(path + "/" + folder)
            for file in dir_contents:
                if file.endswith(tuple([".gmt", ".sig", ".txt"])):
                    print(os.path.join(path, folder, file))
                    try:
                        # TF table
                        if dataType == "TF":
                            num_TF_Kinase_Sets, num_unique_TF_Ki, num_targets_substrates, avg_set_size = extract_TF_KINASE_stats(os.path.join(path, folder, file))
                            cols = ['Data Type', 'Database Name', 'File Name', 'TF Sets', 'Unique TFs','Unique Targets', 'Average Targets per TF']
                            TF_summaryTable = TF_summaryTable.append( \
                                pd.DataFrame( np.column_stack([dataType,folder,file,num_TF_Kinase_Sets,num_unique_TF_Ki, num_targets_substrates, avg_set_size]),columns=cols ) )
                        # KINASE table
                        elif dataType == "KINASE":
                            num_TF_Kinase_Sets, num_unique_TF_Ki, num_targets_substrates, avg_set_size = extract_TF_KINASE_stats(os.path.join(path, folder, file))
                            cols=['Data Type','Database Name','File Name','Kinase Sets', 'Unique Kinases','Unique Substrates', 'Average Substrates per Kinase']
                            KINASE_summaryTable = KINASE_summaryTable.append( \
                                pd.DataFrame(np.column_stack([dataType,folder,file,num_TF_Kinase_Sets,num_unique_TF_Ki,num_targets_substrates,avg_set_size]),columns=cols))
                        # PPI table
                        elif dataType == "PPI":
                            totalInteractions, uniqueInteractions, uniqueProteins, avgInteractionsPerProtein = extract_PPI_stats(os.path.join(path, folder, file))
                            cols=['Data Type','Database Name','File Name','Total Interactions','Unique Interactions','Unique Proteins', 'Average Interactions Per Protein']
                            PPI_summaryTable = PPI_summaryTable.append( \
                                pd.DataFrame(np.column_stack([dataType,folder,file,totalInteractions,uniqueInteractions,uniqueProteins,avgInteractionsPerProtein]),columns=cols))
                    except:
                        print("Could not summarize:"+ os.path.join(path, folder, file))

    if writeExcel != False:
        print("Writing table to excel file...")
        frames = [TF_summaryTable, PPI_summaryTable, KINASE_summaryTable]
        start_row = 1
        writer=  pd.ExcelWriter(writeExcel)
        for df in frames:  # assuming they're already DataFrames
            df.to_excel(writer, 'Sheet1', startrow=start_row, index=False)
            start_row += len(df) + 1  # add a row for the column header?
        writer.save()

    return TF_summaryTable, KINASE_summaryTable, PPI_summaryTable


TF_summaryTable, PPI_summaryTable, KINASE_summaryTable= getDatasetsSummary(writeExcel="X2K_Databases/X2K_Database_Summary.xlsx")

# HomoloGene
# XXXXX DOES NOT HAVE WITHIN-SPECIES GENE SYNONYMS. USE MOSHE'S MOUSE-HUMAN HOMLOG MAPPING INSTEAD XXXXX
taxid = pd.read_table("X2K_Databases/General_Resources/HomoloGene/taxid_taxname.dms", header=None)
homolo = pd.read_table("X2K_Databases/General_Resources/HomoloGene/homologene.data", '\t')
homolo.columns = ["HID", "TaxonID", "GeneID", "GeneSymbol", "ProteinGI", "ProteinAC"]
homolo.head()
human = homolo[homolo.TaxonID == 9606]
mouse = homolo[homolo.TaxonID == 10090]


# If the gene is from mouse, convert to human version
def geneConverter(gene, species='?'):
    def recordMissingGene(gene):
        if gene not in missingGenes:
            missingGenes.append(gene)
        print("Could not find gene: (" + gene + ", " + str(len(missingGenes)) + ")")

    missingGenes = []
    # If it's a human gene, use original
    if gene in list(human.GeneSymbol) or species.upper() == 'HUMAN':
        newGene = gene
    # If it's a mouse gene, convert to human
    elif gene in list(mouse.GeneSymbol) and species.upper() == 'MOUSE':
        HID = mouse[mouse.GeneSymbol == gene].HID.values[0]
        if HID in human[human.HID == HID]:
            newGene = human[human.HID == HID].GeneSymbol.values[0]
        else:
            newGene = ''
    # If it's not human or mouse, check if it's another animal and convert to human
    elif gene in list(homolo.GeneSymbol):
        HID = homolo[homolo.GeneSymbol == gene].HID.values[0]
        if HID in list(human.HID):
            newGene = human[human.HID == HID].GeneSymbol.values[0]
        else:
            newGene = ''
    else:
        newGene = ''
        recordMissingGene(gene)
    return newGene


minNumberInteractors = 5

# ----------- ARCHS4 -----------
with open("X2K_Databases/TF/ARCHS4/Processing/ARCHS4-08-2017_Human_transcription-factor-gene-set.gmt") as human, \
        open("X2K_Databases/TF/ARCHS4/Processing/ARCHS4-02-2018_Mouse_transcription-factor-gene-set.gmt") as mouse, \
        open("X2K_Databases/TF/ARCHS4/ARCHS4-02-2018_Mouse-Human_TF.gmt", "w") as combined:
    for line in human:
        lineSp = line.split("\t")
        geneName = lineSp[0]
        newName =  "ARCHS4" + "_" + "Human"
        combined.write(geneName + "\t"+newName+"\t" + "\t".join(lineSp[2:]))
    for line in mouse:
        lineSp = line.split("\t")
        geneName = lineSp[0]
        newName =  "ARCHS4" + "_" + "Mouse"
        if len(lineSp[2:]) >= minNumberInteractors:
            combined.write(geneName + "\t"+newName+"\t" + "\t".join(lineSp[2:]))

# ----------- ChEA 2016 -----------
with open("X2K_Databases/TF/ChEA_2016/Processing/ChEA_2016.txt") as CHEA, \
        open("X2K_Databases/TF/ChEA_2016/CHEA-2016_Human_TF.gmt", "w") as humanFile,\
        open("X2K_Databases/TF/ChEA_2016/CHEA-2016_Mouse_TF.gmt", "w") as mouseFile:
    for line in CHEA:
        lineSp = line.split("\t")
        gene = lineSp[0].split("_")[0]
        species = lineSp[0].split("_")[-1]
        if species in ['Neurons', 'Ovary', 'Hela', 'Gbm']:
            species = 'Human'
        if species != "Rat":
            info = "-".join(lineSp[0].split("_")[1:-1])
            newName = '_'.join([info, species])
            cleanGenes = [gene.replace(',1.0', '') for gene in lineSp[2:]]
            if len(cleanGenes) >= minNumberInteractors:
                if species == "Human":
                    humanFile.write(gene + "\t" + newName + "\t" + "\t".join(cleanGenes))
                if species == "Mouse":
                    mouseFile.write(gene + "\t" + newName + "\t" + "\t".join(cleanGenes))

# ----------- ENCODE 2015 -----------
with open("X2K_Databases/TF/ENCODE_2015/Processing/ENCODE_TF_ChIP-seq_2015.txt") as ENCODE_2015, \
        open("X2K_Databases/TF/ENCODE_2015/ENCODE-2015_Human_TF.gmt", "w") as humanFile,\
        open("X2K_Databases/TF/ENCODE_2015/ENCODE-2015_Mouse_TF.gmt", "w") as mouseFile:
    for line in ENCODE_2015:
        lineSp = line.split("\t")
        gene = lineSp[0].split("_")[0]
        gene = gene.replace("phosphoS5","")
        gene = gene.replace("phosphoS2", "")
        info = "-".join(lineSp[0].split("_")[1:-1])
        genome = lineSp[0].split("_")[-1]
        speciesDict = {"hg19": "Human", "mm9": "Mouse"}
        species = speciesDict[genome]
        newName = "_".join([info, species])
        if len(lineSp[2:]) >= minNumberInteractors:
            if species=="Human":
                humanFile.write(gene + "\t" + newName + "\t" + "\t".join(lineSp[2:]))
            if species=="Mouse":
                mouseFile.write(gene + "\t" + newName + "\t" + "\t".join(lineSp[2:]))

# ----------- ENCODE 2017 ----------- # DO NOT USE!!!! Moshe made this version, but there may be problems with it. Ali is currently working on an improved version.
with open("X2K_Databases/TF/ENCODE_2017/Processing/encode_transcription_factors_gene_set_2017_08.gmt") as ENCODE_2017, \
        open("X2K_Databases/TF/ENCODE_2017/ENCODE-08-2017_Human?_TF.gmt", "w") as newFile:
    for line in ENCODE_2017:
        lineSp = line.split("\t")
        gene = lineSp[0].split("_")[0]
        info = "-".join(lineSp[0].split("_")[1:-1])
        species = "Human"
        newName = "_".join([gene, info, species]) + "_"
        if len(lineSp[2:]) >= minNumberInteractors:
            newFile.write(newName + "\t\t" + "\t".join(lineSp[2:]))

# ----------- Enrichr -----------
with open("X2K_Databases/TF/Enrichr/Processing/Enrichr_Submissions_TF-Gene_Coocurrence.txt") as Enrichr, \
        open("X2K_Databases/TF/Enrichr/Enrichr-Submissions-TF-Gene-Cooccurence-2018_UnknownSpecies_TF.gmt", "w") as newFile:
    for line in Enrichr:
        lineSp = line.split("\t")
        gene = lineSp[0].split("_")[0]
        species = "UnknownSpecies"
        newName = "_".join(["Enrichr-Submissions", species])
        cleanGenes = [gene.replace(',1.0', '') for gene in lineSp[2:]]
        if len(cleanGenes) >= minNumberInteractors:
            newFile.write(gene + "\t"+ newName +"\t" + "\t".join(cleanGenes))

# ----------- huMAP -----------
with open("X2K_Databases/TF/huMAP/Processing/huMAP_gene_set_2017_07.gmt") as huMAP, \
        open("X2K_Databases/TF/huMAP/huMAP_07-2017_Human_TF.gmt", "w") as newFile:
    for line in huMAP:
        lineSp = line.split("\t")
        gene = lineSp[0]
        species = "Human"
        newName = "huMAP-TF_" + species
        if len(lineSp[2:]) >= minNumberInteractors and gene in TFs:
            newFile.write(gene + "\t" + newName + "\t" + "\t".join(lineSp[2:]))

# ----------- TF PPIs -----------
with open("X2K_Databases/TF/TF_PPIs/Processing/Transcription_Factor_PPIs.txt") as TF_PPIs, \
        open("X2K_Databases/TF/TF_PPIs/TF-PPIs-GENES2FANS-2012_UnknownSpecies_TF.gmt", "w") as newFile:
    for line in TF_PPIs:
        lineSp = line.split("\t")
        gene = lineSp[0]
        species = "UnknownSpecies"
        newName = "TF-PPIs-GENES2FANS_" + species
        if len(lineSp[2:]) >= minNumberInteractors:
            newFile.write(gene + "\t" + newName + "\t" + "\t".join(lineSp[2:]))

# ----------- TF-LOF-Expression-from-GEO -----------
with open("X2K_Databases/TF/TF-LOF_Expression_GEO/Processing/TF-LOF_Expression_from_GEO.txt") as TF_PPIs, \
        open("X2K_Databases/TF/TF-LOF_Expression_GEO/TF-LOF-Expression-GEO_07-2017_Human_TF.gmt", "w") as humanFile,\
        open("X2K_Databases/TF/TF-LOF_Expression_GEO/TF-LOF-Expression-GEO_07-2017_Mouse_TF.gmt", "w") as mouseFile:
    for line in TF_PPIs:
        lineSp = line.split("\t")
        nameSp = lineSp[0].split("_")
        gene = nameSp[0].upper()
        up_dn = nameSp[-1]
        if "human" in nameSp:
            species = 'Human'
            nameSp.remove("human")
        elif "mouse" in nameSp:
            species = 'Mouse'
            nameSp.remove("mouse")
        info = "-".join(nameSp[:-1])
        cleanGenes = [gene.split(",")[0].upper() for gene in lineSp[2:]]
        newName = "_".join([info, species])
        if len(cleanGenes) >= minNumberInteractors:
            if species=="Human":
                humanFile.write(gene + "\t" + newName + "\t" + "\t".join(cleanGenes))
            if species=="Mouse":
                mouseFile.write(gene + "\t" + newName + "\t" + "\t".join(cleanGenes))


# -----------CREEDS Manual TF Perturbations-----------
# Get just the TF experiments
TFs = pd.read_csv("X2K_Databases/General_Resources/compiled-TFs_Mouse-Human.csv", header=None, index_col=False)
TFs = list(TFs[0])

def reformatAndConvert(line, up_dn):
    lineSp = line.split("\t")
    name = lineSp[0].split(" ")
    gene = name[0]
    if "human" in name:
        species = "Human"
    elif "mouse" in name:
        species = 'Mouse'
    else:
        species = "UnknownSpecies"
    newName = "-".join(name[1:] + ["CREEDS-ManualSingleGenePerturbations"]) + "_" + species + "_" + up_dn
    ## Get rid of '1.0' ib DEGs
    DEGs = lineSp[2:][:-1]
    for g in DEGs:
        DEGs[DEGs.index(g)] = g.strip(",1.0").upper()
    newLine = gene.upper() + "\t" + newName + "\t" + "\t".join(DEGs) + '\n'
    return newLine, species


with open("X2K_Databases/General_Resources/CREEDS/Single_Gene_Perturbations_from_GEO_up.txt") as up, \
        open("X2K_Databases/General_Resources/CREEDS/Single_Gene_Perturbations_from_GEO_down.txt") as dn, \
        open("X2K_Databases/TF/CREEDS/CREEDS-manual-single-gene-perturbations-from-GEO_Human_TF.txt", "w") as humanGMT_TF,\
        open("X2K_Databases/TF/CREEDS/CREEDS-manual-single-gene-perturbations-from-GEO_Mouse_TF.txt", "w") as mouseGMT_TF:
    # Run through each line of both up and dn files
    dnLines = dn.readlines()
    for i, line in enumerate(up):
        gene = line.split(" ")[0]
        # Include only the TF experiments
        if gene in TFs or gene.upper() in TFs:
            # Format Up line
            newLineUp, species = reformatAndConvert(line, 'up')
            # Format Dn line
            newLineDn, species = reformatAndConvert(dnLines[i], 'dn')
            if species=="Human":
                humanGMT_TF.write(newLineUp)
                humanGMT_TF.write(newLineDn)
            if species=="Mouse":
                mouseGMT_TF.write(newLineUp)
                mouseGMT_TF.write(newLineDn)

# -----------CREEDS Manual Kinase Perturbations -----------
KINASES = pd.read_csv("X2K_Databases/General_Resources/compiled-Kinases_Mouse-Human.csv", header=None, index_col=False)
KINASES = list(KINASES[0])

with open("X2K_Databases/General_Resources/CREEDS/Single_Gene_Perturbations_from_GEO_up.txt") as up, \
        open("X2K_Databases/General_Resources/CREEDS/Single_Gene_Perturbations_from_GEO_down.txt") as dn, \
        open("X2K_Databases/KINASE/CREEDS/Processing/CREEDS-manual-single-gene-perturbations-from-GEO_Human_KINASE.txt", "w") as humanGMT_KIN,\
        open("X2K_Databases/KINASE/CREEDS/Processing/CREEDS-manual-single-gene-perturbations-from-GEO_Mouse_KINASE.txt", "w") as mouseGMT_KIN,\
        open("X2K_Databases/KINASE/CREEDS/CREEDS-manual-single-gene-perturbations-from-GEO_Mouse-Human_KINASE.txt", "w") as GMT_KIN:
    # Run through each line of both up and dn files
    dnLines = dn.readlines()
    for i, line in enumerate(up):
        gene = line.split(" ")[0]
        # Include only the TF experiments
        if gene in KINASES or gene.upper() in KINASES:
            # Format Up line
            newLineUp, species = reformatAndConvert(line, 'up')
            # Format Dn line
            newLineDn, species = reformatAndConvert(dnLines[i], 'dn')
            if species == "Human":
                humanGMT_KIN.write(newLineUp)
                humanGMT_KIN.write(newLineDn)
                # Combined file
                GMT_KIN.write(newLineUp)
                GMT_KIN.write(newLineDn)
            if species == "Mouse":
                mouseGMT_KIN.write(newLineUp)
                mouseGMT_KIN.write(newLineDn)
                # Combined file
                GMT_KIN.write(newLineUp)
                GMT_KIN.write(newLineDn)




# -----------[2] BIOGRID -----------
# XXXXXXXXX USE MOSHE'S  JUPYTER PIPELINE INSTEAD XXXXXXXXX

# GMT format
# Reformat and Separate high- vs. low-throughput
BIOGRID = pd.read_table("X2K_Databases/General_Resources/BioGRID/BIOGRID-MV-Physical-3.4.154.tab2.txt", header=0)
BIOGRID.head()
BIOGRID.shape
BIOGRID["Official Symbol Interactor A"].nunique()

with open("X2K_Databases/General_Resources/BioGRID/BIOGRID_LowThroughput.txt", "w") as BIOGRID_Low, \
        open("X2K_Databases/General_Resources/BioGRID/BIOGRID_HighThroughput.txt", "w") as BIOGRID_High:
    interactors = list(BIOGRID["Official Symbol Interactor A"].unique())


    def getSortedLine(inter, subset="All"):
        interSub = BIOGRID[BIOGRID["Official Symbol Interactor A"] == inter]
        # Get just subset
        if subset == "Low":
            interSub_thru = interSub[interSub.Throughput == 'Low Throughput']
        elif subset == "High":
            interSub_thru = interSub[interSub.Throughput == 'High Throughput']
        else:
            interSub_thru = interSub
        Blist = interSub_thru[["Official Symbol Interactor B", "Score"]]
        try:
            # Get mean score per B gene
            pd.to_numeric(Blist.Score, errors='coerce')
            Blist_avg = Blist.groupby("Official Symbol Interactor B").mean()
            # Sort df by avg Score
            Blist_avg_sort = Blist_avg.sort_values(['Score'], ascending=[0])
            geneList = list(Blist_avg_sort.index)
        except:  # Avoid errors when there's no scores and just use unsorted list
            geneList = list(Blist["Official Symbol Interactor B"])
        newLine = inter + "\t\t" + "\t".join(geneList) + "\n"
        return newLine


    for inter in interactors:
        # Low Throughout
        BIOGRID_Low.write(getSortedLine(inter, "Low"))
        # High Throughput
        BIOGRID_High.write(getSortedLine(inter, "High"))
# Get low-thru TFs (GMT)
with open("X2K_Databases/General_Resources/BioGRID/BIOGRID_LowThroughput.txt") as BIOGRID_Low, \
        open("X2K_Databases/TF/BioGRID/BIOGRID_LowThroughput_TF.txt", "w") as BIOGRID_Low_TF:
    for line in BIOGRID_Low:
        lineSP = line.strip("\n")
        gene = lineSP.split("\t")[0]
        newGene = geneConverter(gene)
        DEGs = lineSP.split('\t')[2:]

        if newGene in TFs and len(DEGs) > 0 and gene != newGene:
            DEGs_conv = []
            for g in DEGs:
                DEGs_conv.append(geneConverter(g))

# Get high-thru TFs (GMT)
with open("X2K_Databases/General_Resources/BioGRID/BIOGRID_HighThroughput.txt") as BIOGRID_High, \
        open("X2K_Databases/TF/BioGRID/BIOGRID_HighThroughput.txt") as BIOGRID_High_TF:
    BIOGRID_Low

# -----------Split JASPAR-TRANSFAC TF -----------
# Create separate background files
JASP_TRANS = pd.read_csv("X2K_Databases/TF/JASPAR-TRANSFAC/JASPAR-TRANSFAC_background.csv", header=None)  # "transfac_background" actually has both transfac AND Jaspar
JASP = JASP_TRANS[JASP_TRANS[5] == "jaspar"]
TRANS = JASP_TRANS[JASP_TRANS[5].isin(["transfac", "transfacs"])]


def makeJaspTransBackground(df, species, databaseName):
    if databaseName=="JASPAR-TRANSFAC":
        subset = df[df[7] == species.lower()]
        subset.to_csv("X2K_Databases/TF/JASPAR-TRANSFAC/" + databaseName + "_" + species + "_background.csv",\
                      header=["", "", "GeneA", "GeneB", "", "Database", "Filler", "Species", "DateTime"], index=False, index_label=False)
    else:
        sub = df[df[5].isin(databaseName)]
        subset = sub[sub[7]==species.lower()]
        subset.to_csv("X2K_Databases/TF/JASPAR-TRANSFAC/"+databaseName[0].upper()+"_"+species+"_background.csv",\
                      header=["","","GeneA","GeneB","","Database", "Filler", "Species", "DateTime"], index=False, index_label=False)
for s in ["Human", "Mouse"]:
    #makeJaspTransBackground(JASP_TRANS, s, databaseName=["jaspar"])
    #makeJaspTransBackground(JASP_TRANS, s, databaseName=["transfac","transfacs"])
    makeJaspTransBackground(JASP_TRANS, s, databaseName="JASPAR-TRANSFAC")



# Create separate GMTs
def backgroundToGMT(database, speciesList=["Mouse","Human"]):
    for species in speciesList:
        background = pd.read_csv("X2K_Databases/TF/JASPAR-TRANSFAC" + "/" + database.upper()+"_"+species+"_background.csv", names=["","","GeneA","GeneB","","Database", "Filler", "Species", "DateTime"]) # names=["","","","GeneA","GeneB","GeneC","", "Filler", "Species", "Date", "Time"]
        if species == "mouse":
            background  =  background.iloc[1:]
        subset = background[background.Species==species.lower()]
        tfList = subset["GeneA"].unique()
        with open("X2K_Databases/TF/" + database + "/"+ database.upper()+"_"+species+".gmt","w") as newFile:
            for tf in tfList:
                tf_targets = [x.upper() for x in list(subset[subset["GeneA"] == tf]["GeneB"])]
                newName = database + "_" + species
                if len(tf_targets) > 0:
                    newFile.write(tf.upper() + "\t" + newName + "\t" + "\t".join(tf_targets) + "\n")
#backgroundToGMT(database="JASPAR", speciesList=["Mouse","Human"])
#backgroundToGMT(database="TRANSFAC", speciesList=["Mouse","Human"])
backgroundToGMT(database="jaspar-transfac", speciesList=["Mouse","Human"])


# -----------iREF-----------

def mitabGeneNames(mitab_path, output_path=False):
    import pandas as pd
    iref = pd.read_table(mitab_path)
    iref["hgnc_A"] = ""
    iref["hgnc_B"] = ""
    for i, row in iref.iterrows():
        # Print progress every 100th iteration
        if (i + 1) % 100 == 0:
            print(str(round((i + 1) / iref.shape[0] * 100, 2)) + "% complete")
        # Extract geneNames
        for item in row["aliasA"].split("|"):
            if item.startswith("hgnc"):
                geneName = item.split(":")[1]
                iref.at[row.name, "hgnc_A"] = geneName.strip().upper()
        for item in row["aliasB"].split("|"):
            if item.startswith("hgnc"):
                geneName = item.split(":")[1]
                iref.at[row.name, "hgnc_B"] = geneName.strip().upper()
        iref_filt = iref[(iref["hgnc_A"] != '') & (iref["hgnc_B"] != '')]
    # Write file
    if output_path != False:
        iref_filt.to_csv(output_path, sep="\t")
    return iref_filt

def mitabGeneNames(mitab_path, output_path=False):
    import pandas as pd
    iref = pd.read_table(mitab_path)
    iref = iref[(iref.aliasA.str.contains("hgnc:"))&(iref.aliasB.str.contains("hgnc:"))]
    # Strip down to get geneSymbols
    iref["hgnc_A"] = iref.aliasA.str.split("|").apply(lambda parts: [x for x in parts if x.startswith("hgnc:") ][0]).str.replace("hgnc:","").str.upper()
    iref["hgnc_B"] = iref.aliasB.str.split("|").apply(lambda parts: [x for x in parts if x.startswith("hgnc:") ][0]).str.replace("hgnc:","").str.upper()
    # Filter NAs and blanks
    iref = iref.dropna(subset=['hgnc_A','hgnc_B'])
    iref_filt = iref[(iref["hgnc_A"] != '') & (iref["hgnc_B"] != '')]
    # Write file
    if output_path != False:
        iref_filt.to_csv(output_path, sep="\t")
    return iref_filt

mouse_df = mitabGeneNames("X2K_Databases/PPI/iREF/Processing/10090.mitab.01-22-2018.txt","X2K_Databases/PPI/iREF/Processing/mouse_mittab.txt")
human_df = mitabGeneNames("X2K_Databases/PPI/iREF/Processing/9606.mitab.01-22-2018.txt","X2K_Databases/PPI/iREF/Processing/human_mittab.txt")

# mouse_df = pd.read_table("X2K_Databases/PPI/iREF/Processing/mouse_mittab.txt")
# human_df = pd.read_table("X2K_Databases/PPI/iREF/Processing/human_mittab.txt")
iREF_df = pd.concat([mouse_df, human_df])


def mitabToSig(iREF_df, output_path):
    # Remove any interactions without PMIDS
    iREF_df = iREF_df[iREF_df.pmids!="-"]
    iREF_df['Interaction'] = iREF_df.hgnc_A + "@" + iREF_df.hgnc_B
    allInteractions = list( iREF_df['Interaction'].unique())
    for i,inter in enumerate(allInteractions):
        # Progress bar
        print(str(i) + " : " + str(round((i + 1) / iREF_df.shape[0] * 100, 2)) + "% complete")
        geneA=[]; geneB=[]; PMIDS=[]
        inter = str(inter)
        reverseInter = inter.split("@")[1]+"@"+inter.split("@")[0]
        subset = iREF_df[(iREF_df.Interaction==inter)|(iREF_df.Interaction==reverseInter)]
        pmidList = subset.pmids.str.replace("pubmed:","").unique().tolist()
        if len(pmidList)<=15:
            geneA.append( inter.split("@")[0] )
            geneB.append( inter.split("@")[1] )
            PMIDS.append(",".join(pmidList))
    # Make dataframe
    pd.DataFrame(np.column_stack([geneA, geneB, ]))
    with open(output_path, "w") as sig:
        # Only include interactions with no more than 5 PubMedIDs (to exclude any MassSpec experiments)
            newLine = geneA + " NA NA NA NA " + geneB + " NA NA NA NA NA " + PMIDS + '\n'
            sig.write(newLine)


#mitabToSig(mouse_df, "X2K_Databases/PPI/iREF/Processing/iREF_02-2018_mouse_PPI.sig")
#mitabToSig(human_df, "X2K_Databases/PPI/iREF/Processing/iREF_02-2018_human_PPI.sig")
mitabToSig(iREF_df, output_path="X2K_Databases/PPI/iREF/Processing/iREF_02-2018_Mouse-Human_PPI.sig")


# iREF GMT
def mitabToGMT(iREF_df, output_path, subsetList="", speciesList=["Human", "Mouse"], mergeSpecies=False):
    with open(output_path, "w") as gmt:
        # Include only TFs
        if subsetList != "":
            df = iREF_df[iREF_df.hgnc_A.isin(subsetList)]
        else:
            df = iREF_df
        for i, gene in enumerate(df.hgnc_A.unique()):
            print("Processing: "+str(gene))
            # Keep species as separate lines
            def writeSpeciesLine(species):
                speciesDict = {"Mouse": "taxid:10090(Mus musculus)", "Human": "taxid:9606(Homo sapiens)"}
                sub = df[(df.hgnc_A == gene) & (df.taxa == speciesDict[species])]
                Bgenes = list(sub.hgnc_B.unique())
                if len(Bgenes) >= 5:  # Ony include genes with at least 5 interactors
                    newLine = gene + "\t" + "iREF_" + species + "\t" + "\t".join(Bgenes) + "\n"
                    gmt.write(newLine)
            if mergeSpecies == False:
                for spec in speciesList:
                    writeSpeciesLine(spec)
            # Merge each gene into one line across all species
            elif mergeSpecies == True:
                sub = df[df.hgnc_A == gene]
                Bgenes = list(sub.hgnc_B.unique())
                Bgenes = [str(x) for x in Bgenes]
                species = "Mouse-Human"
                Group, Family = getKinaseGroupAndFamily([str(gene)])
                if len(Bgenes) >= 5:  # Ony include genes with at least 5 interactors
                    newLine = str(gene) + "\t" + "-".join(["KinaseGroup:"+Group,"KinaseFamily:"+Family, "iREF_"]) + species + "\t" + "\t".join(Bgenes) + "\n"
                    gmt.write(newLine)


# iREF GMT: Species as separate files
## TFs
TFs = list(pd.read_csv("X2K_Databases/General_Resources/compiled-TFs_Mouse-Human.csv", header=None)[0])
mitabToGMT(mouse_df, "X2K_Databases/TF/iREF/iREF_02-2018_Mouse_TF.gmt", subsetList=TFs, speciesList=["Mouse"])
mitabToGMT(human_df, "X2K_Databases/TF/iREF/iREF_02-2018_Human_TF.gmt", subsetList=TFs, speciesList=["Human"])

# iREF GMT: Species as separate lines
## KINASES
mitabToGMT(iREF_df, "X2K_Databases/KINASE/iREF/iREF_02-2018_Mouse-Human_KINASE.gmt", subsetList=KINASES, speciesList=["Human", "Mouse"])
## All Genes  (takes a long time)
mitabToGMT(iREF_df, "X2K_Databases/PPI/iREF/Processing/iREF_02-2018_Mouse-Human_AllGenes.gmt", subsetList="", speciesList=["Human", "Mouse"])

# iREF GMT: Species merged into same line
## TFs
mitabToGMT(iREF_df, "X2K_Databases/TF/iREF/Processing/iREF_02-2018_Mouse-Human_AllGenes-Merged.gmt", subsetList=TFs, mergeSpecies=True)
## KINASES
mitabToGMT(iREF_df, "X2K_Databases/KINASE/iREF/Processing/iREF_02-2018_Mouse-Human_AllGenes-Merged.gmt", subsetList=KINASES, mergeSpecies=True)
# ## All Genes  (takes a long time)
mitabToGMT(iREF_df, "X2K_Databases/PPI/iREF/Processing/iREF_02-2018_Mouse-Human_AllGenes-Merged.gmt", subsetList="", mergeSpecies=True)






# -----------iPTMnet-----------
with open("X2K_Databases/KINASE/iPTMnet/Processing/readme.txt") as file:
    readme = file.readlines()
# Data
score = pd.read_table("X2K_Databases/KINASE/iPTMnet/Processing/score.txt", header=None)
score.columns = ["substrate_AC", "site", "enzyme_AC", "ptm_type", "score"]

score = score.dropna(subset=["enzyme_AC"])
score.groupby("enzyme_AC").mean()
score.head()


# Protein metadata
cols = ["ptm_type", "source", "substrate_AC", "substrate_genename", "organism", "site", "enzyme_AC","enzyme_genename", "note", "pmid"]
ptm = pd.read_table("X2K_Databases/KINASE/iPTMnet/Processing/ptm.txt", header=None, names=cols, dtype=dict(zip(cols, [str]*len(cols))))
ptm.head()
# Protein metadata 2?
protein = pd.read_table("X2K_Databases/KINASE/iPTMnet/Processing/protein.txt", header=None)
protein.columns = ["UniProtAC", "UniProtID", "protein_name", "genename", "organism", "PRO_id", "SwissProt/TrEmbl"]
protein.head()

# Convert substrates
substrateDict = dict(zip(ptm.substrate_AC, zip(ptm.substrate_genename, ptm.pmid, ptm.organism)))
enzymeDict = dict(zip(ptm.enzyme_AC, zip(ptm.enzyme_genename, ptm.organism)))

# Append synonyms as entries in kgroupDict


enzyme_Genes = [str(enzymeDict[x][0]) for x in score.enzyme_AC.tolist() if x in enzymeDict.keys()]
kgroupDict = addSynonymsToKinaseGroups(kgroupDict, enzyme_Genes)

# iterate over all interactions in iPTMnet
substrateGenes=[]; enzymeGenes=[]; Scores=[]; PMIDS = []; Species=[]; KinaseGroup=[]; KinaseFamily=[]
for i, row in score.iterrows():
    if row.substrate_AC in substrateDict.keys() and row.enzyme_AC in enzymeDict.keys():
        substrateGenes.append(substrateDict[row.substrate_AC][0])
        enzymeGenes.append(enzymeDict[row.enzyme_AC][0])
        if enzymeDict[row.enzyme_AC][0] in kgroupDict.keys():
            KinaseGroup.append( kgroupDict[enzymeDict[row.enzyme_AC][0]][0] )
            KinaseFamily.append( kgroupDict[enzymeDict[row.enzyme_AC][0]][1] )
        else:
            KinaseGroup.append("NA")
            KinaseFamily.append("NA")
        Scores.append(row.score)
        PMIDS.append(substrateDict[row.substrate_AC][1])
        Species.append(enzymeDict[row.enzyme_AC][1])
    else:
        print(">Can't find protein<")

# DataFrame
import numpy as np
convertedScores = pd.DataFrame(np.column_stack([substrateGenes, enzymeGenes, KinaseGroup, KinaseFamily, Species, Scores, PMIDS]), \
                               columns=["substrateGenes", "enzymeGenes", "KinaseGroup", "KinaseFamily","Species", "Scores", "PMIDS"])
convertedScores = convertedScores.apply(lambda x: x.astype('category'))
convertedScores['Scores'] = pd.to_numeric(Scores, errors='ignore')*1.00


# Get the mean for each enzyme
keepCols = list(convertedScores.columns); keepCols.remove("Scores")
groupedScores = convertedScores.groupby(keepCols)['Scores'].mean().reset_index()
groupedScores.to_csv("X2K_Databases/KINASE/iPTMnet/Processing/score_PPI.txt", sep="\t")

groupedScores = pd.read_table("X2K_Databases/KINASE/iPTMnet/Processing/score_PPI.txt")
groupedScores['substrateGenes'] =  groupedScores['substrateGenes'].astype(str)
groupedScores['enzymeGenes'] =  groupedScores['enzymeGenes'].astype(str)

# Prepare PPI in SIG format
## Subset humans and mice
mhScores = groupedScores[(groupedScores.Species == "Homo sapiens (Human)") | (groupedScores.Species == "Mus musculus (Mouse)")]
## Convert to simple Species names
mhScores['Species'] = mhScores['Species'].map({"Homo sapiens (Human)":"Human","Mus musculus (Mouse)":"Mouse"})
# Make genes all uppecase
mhScores['substrateGenes'] = mhScores['substrateGenes'].str.upper()
mhScores['enzymeGenes'] = mhScores['enzymeGenes'].str.upper()
mhScores.shape
# Remove entries without PubMedIDs
mhScores = mhScores.dropna(subset=["PMIDS","enzymeGenes","substrateGenes"])
mhScores.shape
# Remove any duplicates
mhScores = mhScores.drop_duplicates()
mhScores.shape


# Write SIG for all genes in iPTMnet
NAs = ["NA"]*len(mhScores.enzymeGenes)
# NOTE: reordered columns so enzyme (e.g. kinase) is in the first col
iPTMnet_SIG = pd.DataFrame(np.column_stack([mhScores.enzymeGenes,NAs, NAs, NAs, NAs,mhScores.substrateGenes, \
                                            NAs, NAs, mhScores.KinaseGroup, mhScores.KinaseFamily, \
                                            mhScores.Species, round(mhScores.Scores,3), mhScores.PMIDS]) )
iPTMnet_SIG.to_csv("X2K_Databases/KINASE/iPTMnet/Processing/iPTMnet-2018_Mouse-Human_AllGenes_PPI.sig", sep=" ", header=False, index=False)

# Make GMT for KINASES in iPTMnet
kinaseScores = iPTMnet_SIG[iPTMnet_SIG[0].isin(KINASES)]
kinaseScores.columns = ['enzymeGenes','NA','NA','NA','NA','substrateGenes','NA','NA', 'KinaseGroup', 'KinaseFamily', 'Species', 'Scores', 'PMIDS']

with open("X2K_Databases/KINASE/iPTMnet/iPTMnet-2018_Mouse-Human_KINASES.gmt","w") as GMT:
    uniqueKinases = kinaseScores['enzymeGenes'].unique().tolist()
    len(uniqueKinases)
    for ki in uniqueKinases:
        kiSub = kinaseScores[kinaseScores['enzymeGenes'] == ki]
        for species in list(kiSub['Species'].unique()):
            spSub = kiSub[kiSub['Species']==species]
            interactors = spSub.sort_values(by=['Scores'], ascending=False)['substrateGenes'].values
            newName = "iPTMnet-"+"KinaseGroup:"+str(kiSub['KinaseGroup'].unique()[0])+"-KinaseFamily:"+str(kiSub['KinaseFamily'].unique()[0])+"_"+str(species)
            if len(interactors)>=1:
                GMT.write(ki +"\t"+ newName+"\t"+"\t".join(interactors)+"\n")


# -----------PHOSPHOSITE_PLUS-----------
phos = pd.read_table("X2K_Databases/KINASE/PhosphositePlus/Processing/PhosphositePlus.Kinase_Substrate_Dataset.02-2018.txt", skiprows=[0,1,2])
phos.head()
noRXN = phos[(phos.IN_VIVO_RXN!="X") & (phos.IN_VITRO_RXN!="X")]
noRXN.shape[0] # Just double checking that all kinase-substrate pairings have at least one kind of interaction

## Convert to KINASE GMT
with open("X2K_Databases/KINASE/PhosphositePlus/PhosphositePlus-02-2018_Mouse-Human_KINASES.gmt","w") as PHOS:
    phosHM = phos[phos.KIN_ORGANISM.isin(["human","mouse"]) & phos.SUB_ORGANISM.isin(["human","mouse"])]
    for k in list(phosHM.GENE.unique()):
        print("Processing: "+k)
        phosSub = phosHM[phosHM.GENE==k]
        Group, Family = getKinaseGroupAndFamily([k])
        for species in list(phosSub.KIN_ORGANISM.unique()):
            Species = species[0].upper()+species[1:]
            phosSpec = phosSub[(phosSub.KIN_ORGANISM == species) & (phosSub.SUB_ORGANISM == species)]
            substrates = list(set(phosSub.SUB_GENE.values))
            cleanSubstrates = [x for x in substrates if str(x) != 'nan']
            cleanSubstrates = [x.upper() for x in cleanSubstrates]
            newName = "-".join(["KinaseGroup:"+Group,"KinaseFamily:"+Family])+"-PHOSPHOSITE_"+Species
            if len(substrates)>=1:
                PHOS.write(k.upper()+"\t"+newName+"\t"+"\t".join(cleanSubstrates)+"\n")


# -----------NetworkIN-----------
## Kinase-substrate predictions from neural networks
## "id" is the kinase name
net = pd.read_table("X2K_Databases/KINASE/NetworkIN_2017/Processing/networkin-05-2017_human_predictions_3.1.tsv")
net.head()

## Make GMT
### Need to apply a cutoff for 'networkin_score', otherwise it lists too many substrates
with open("X2K_Databases/KINASE/NetworkIN_2017/NetworkIN-05-2017_UnknownSpecies_KINASES.gmt","w") as newNet:
    for k in list(net.id.unique()):
        print("Processing: "+k)
        # Convert annoying greek alphabet to actual GeneSymbol
        for substring in ['alpha','beta','gamma','delta','epsilon','iota','theta','zeta']:
            if substring in k:
                newK = k.replace(substring,substring[0]).upper()
            else:
                newK = k.upper()
        netSub = net[net.id==k]
        # Alternatively could use "netphorest_score"
        netMean = netSub.groupby('substrate_name').mean().reset_index()
        substrates = list( netMean.sort_values(by=['networkin_score'], ascending=False).substrate_name )
        substrates = [x.upper() for x in substrates]
        group,family = getKinaseGroupAndFamily([newK])
        newName = "KinaseGroup:"+group+"-KinaseFamily:"+family+"-NetworkIN_UnknownSpecies_"
        newNet.write(newK+"\t"+newName+"\t"+"\t".join(substrates)+"\n")


# -----------HMS LINCS KINOMEscan-----------
kscan = pd.read_excel("X2K_Databases/General_Resources/KINOMEscan/Processing/HMS-LINCS_KinomeScan_Datasets_2017-10-23.xlsx")
kscanProts = pd.read_excel("X2K_Databases/General_Resources/KINOMEscan/Processing/proteins_20180116225538.xlsx")
#kscanProts = kscanProts.dropna(subset=['Gene Symbol'])
kscanProts.head()
geneDict = dict(zip(kscanProts['Name'], kscanProts['Gene Symbol']))

# [1] Collect all data first
superData = pd.DataFrame(columns=['controlType', 'datapointName', 'datapointUnit', 'datapointValue',\
   'datarecordID', 'hmsDatasetID', 'protein_ppCenterBatchID',\
   'protein_ppCenterCanonicalID', 'protein_ppName', 'protein_pplincsid',\
   'recordedPlate', 'recordedWell', 'smallmolecule_smCenterBatchID',\
   'smallmolecule_smCenterCompoundID', 'smallmolecule_smLincsID',\
   'smallmolecule_smName', 'smallmolecule_smSalt', 'GeneSymbol'] )
for id in kscan.dataset_id:
    print('Adding Data to SuperData: '+str(id))
    # Download data from json
    url = "http://lincs.hms.harvard.edu/db/api/v1/"+"datasetdata/"+str(id)
    data = pd.read_json(url)
    if len(data)>0:
        superData = superData.append(data)
    else:
        print('Missing Data: ' + str(id))

# [2] Convert superData into GMT (with each drug combined into one line)
with open("X2K_Databases/General_Resources/KINOMEscan/KINOMEscan-01-2018_Small-Molecule-Kinase-Perturbations.gmt","w") as newFile:
    for drug in list(superData.smallmolecule_smName.unique()):
        data = superData[superData.smallmolecule_smName==drug]
        data = data[data.datapointName.isin(["percentControl","dissociationConstant","percentInhibition"])]
        # Convert genes to GeneSymbols
        data['GeneSymbol'] = [geneDict[x] for x in data.protein_ppName]
        # Get mean value of each gene
        data['datapointValue'] = pd.to_numeric(data['datapointValue'], errors='ignore') * 1.00
        dataGrouped = data.groupby('GeneSymbol').mean().reset_index()
        dataGrouped = dataGrouped.dropna(subset=['GeneSymbol'])
        dataGrouped['GeneSymbol'] = [x.upper() for x in dataGrouped['GeneSymbol']]
        # Order genes by datapointValue (direction of order depends on metric?)
        if str(dataGrouped.datapointValue[0])=="percentControl":
            geneList = list( dataGrouped.sort_values(by=['datapointValue'], ascending=False).GeneSymbol )
        elif str(dataGrouped.datapointValue[0])=="percentInhibition":
            geneList = list( dataGrouped.sort_values(by=['datapointValue'], ascending=False).GeneSymbol)
        else:
            geneList = list(dataGrouped.sort_values(by=['datapointValue'], ascending=False).GeneSymbol)
        # Write to file
        drugInfo = "HMS-LINCS-KINOMEscan_smCenterCompoundID:"+str(data.smallmolecule_smCenterCompoundID.values[0])+"_smLincsID:"+str(data.smallmolecule_smLincsID.values[0])+"_smName:"+str(data.smallmolecule_smName.values[0])
        newFile.write(drugInfo+"\t\t"+"\t".join(geneList)+"\n")



# -----------Split KEA file-----------

KEA = pd.read_csv("X2K_Databases/KINASE/KEA_datasets/Processing/kinase-protein_interactions.csv", header=None, names=['KinaseFamily','KinaseGroup','Kinase','Substrate','PMIDS','Databases'])
KEA.head()
# Get complete list of databases used in
## 'SAVI' == 'SNAVI'?...
allDbs=[]
for dbs in KEA.Databases.unique():
    dbsSp = dbs.split(";")
    for db in dbsSp:
        if db not in allDbs:
            allDbs.append(db)
## Split file into individual csv's
for db in allDbs:
    print("Processing CSV: "+db)
    newDF=pd.DataFrame(columns=['KinaseFamily','KinaseGroup','Kinase','Substrate','PMIDS','Databases'])
    newPath = "X2K_Databases/KINASE/KEA_datasets/Processing/"+db+"_UnknownSpecies_KINASES.csv"
    for i,row in KEA.iterrows():
        if db in row.Databases.split(";"):
            newRow = row
            newRow['Databases'] = db
            newDF = newDF.append(newRow)
    newDF.to_csv(newPath, header=None, index=False)
## Convert each csv to GMT
for db in allDbs:
    print("Processing GMT: "+db)
    path = "X2K_Databases/KINASE/KEA_datasets/Processing/" + db + "_UnknownSpecies_KINASES.csv"
    gmtPath = "X2K_Databases/KINASE/KEA_datasets/"+db+"_UnknownSpecies_KINASES.gmt"
    df = pd.read_csv(path, header=None, names=['KinaseFamily','KinaseGroup','Kinase','Substrate','PMIDS','Databases'])
    with open(gmtPath,"w") as newGMT:
        for k in list(df.Kinase.unique()):
            dfSub = df[df.Kinase==k]
            substrates = list(dfSub.Substrate.unique())
            group = dfSub.KinaseGroup.values[0]
            family = dfSub.KinaseFamily.values[0]
            newName = "-".join(["KinaseGroup:"+group,"KinaseFamily:"+family,db])+"_UnknownSpecies"
            if len(substrates)>=1:
                newGMT.write(k+"\t"+newName+"\t"+"\t".join(substrates)+"\n")

############ Recombine new KEA file ##########
files = os.listdir("X2K_Databases/KINASE/")
files.remove("KEA_datasets")
for db in files:
    Dir = "X2K_Databases/KINASE/"+db+"/"
    for file in os.listdir():
        if file.endswith(".gmt"):
            print("F")






# -----------MINT 2018-----------
import pandas as pd
## colNames info found here: https://psicquic.github.io/MITAB27Format.html
colNames = ["uniprotA", "uniprotB", "InteractionID", "InteractionID2","GeneInfoA","GeneInfoB","InteractionMethods","Authors","PMID","SpeciesA", "SpeciesB", "InteractionType", "SourceDatabase", "InteractionID3", "ConfidenceScore"]
hMINT = pd.read_table("X2K_Databases/PPI/MINT/Processing/MINT-03-2018_Human_raw.txt", names=colNames)

# Show contents example
for i in range(0,14):
    print(str(i) + ":    "+ hMINT.iloc[:,i][0] )

mMINT = pd.read_table("X2K_Databases/PPI/MINT/Processing/MINT-03-2018_Mouse_raw.txt", names=colNames)
#MINT = hMINT.append(mMINT)

def extractMINT(MINTdf, SPECIES="", writeCSV=True):
    newCols = ["UniprotA", "UniprotB", "GeneSymbolA", "GeneSymbolB", "SpeciesA", "SpeciesB", "PMID", "ConfidenceScore","InteractionType", "SPECIES"]
    newDF=pd.DataFrame(columns=newCols)
    import numpy as np

    def getGeneSymbol(raw):
        geneInfo = raw.split("|")
        for item in geneInfo:
            if item.endswith("(gene name)"):
                GeneSymbol = item.split(":")[1].replace("(gene name)", "")
                break
            elif item != "-":
                GeneSymbol = item.split(":")[1].split("(")[0].split("_")[0]
            else:
                GeneSymbol = "NA"
        return GeneSymbol
    def getSpecies(raw):
        if raw=="-":
            species = "NA"
        else:
            species = raw.split("|")[0].split("(")[1].strip(')')
        return species
    def getPMID(raw):
        pmidInfo = raw.split("|")
        for item in pmidInfo:
            if item.startswith("pubmed"):
                pmid = item.replace("pubmed:","")
                break
            else:
                pmid = "NA"
        return pmid
    def getConfidence(raw):
        if raw=="-":
            CS = "NA"
        else:
            CS = raw.split(":")[1]
        return CS
    def getInteraction(raw):
        if raw=="-":
            inter = "NA"
        else:
            inter = "_".join( raw.split("(")[1].strip(")").split(" ") )
        return inter
    for i,row in MINTdf.iterrows():
        print(str(i)+" : "+str(round(i/MINTdf.shape[0]*100,2)),"% complete")
        UniprotA = row.uniprotA.replace('uniprotkb:', '')
        UniprotB = row.uniprotB.replace('uniprotkb:', '')
        GeneSymbolA = getGeneSymbol(row.GeneInfoA)
        GeneSymbolB = getGeneSymbol(row.GeneInfoB)
        SpeciesA = getSpecies(row.SpeciesA)
        SpeciesB = getSpecies(row.SpeciesB)
        PMID = getPMID(row.PMID)
        ConfidenceScore  = getConfidence(row.ConfidenceScore)
        InteractionType = getInteraction(row.InteractionType)
        newDF = newDF.append( pd.DataFrame(np.column_stack([UniprotA, UniprotB, GeneSymbolA, GeneSymbolB, SpeciesA, SpeciesB, PMID, ConfidenceScore, InteractionType, SPECIES]), columns=newCols) )
    if writeCSV==True:
        newDF.to_csv("X2K_Databases/PPI/MINT_2018/MINT_"+SPECIES+"_table.csv")
    return newDF

# mouseDF = extractMINT(mMINT, "Mouse")
mouseDF = pd.read_csv("X2K_Databases/PPI/MINT_2018/Processing/MINT_Mouse_table.csv")
# humanDF = extractMINT(hMINT, "Human")
humanDF = pd.read_csv("X2K_Databases/PPI/MINT_2018/Processing/MINT_Human_table.csv")

MINTdf2 = mouseDF.append(humanDF)
MINTdf2['ConfidenceScore'] = pd.to_numeric(MINTdf2.ConfidenceScore, errors='coerce')*1.00
MINTdf2.head()

def MINTtoKinaseGMT(MINTdf2, species=""):
    with open("X2K_Databases/KINASE/MINT_2018/MINT-03-2018_"+species+"_KINASE.gmt", "w") as GMT:
        MINTkin = MINTdf2[MINTdf2.GeneSymbolA.isin(KINASES)]
        for gene in list(MINTkin.GeneSymbolA.unique()):
            df = MINTkin[MINTkin.GeneSymbolA==gene]
            dfGrouped  = df.groupby('GeneSymbolB').mean()
            interactors = list( dfGrouped.sort_values(by=['ConfidenceScore'], ascending=False).index )
            newName = "MINT-03-2018_Kinases_"+species
            if len(interactors)>=5:
                GMT.write(gene+"\t"+newName+"\t"+"\t".join(interactors)+"\n")
MINTtoKinaseGMT(MINTdf2, "Mouse-Human")

def MINTtoPPIsig(MINTdf2, species=""):
    sigDF = MINTdf2.loc[:,["GeneSymbolA", "GeneSymbolB", "SPECIES", "ConfidenceScore", "InteractionType", "PMID"]]
    sigDF['GeneSymbolA'] = sigDF['GeneSymbolA'].str.upper()
    sigDF['GeneSymbolB'] = sigDF['GeneSymbolB'].str.upper()
    # Insert NA cols
    NAs = pd.Series(["NA"]*MINTdf2.shape[0])
    for ind in range(4):
        sigDF.insert(loc=1, column="NA"+str(ind), value=NAs)
    for ind in range(3):
        sigDF.insert(loc=6, column="NA"+str(ind+4), value=NAs)
    sigDF.to_csv("X2K_Databases/PPI/MINT_2018/MINT-03-2018_" + species + "_PPI.sig", sep=" ", header=False, index=False)
MINTtoPPIsig(MINTdf2, species='Mouse-Human')

