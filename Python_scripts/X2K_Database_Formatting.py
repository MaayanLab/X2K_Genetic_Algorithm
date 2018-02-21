## ********** FORMATTING DATABASE GMT FILES FOR X2K ********** ##

import os
import pandas as pd
os.chdir("../")

# HomoloGene
taxid = pd.read_table("General_Resources/Homologue_Mapping/taxid_taxname.dms")
homolo = pd.read_table("General_Resources/Homologue_Mapping/homologene.data",'\t')
homolo.head()
homolo.columns = ["HID","TaxonID","GeneID","GeneSymbol","ProteinGI","ProteinAC"]
human = homolo[homolo.TaxonID==9606]


# -----------[1] CREEDS Manual -----------
# Get list of TFs
with open("TF_datasets/ARCHS4_human/archs4_transcription_factor_gene_set_2017_08.gmt") as ARCHS4_TF_human:
    TFs=[]
    for line in ARCHS4_TF_human:
        TFs.append(line.split('\t')[0])
# Filter to just TFs perturbation experiments and write to new file
with open("General_Resources/CREEDS/Gene_perturbations/Preprocessed_manual_from_Enrichr/Single_Gene_Perturbations_from_GEO_up.txt") as up,\
 open("General_Resources/CREEDS/Gene_perturbations/Preprocessed_manual_from_Enrichr/Single_Gene_Perturbations_from_GEO_down.txt") as dn,\
 open("General_Resources/CREEDS/Gene_perturbations/manual_single_gene_perturbations_TF.gmt", "w") as newGMT_TF:
    dnLines = dn.readlines()
    for i,line in enumerate(up):
        lineSp = line.split("\t")
        name = lineSp[0].split(" ")
        gene = name[0]
        DEGs = lineSp[2:]
        if gene in TFs:
            # Up line
            newLine_up = "_".join(name)+"-up"+"\t\t"+"\t".join(DEGs)
            newGMT_TF.write(newLine_up)
            # Dn line
            lineSp_dn= dnLines[i].split('\t')
            newName_dn = "_".join( lineSp_dn[0].split(" "))+"-dn"
            DEGs_dn = lineSp_dn[2:]
            newGMT_TF.write(newName_dn+'\t'+"\t".join(DEGs_dn))

# Filter CREEDS automatic experiments by KINASES
with open("data/KEA/kea_ranks.txt") as KEAfile:
    KINASES = []
    KEAlines = KEAfile.readlines()
    for line in KEAlines:
        kinase = line.split("\t")[0]
        KINASES.append(kinase)
with open("General_Resources/CREEDS/Manual_single_gene_perturbations-v1.0.csv") as CREEDS_auto_gene_CONV, \
        open("Kinase_datasets/CREEDS_Kinases/CREEDS_manual_single_gene_pert_KINASES.txt", "w") as newFile:
    for line in CREEDS_auto_gene_CONV:
        lineSp = line.split('\t')
        if type(lineSp) == str:
            DEGS = lineSp[2:]
        gene = lineSp[0].split('-')[0]
        if gene in KINASES:
            newFile.write(line)

# -----------[2] CREEDS Automatic TF-----------
"""
# Manual Single-gene perturbations (Biased)
# Import metadata (have to use MANUAL since the Automatic metadata is missing human homologues)
CREEDS_gene_meta = pd.read_csv("General_Resources/CREEDS/Gene_perturbations/Manual_single_gene_perturbations-v1.0.csv")
CREEDS_gene_meta.head()
# Import data
#CREEDS_auto_gene = pd.read_table("Gene_perturbations/manual_single_gene_perturbations-v1.0.gmt", sep='\t', lineterminator='\r', header=None)

with open("General_Resources/CREEDS/Gene_perturbations/automatic_single_gene_perturbations-p1.0.gmt") as file:
    CREEDS_auto_gene = file.readlines()

## Write new converted file with only human GeneSymbols
### Use CREEDS api to get addiitonal metainformation. Due to the automated nature of the dataset, there's no consistent gene symboels
### However there are Uniprot IDs for many entries, which can be extracted using the following foramt:
### http://amp.pharm.mssm.edu/CREEDS/api?id=gene:P12317
## To convert any mouse IDs, need to use the complete NCBI HomloGene data, because the Manual CREEDS metadata only has several thousand genes
with open("General_Resources/CREEDS/Gene_perturbations/automatic_single_gene_perturbations.CONVERTED.txt", "w") as file:
    missingGenes=[]; editedLines=[]
    def addMissingGene(gene):
        if gene not in missingGenes:
            missingGenes.append(gene)
            print("Could not find gene conversion. (" + gene + " : # " + str(len(missingGenes)) + ")")
    # If the gene is from mouse, convert to human version
    def geneConverter(gene):
        if gene.startswith("gene:") and gene in list(CREEDS_gene_meta.id):
            newGene = CREEDS_gene_meta[CREEDS_gene_meta.id==gene].hs_gene_symbol.values[0]
        # If it's a human gene, use original
        elif gene in list(CREEDS_gene_meta.hs_gene_symbol):# or gene.upper() == gene:  # Add .upper comparison in case there's some human genes missing from HomoloGene
               #!!!!!!! newGene = gene
        # If it's an animal gene, see if there's a human homologue
        elif gene in list(CREEDS_gene_meta.mm_gene_symbol):
               #!!!!!!!!! newGene = homolo[homolo.mm_gene_symbol==gene].hs_gene_symbol.values[0]
        else:
            newGene = ""
            addMissingGene(gene)
        return newGene
    
    # Loop through all lines in CREEDS file
    for line in CREEDS_auto_gene:
        lineSp = line.split("\t")
        geneP = lineSp[1]
        up_dn = lineSp[0].split("-")[-1]
        ## Add only the first target to new line name
        targetGene_conv = geneConverter(geneP)
        # Convert Differentially Expressed Genes (DEGs)
        DEGs =lineSp[2:]
        if len(DEGs)>0:
            DEGs[-1] = DEGs[-1].strip()
        DEGs_conv=[]
        for d in DEGs:
            DEGs_conv.append( geneConverter(d) )
            DEGs_conv = list(filter(None, DEGs_conv))  # Remove any genes that couldn't be converted ("")
        # Write line to new file only if there's a converted target gene symbol
        if targetGene_conv != "":
            newLine = targetGene_conv + "-" + up_dn + '\t\t' + '\t'.join(DEGs_conv) + '\n'
            editedLines.append(newLine)
            file.write(newLine)

# Filter CREEDS experiments by TFs
with open("General_Resources/CREEDS/Gene_perturbations/automatic_single_gene_perturbations.CONVERTED.txt") as CREEDS_auto_gene_CONV,\
 open("TF_datasets/CREEDS_TF/CREEDS_manual_single_gene_pert_TFs.txt", "w") as newFile:
    for line in CREEDS_auto_gene_CONV:
        lineSp = line.split('\t')
        if type(lineSp)==str:
            DEGS = lineSp[2:]
        gene = lineSp[0].split('-')[0]
        if gene in TFs:
            newFile.write(line)
"""

