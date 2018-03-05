import sys
import datetime
import numpy as np
import pandas as pd

def mapgenesymbols(inputDF, mappingDF, Column1, Column2):

    lst1 = []
    lst2 = []

    for i, index in enumerate(inputDF.index):

        progressPercent = ((i+1)/len(inputDF.index))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.index)))
        sys.stdout.flush()

        if inputDF.ix[index, Column1] in mappingDF.index:
            lst1.append(mappingDF.ix[inputDF.ix[index, Column1], 1])
        else:
            lst1.append(np.nan)

        if inputDF.ix[index, Column2] in mappingDF.index:
            lst2.append(mappingDF.ix[inputDF.ix[index, Column2], 1])
        else:
            lst2.append(np.nan)



    inputDF[Column1] = lst1
    inputDF[Column2] = lst2

    inputDF.dropna(inplace=True, subset=[Column1,Column2])


def getPubMedIds(inputDF, Column):

    lst = []

    for pub in inputDF[Column]:

        pub = str(pub)

        if '|' in pub:
            A = pub.split('|')
            for B in A:
                if 'pubmed' in B and 'unassigned' not in B:
                    C = B.split(':')[1]
                    break

        else:
            C = pub.split(':')[1]
        lst.append(C)
    inputDF[Column] = lst

def combineDupPPIs(inputDF):
    lst = []
    for index in inputDF.index:
        lst.append(tuple(sorted(tuple((inputDF.ix[index,'Protein A (gene name)'], inputDF.ix[index, 'Protein B (gene name)'])))))
    inputDF['ppi'] = lst

    inputDF.drop_duplicates(inplace=True)
    inputDF = inputDF.reset_index().drop('index', axis=1)

    if type(inputDF.ix[0, 'PubMed ID']) != str:
        inputDF['PubMed ID'] = inputDF['PubMed ID'].astype('str')

    for i,ppi in enumerate(inputDF['ppi'].unique()):

        progressPercent = ((i+1)/len(inputDF['ppi'].unique()))*100

        sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF['ppi'].unique())))
        sys.stdout.flush()

        ref = ('|').join(inputDF.ix[inputDF[inputDF['ppi'] == ppi].index, 'PubMed ID'])
        inputDF.ix[inputDF[inputDF['ppi'] == ppi].index[0], 'PubMed ID'] = ref
        inputDF.drop(inputDF[inputDF['ppi'] == ppi].index[1:], inplace=True)

def createSigFile(inputDB, filePath, name, filtered):

    sig_col = ['Source Name', 'Source Human Accession', 'Source Mouse Accession',
          'Source Type', 'Source Location', 'Target Name',
          'Target Human Accession', 'Target Mouse Accession', 'Target Type',
          'Target Location', 'Effect', 'Type of Interaction',
          'PubMed IDs']

    outputDB = pd.DataFrame(columns=sig_col)

    outputDB['Source Name'] = inputDB['Protein A (gene name)']

    outputDB['Target Name'] = inputDB['Protein B (gene name)']

    outputDB['PubMed IDs'] = inputDB['PubMed ID']

    outputDB.replace(np.nan, 'NA', inplace=True)

    if filtered:

        outputDB_ppiSIG = filePath+name+'_'+'filtered_ppi_%s.sig.zip'% str(datetime.date.today()).replace('-', '_')
        outputDB.to_csv(outputDB_ppiSIG, index=None, header=None, sep='\t', compression='gzip')

    else:

        outputDB_ppiSIG = filePath+name+'_'+'unfiltered_ppi_%s.sig.zip'% str(datetime.date.today()).replace('-', '_')
        outputDB.to_csv(outputDB_ppiSIG, index=None, header=None, sep='\t', compression='gzip')


def filterPPIbyPubmed(inputDB, level):

    for i,pub in enumerate(inputDB['PubMed ID'].unique()):

        progressPercent = ((i+1)/len(inputDB['PubMed ID'].unique()))*100

        sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDB['PubMed ID'].unique())))
        sys.stdout.flush()

        if inputDB[inputDB['PubMed ID'] == pub].shape[0] > level or 'unassigned' in inputDB[inputDB['PubMed ID'] == pub].values :

            inputDB.drop(inputDB[inputDB['PubMed ID'] == pub].index, inplace=True)


    # for i,value in enumerate(inputDB['PubMed ID']):
    #
    #     progressPercent = ((i+1)/len(inputDB['PubMed ID']))*100
    #
    #     sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDB['PubMed ID'])))
    #     sys.stdout.flush()
    #
    #     if 'unassigned' in value:
    #         lst.append(True)
    #     else:
    #         lst.append(False)
    #
    # inputDB = inputDB.drop(inputDB[lst].index)


def changePPIToShowGeneName(inputDB):

    lstA = []
    lstB = []

    for i,index in enumerate(inputDB.index):

        progressPercent = ((i+1)/len(inputDB['Protein A (gene name)']))*100

        sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDB['Protein A (gene name)'])))
        sys.stdout.flush()


        nameA = inputDB.ix[index, 'Protein A (gene name)']
        nameA = str(nameA)

        if '(gene name)' in nameA:
            temp = nameA.split('|')
            for ID in temp:
                if '(gene name)' in ID:
                    if '"' in ID:
                        lstA.append(ID.split(':')[1].split('(')[0].split('"')[1])
                        break
                    else:
                        lstA.append(ID.split(':')[1].split('(')[0])
                        break

        elif '(orf name)' in nameA:
            temp = nameA.split('|')
            for ID in temp:
                if '(orf name)' in ID:
                    lstA.append(ID.split(':')[1].split('(')[0])
                    break
        elif '(locus name)' in nameA:
            temp = nameA.split('|')
            for ID in temp:
                if '(locus name)' in ID:
                    lstA.append(ID.split(':')[1].split('(')[0])
                    break
        elif '(display_long)' in nameA:
            temp = nameA.split('|')
            for ID in temp:
                if '(display_long)' in ID:
                    lstA.append(ID.split(':')[1].split('(')[0].split('_')[0])
                    break
        elif nameA != 'nan':
            lstA.append(nameA.split('|')[-1].split(":")[1].split('(')[0])
        else:
            lstA.append(np.nan)

        nameB = inputDB.ix[index, 'Protein B (gene name)']
        nameB = str(nameB)

        if '(gene name)' in nameB:
            temp = nameB.split('|')
            for ID in temp:
                if '(gene name)' in ID:
                    if '"' in ID:
                        lstB.append(ID.split(':')[1].split('(')[0].split('"')[1])
                        break
                    else:
                        lstB.append(ID.split(':')[1].split('(')[0])
                        break
        elif '(orf name)' in nameB:
            temp = nameB.split('|')
            for ID in temp:
                if '(orf name)' in ID:
                    lstB.append(ID.split(':')[1].split('(')[0])
                    break
        elif '(locus name)' in nameB:
            temp = nameB.split('|')
            for ID in temp:
                if '(locus name)' in ID:
                    lstB.append(ID.split(':')[1].split('(')[0])
                    break
        elif '(display_long)' in nameB:
            temp = nameB.split('|')
            for ID in temp:
                if '(display_long)' in ID:
                    lstB.append(ID.split(':')[1].split('(')[0].split('_')[0])
                    break
        elif nameB != 'nan':
            lstB.append(nameB.split('|')[-1].split(":")[1].split('(')[0])
        else:
            lstB.append(np.nan)

    inputDB['Protein A (gene name)'] = lstA
    inputDB['Protein B (gene name)'] = lstB
