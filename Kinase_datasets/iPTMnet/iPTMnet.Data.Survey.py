
#iPTMnet
with open("protein.txt") as prot_file:
    prot = prot_file.readlines()


new_genenames = []
for line in prot:
    genenames = line.split("\t")[3].split(" ")
    geneList = []
    for item in genenames:
        if item.endswith(";"):
            geneList.append(item.strip(";"))
    new_genenames.append(",".join(geneList))


with open("protein_edit.txt",'w') as prot_edit:
    for i,geneList in enumerate(new_genenames):
        prot_edit.write(prot[i].strip()+'\t'+geneList+'\n')


