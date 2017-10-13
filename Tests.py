import mdtraj as md

Data_pathP = 'D:\Cours\Master 2\Modelisation Bioinformatique\p3_p4_p5_p8'
Data_pathD = 'C:\Users\Asus\Documents\Master 2 GPhy\Semestre 1\Modelisation'
t = md.load(Data_pathP+'\md_200ns_OK.xtc', top=Data_pathP+'\start.pdb')
#print(t)

'''
residu = []
for res in t.topology.residues:
    if (str(res)[0:3] == 'DOP'):            // changer par res.name == 'DOP'
        residu.append(res)
        
>>> listElement = []
>>> for i in t.topology.atoms:
...     if i.element not in listElement:
...             listElement.append(i.element)
...
>>> listElement
[(6, 'carbon', 'C', 12.01078, 0.17), 
 (7, 'nitrogen', 'N', 14.00672, 0.155),
 (8, 'oxygen', 'O', 15.99943, 0.152), 
 (15, 'phosphorus', 'P', 30.9737622, 0.18), 
 (1, 'hydrogen', 'H', 1.007947, 0.12)]


for i in t.topology.atoms:
     if(i.name == 'P1'):
            print t.xyz[0,i.index]

'''

print("tototest2")
