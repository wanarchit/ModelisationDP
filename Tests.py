import mdtraj as md

Data_pathP = 'D:\Cours\Master 2\Modelisation Bioinformatique\p3_p4_p5_p8'
Data_pathD = ''
t = md.load(Data_pathP+'\md_200ns_OK.xtc', top=Data_pathP+'\start.pdb')
print(t)
print("tototest2")