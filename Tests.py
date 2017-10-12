import mdtraj as md

t = md.load('Data\md_200ns_OK.xtc', top='Data\start.pdb')
print(t)
print("tototest")