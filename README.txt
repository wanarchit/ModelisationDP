			########################################
			Execution du script DiffusionMolecule.py
			########################################

Auteurs : Delphine Rousse et Paul Gand
Date : 16/10/2017


Avant-propos:
Ce script permet de calculer les distances parcourus par un atome de r�f�rence
dans diff�rents intervalles de temps.
Ces distances sont ensuite clusteris�es par groupe de distances choisis par l'utilisateur.

Ce READEME a pour but d'expliquer comment utiliser le script Python.


Pr�requis et installation:
Consulter la documentation "Manuel d'installation et pr�requis".


Execution du script:
USAGE = python.exe DiffusionMolecule.py [-h] fichier.pdb fichier.xtc Rscript.exe
- Ouvrir le PowerShell Windows
- Taper la commande en conformit� avec USAGE
- python.exe doit avoir le module MDtraj
- Les chemins sont tous des chemins absolus (PATH)
- Rscript.exe doit avoir les librairies GGPLOT2 et STRINGR d'install�es

R�sultats en sortie:
Les r�sultats sont enregistr�s dans le plus grand r�pertoire commun des deux fichiers
sources .pdb et .xtc.
- Le fichier CSV contenant la table de comptage (clusters de distances parcourues) � 
diff�rents intervalles de temps pour l'atome de r�f�rence s�lectionn�.
- Le graphique au format PNG