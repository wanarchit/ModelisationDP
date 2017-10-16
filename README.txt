			########################################
			Execution du script DiffusionMolecule.py
			########################################

Auteurs : Delphine Rousse et Paul Gand
Date : 16/10/2017


Avant-propos:
Ce script permet de calculer les distances parcourus par un atome de référence
dans différents intervalles de temps.
Ces distances sont ensuite clusterisées par groupe de distances choisis par l'utilisateur.

Ce READEME a pour but d'expliquer comment utiliser le script Python.


Prérequis et installation:
Consulter la documentation "Manuel d'installation et prérequis".


Execution du script:
USAGE = python.exe DiffusionMolecule.py [-h] fichier.pdb fichier.xtc Rscript.exe
- Ouvrir le PowerShell Windows
- Taper la commande en conformité avec USAGE
- python.exe doit avoir le module MDtraj
- Les chemins sont tous des chemins absolus (PATH)
- Rscript.exe doit avoir les librairies GGPLOT2 et STRINGR d'installées

Résultats en sortie:
Les résultats sont enregistrés dans le plus grand répertoire commun des deux fichiers
sources .pdb et .xtc.
- Le fichier CSV contenant la table de comptage (clusters de distances parcourues) à 
différents intervalles de temps pour l'atome de référence sélectionné.
- Le graphique au format PNG