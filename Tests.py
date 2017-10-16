from Tkinter import Tk
from Tkinter import *
from tkFileDialog import askopenfilename
import sys
from subprocess import call


def loadData():
    val = raw_input("You want load a PDB file ? (Y/N) :")
    if val == "Y":
        select_file('pdb')
    val = raw_input("You want load a XTC file ? (Y/N) :")
    if val == "Y":
        select_file('xtc')
def select_file(format):
    root = Tk()
    filepath = askopenfilename(title="Select the "+format.upper()+" file",
                               filetypes=[(
        format.upper()+' files','.'+format.lower())])
    root.destroy()
    print filepath


var = "Toto"
var2 = "Tata"
script = "D:\\Cours\\Master 2\\Modelisation " \
          "Bioinformatique\\ModelisationDP\\testR.r"
Rscript = "D:\\Logiciels\\R\\R-3.3.2\\bin\\Rscript.exe"
toto = call([Rscript,script,var,var2], shell = False)


