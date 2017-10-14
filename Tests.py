from Tkinter import Tk
from Tkinter import *
from tkFileDialog import askopenfilename
import sys


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





print("Welcome in the program !\n")
print("  - First, you will have to choose your data (PDB and XTC files")
print("  - Next, you will have to choose a reference atom")
print("  - Finally, you will have to choose a slice of distance to build the clusters\n")
print("Press a key to continue")
raw_input(" >> ")

