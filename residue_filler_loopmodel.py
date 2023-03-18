import urllib.request
import  glob, os
from modeller import *
from modeller.automodel import *
from Bio.Align.Applications import ClustalOmegaCommandline

def residue_filler_loopmodel(PDB):
    PDB = PDB.lower()
    '''
    Missing Residue Filler By @biomedical_informatics Edris Sharif Rahamni March 18, 2023
    '''
    # Downloading the reference sequence for given PDB ID
    url = f'https://www.rcsb.org/fasta/entry/{PDB}'
    urllib.request.urlretrieve(url, f'{PDB}.fasta.temp')
    with open(f'{PDB}.fasta.temp') as ref, open('alignment', 'w') as align:
        next(ref)
        print(f">P1;{PDB}_fill", file=align)
        for row, line in enumerate(ref):
            line = line.strip()
            if row == 0:
                print(line, file=align)
            else:
                print(line, end="", file=align)
        print("*", file=align)
    
    e = Environ()
    m = Model(e, file = PDB)
    aln = Alignment(e)
    aln.append_model(m, align_codes = PDB)
    aln.write(file = PDB + '.seq')
    # Preparing the sequence of the PDB file
    header = []
    with open(f'{PDB}.seq') as pdseq, open("PDB.temp", "w") as new:
        next(pdseq)
        for line in pdseq:
            line = line.strip()
            line = line.replace('*', '')
            if ":" not in line:
                print(line, file = new)
            else:
                header.append(line)
    # Creating an input for multiple sequence alignment
    files = ['PDB.temp', 'alignment']
    with open("alignment_input.temp", 'w') as outfile:
        for file in files:
            with open(file) as infile:
                outfile.write(infile.read())
    with open('alignment_output.temp', 'w') as al:
        pass
    clustalomega_cline = ClustalOmegaCommandline(infile = "alignment_input.temp",
                                             outfile = "alignment_output.temp",
                                             verbose = True,
                                             auto = True,
                                             force = True)
    clustalomega_cline()
    # Prepearing the alignment file 
    with open('alignment_output.temp') as al:
        for row,line in enumerate(al):
            if f'>P1;{PDB}_fill' in line:
                lastrow = row -1
    with open('alignment_output.temp') as al, open('alignment.align.temp', 'w') as align:
        for row,line in enumerate(al):
            line = line.strip()
            if row == 0:
                print(line, '\n', ''.join(header), file = align)
            if row > 0 and f'>P1;{PDB}_fill' not in line and row != lastrow:
                print(line, file = align)
            if f'>P1;{PDB}_fill' in line:
                print(line, '\n', 'sequence:::::::::', file = align)
            if row == lastrow:
                if '*' not in line:
                    line = line + '*'
                    print(line, file = align)
    # Modelling the missing residues
    log.verbose()
    env = Environ()
    env.io.atom_files_directory = ['.', '/atom_files']
    a = LoopModel(env, alnfile = 'alignment.align.temp',
                  knowns = f'{PDB}', sequence = f'{PDB}_fill')
    a.starting_model = 1
    a.ending_model = 1
    a.loop.starting_model = 1
    a.loop.ending_model = 2
    a.loop.md_level = refine.fast
    a.make()
    
    
    for file in glob.glob("*.temp"):
        os.remove(file)
    os.remove('alignment')
    os.remove(f'{PDB}.seq')
    for file in glob.glob("*"):
        if '.pdb' not in file and 'fill' in file:
            os.remove(file)
