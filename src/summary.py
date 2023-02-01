
from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
import os
import shutil

GENOMES = ["GGA", "ACA", "NNA", "LAG"] # For all genomes
MONO = 12
DI = 7
TRI = 5
TETRA = 4
PENTA = 4
HEXA = 4

def get_motif_standard(motif):
    options = []
    for i in range(len(motif)):
        ex = motif[i:] + motif[:i]
        options.append(ex)
        options.append(Seq(ex).complement())
    def cost(motif):
        output = 0
        for i in range(len(motif)):
            output += 2**(len(motif)-i-1) * "ACTGN".index(motif[i])
        return (motif, output)
    cost_list = [cost(m) for m in options]
    standard = min(cost_list)[0]
    return standard

def search_ssr(records):

    output = []
    for seq_record in records:

        seq = str(seq_record.seq)
        motif = ""

        tmp = None
        rep = [MONO, DI, TRI, TETRA, PENTA, HEXA]

        seq_len = len(seq)
        i = 0
        while i < seq_len:
            if seq[i] == "N":
                i += 1 # important
                continue
            j = 1
            while j <= 6:
                start = i
                length = j
                while start+length<seq_len and seq[i]==seq[i+j] and seq[i]!="N":
                    i+=1
                    length+=1
                repeat = int(length/j)

                if(repeat>=rep[j-1]):
                    motif = seq[start:start+j]
                    length = repeat*j
                    tmp = (len(output)+1, seq_record.id, get_motif_standard(motif), motif, j, repeat, start+1, start+length, length)
                    output.append(tmp)
                    i = start+length
                    j = 0
                else:
                    i = start
                j+= 1
            i+=1
    return output

def plot():
    # In order to parse tsv file and to plot the SSRs ratio
    pass

if __name__ == '__main__':

    for GENOME in GENOMES: # we are here
        # List all chromosome fasta file
        fasta_list = glob(f"../data/fasta/{GENOME}/*.fasta")
        print(f"Running {GENOME}")

        # No need to clear the temporary director
        # Clear temporary tmp
        #if os.path.isdir("../data/tmp"):
        #    shutil.rmtree("../data/tmp")
        #    print("Empty TMP DIRECTORY")
        #os.mkdir("../data/tmp")

        for fasta_file in fasta_list:
            records = SeqIO.parse(fasta_file, "fasta")
            kw = fasta_file[fasta_file.index(GENOME + "")+4:-6]
            if os.path.isfile(f"../data/tmp/{kw}.tsv"):
                print(f"Skip {kw}.tsv, file already exists")
                continue
            res = search_ssr(records)

            # save it to tsv
            with open(f"../data/tmp/{kw}.tsv", "w") as f:
                tsv_content = "\n".join(["\t".join([str(a) for a in b]) for b in res])
                f.write(tsv_content)
                print(f"Write tmp/{kw}.tsv")
            print("Done")
