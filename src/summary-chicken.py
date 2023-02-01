
from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = ['Arial', 'sans-serif']
matplotlib.rcParams['font.family'] = ['Arial', 'serif']

specie_name = "chicken"
GENOMES = ["GGA"]
SSRS = ["AG", "GAG"] # According  to papers
#SSRS = ['AG','AAG'] # According to algorithm


"""
    In order to demonstrate high ratio of AG microsatellite
    As previously studied by Matsubara et al. 2015
"""

MONO = 12
DI = 7
TRI = 5
TETRA = 4
PENTA = 4
HEXA = 4

chromosome_orders = ['01', '02', '03', '04', '05', '06', '07', '08', '09',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99',
                     'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'W', 'X', 'Y', 'Z', 'MT']
from utils import pastel_colors  # One colors for each microsatellites


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

def analyse_ratio(tsv_file): # Get the ratio

    with open(tsv_file) as f:
        data = f.read().split("\n")[1:]
        data = [t.split("\t") for t in data]
        data = [l for l in data if len(l) > 1]

        total_len = len(data)
        # Consider the "Standard" column
        # l[2] standard
        # l[3] motif
        if(total_len == 0):
            return {}
        output = {k: len([l for l in data if l[2]==k])/total_len * 100 for k in SSRS}

    return output

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

        # Plotting Microsatellite of interest
        print("Plotting microsatellite ratio")
        tsv_list = [g for g in glob("../data/tmp/*.tsv") if GENOME in g]
        tsv_list.sort(key=lambda x:chromosome_orders.index(x[x.index(GENOME)+len(GENOME):-4]))

        data = {d:[] for d in SSRS}
        chromosome_list = [x[x.index(GENOME):-4] for x in tsv_list]
        for chromosome_file in tsv_list:
            res = analyse_ratio(chromosome_file)
            for k in res.keys():
                data[k].append(res[k])
        barWidth = 1/(len(SSRS)+2)
        fig = plt.subplots(figsize=(17, 8))

        br = np.arange(len(data[k]))
        for i,k in enumerate(SSRS):
            brx = [x + barWidth * i - 0.5 + (barWidth * 1.5) for x in br]
            plt.bar(brx, data[k], color=pastel_colors[i], width=barWidth, label=k)

        plt.xlabel('Chromosomes')
        plt.ylabel('SSR ratio (%)')
        xticks = [f"{n}" for n in chromosome_list]
        plt.xticks([r for r in range(len(chromosome_list))], xticks, rotation=90)
        plt.title(f"SSRs ratio for {specie_name} {len(chromosome_list)} chromosomes")
        plt.legend(loc='upper left')
        plt.show()
        print("Done")