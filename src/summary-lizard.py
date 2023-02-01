
from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from utils import microsat_length_report
matplotlib.rcParams['font.family'] = ['Arial', 'sans-serif']
matplotlib.rcParams['font.family'] = ['Arial', 'serif']

specie_name = "Sand Lizard"
GENOME = "LAG"
SSRS = [ # According to paper
    "AAGG",
    "ATAG", # Same as AGAT
    "AAAC",
    "ACAG",
    "AATC",
    "AAAAT",
    "AATCCC" # Same as TTAGGG
]
#SSRS = [ # According to algorithm
#    "ACTG","AAGG","C","ACCG","AAAT","AAG","CCG","A"
#]

MONO = 12
DI = 7
TRI = 5
TETRA = 4
PENTA = 4
HEXA = 4

from utils import chromosome_orders, pastel_colors, get_motif_standard, search_ssr, analyse_ratio, analyze_length_ratio


if __name__ == '__main__':

    fasta_list = glob(f"../data/fasta/{GENOME}/*.fasta")

    print(f"Fetching {GENOME} from fasta file")

    # No need to clear the temporary director
    # Clear temporary tmp
    # if os.path.isdir("../data/tmp"):
    #    shutil.rmtree("../data/tmp")
    #    print("Empty TMP DIRECTORY")
    # os.mkdir("../data/tmp")

    for fasta_file in fasta_list:
        records = SeqIO.parse(fasta_file, "fasta")
        kw = fasta_file[fasta_file.index(GENOME + "") + 4:-6]
        if os.path.isfile(f"../data/tmp/{kw}.tsv"):
            print(f"Skip {kw}.tsv, file already exists")
            continue

        # save it to tsv
        res = search_ssr(records, (MONO, DI, TRI, TETRA, PENTA, HEXA))
        with open(f"../data/tmp/{kw}.tsv", "w") as f:
            tsv_content = "\n".join(["\t".join([str(a) for a in b]) for b in res])
            f.write(tsv_content)
            print(f"Write tmp/{kw}.tsv")
    print("Done")

    # Plotting Microsatellite of interest
    print("Plotting microsatellite ratio")
    tsv_list = [g for g in glob("../data/tmp/*.tsv") if GENOME in g]
    tsv_list.sort(key=lambda x: chromosome_orders.index(x[x.index(GENOME) + len(GENOME):-4]))

    data = {d: [] for d in SSRS}
    chromosome_list = [x[x.index(GENOME):-4] for x in tsv_list]
    for chromosome_file in tsv_list:
        if True: # Ratio analysis
            res = analyse_ratio(chromosome_file, SSRS)
            for k in res.keys():
                data[k].append(res[k])
        #if True: # Length analysis
        #    res = analyze_length_ratio(chromosome_file, SSRS, "LAG")
        #    for k in res.keys():
        #        data[k].append(res[k])
    barWidth = 1/(len(SSRS)+2)
    fig = plt.subplots(figsize=(17, 8))

    br = np.arange(len(data[k]))
    for i, k in enumerate(SSRS):
        brx = [x + barWidth*i - 0.5+(barWidth*1.5) for x in br]
        plt.bar(brx, data[k], color=pastel_colors[i], width=barWidth, label=k)

    plt.xlabel('Chromosomes')
    plt.ylabel('SSR ratio (%)')
    xticks = [f"{n}" for n in chromosome_list]
    plt.xticks([r for r in range(len(chromosome_list))], xticks, rotation=90)
    plt.title(f"SSRs ratio for {specie_name} {len(chromosome_list)} chromosomes")
    plt.legend(loc='upper left')
    plt.show()
    print("Done")