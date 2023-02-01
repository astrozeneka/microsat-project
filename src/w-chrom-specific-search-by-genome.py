from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
matplotlib.rcParams['font.family'] = ['Arial', 'sans-serif']
matplotlib.rcParams['font.family'] = ['Arial', 'serif']
from utils import chromosome_orders, pastel_colors, get_motif_standard, search_ssr, analyse_ratio, \
    microsat_length_report, microsat_occ_report

PLOTTED_SSR_NUMBER = 8
GENOMES = ["LAG", "GGA"]
GENOMES =  ["NNA"]
CHROMOSOME_OF_INTEREST = "Z"
MAX_REPEAT_LEN = 4 # For only tetranucleotide repeats

MONO = 12
DI = 7
TRI = 5
TETRA = 4
PENTA = 4
HEXA = 4

def ssrs_report(chromosome, ssr_number):
    with open(f"../data/tmp/{chromosome}.tsv") as f:
        microsat_list = f.read().split("\n")
        microsat_list = [t.split("\t") for t in microsat_list]
        microsat_list = [t for t in microsat_list if len(t) > 1]
        microsat_list = [(int(t[0]), t[1], t[2], t[3], int(t[4]), int(t[5]), int(t[6]), int(t[7]), int(t[8])) for t in
                     microsat_list]


if __name__ == '__main__':
    for GENOME in GENOMES:
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
        print(f"{GENOME} Done")

    for GENOME in GENOMES:
        tsv_list = [g for g in glob("../data/tmp/*.tsv") if GENOME in g]
        tsv_list.sort(key=lambda x: chromosome_orders.index(x[x.index(GENOME) + len(GENOME):-4]))
        chromosome_list = [x[x.index(GENOME):-4] for x in tsv_list]

        data_indexes = {"A", "C"}
        reports = {}
        microsat_stats = {} # In each chromosome
        for chromosome_file in chromosome_list:
            report = microsat_occ_report(chromosome_file)
            #report = microsat_length_report(chromosome_file)
            reports[chromosome_file] = report
            data_indexes = data_indexes.union(set([x[0] for x in report]))
            data_indexes = {x for x in data_indexes if len(x) <= MAX_REPEAT_LEN}

        # TODO NEXT: Data frame manipulation
        # Transpose the previous matrix collected
        microsat_stats = {x:[] for x in data_indexes}
        for chrom in reports:
            report = reports[chrom]
            chrom_microsat = {a[0]:a[1] for a in report}
            for microsat in data_indexes:
                if microsat in chrom_microsat.keys():
                    microsat_stats[microsat].append(chrom_microsat[microsat])
                else:
                    microsat_stats[microsat].append(0)

        dataframe = pd.DataFrame(microsat_stats, index=chromosome_list) # No need to transpose

        # Select which microsatellite has highest percentage in W chromosome
        highest_ratio_microsatellites = [] # The final result
        for microsat in data_indexes:
            if max(dataframe[microsat]) == dataframe[microsat][f"{GENOME}{CHROMOSOME_OF_INTEREST}"]:
                # Should add value for sorting
                highest_ratio_microsatellites.append((microsat, dataframe[microsat][f"{GENOME}{CHROMOSOME_OF_INTEREST}"] / dataframe[microsat].mean()))
        highest_ratio_microsatellites.sort(key=lambda x:x[1], reverse=True)
        print(f"\nAnalysis done for {GENOME}:")
        print(f"{len(highest_ratio_microsatellites)} high ratio microsatellites discovered in {GENOME} :"
              f"{','.join([x[0] for x in highest_ratio_microsatellites])}")
        print("Plotting")

        ssrs = [x[0] for x in highest_ratio_microsatellites]
        ssrs = ssrs[:min(PLOTTED_SSR_NUMBER, len(ssrs))]
        barWidth = 1/(len(ssrs)+2)
        dataframe_t = dataframe.transpose()
        fig = plt.subplots(figsize=(17,8))
        br = np.arange(len(chromosome_list))
        for i,k in enumerate(ssrs):
            brx = [x + barWidth*i - 0.5+(barWidth*1.5) for x in br]
            _d = dataframe[k]
            _c = pastel_colors[i]
            plt.bar(brx, dataframe[k], color=pastel_colors[i], width=barWidth, label=k)

        plt.xlabel('Chromosomes')
        plt.ylabel('SSR ratio (%)')
        xticks = [f"{n}" for n in chromosome_list]
        plt.xticks([r for r in range(len(chromosome_list))], xticks, rotation=90)
        plt.title(f"SSRs ratio for {GENOME} {len(chromosome_list)} chromosomes")
        plt.legend(loc='upper left')
        txt = "I need the caption to be present a little below X-axis"
        plt.show()

