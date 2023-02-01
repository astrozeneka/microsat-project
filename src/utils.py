import math
from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt


chromosome_orders = ['01', '02', '03', '04', '05', '06', '07', '08', '09',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99',
                     'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'W', 'X', 'Y', 'Z', 'MIC01', 'MIC02', 'MIC03', 'MIC04', 'MIC05', 'MIC06', 'MIC07', 'MIC08', 'MIC09',
    'MIC10', 'MIC11', 'MIC12', 'MIC13', 'MT']
pastel_colors = ['#FD8A8A', '#A8D1D1', '#9EA1D4', '#ADA2FF', '#FFCAC8',
          "#DBA39A", '#F1F7B5', '#FCDDB0']  # One colors for each microsatellites


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

def search_ssr(records, config):
    MONO, DI, TRI, TETRA, PENTA, HEXA = config

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

def analyse_ratio(tsv_file, SSRS): # Get the ratio

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


def analyze_length_ratio(tsv_filename, SSRS, GENOME):
    chrom_name = tsv_filename[tsv_filename.rfind(GENOME):-4]
    fasta_file = f"../data/fasta/{GENOME}/{chrom_name}.fasta"
    chrom_len = len(SeqIO.read(fasta_file, "fasta"))
    with open(f"{tsv_filename}") as f:
        data = f.read().split("\n")
        data = [t.split("\t") for t in data]
        data = [t for t in data if len(t) > 1]
        total_len = len(data)
        output = {k: len([l for l in data if l[2]==k])/chrom_len * 100 for k in SSRS}
    return output

def microsat_length_report(chromosome, ssr_number=math.inf, start=None, end=None):
    with open(f"../data/tmp/{chromosome}.tsv") as f:
        microsat_list = f.read().split("\n")
        microsat_list = [t.split("\t") for t in microsat_list]
        microsat_list = [t for t in microsat_list if len(t) > 1]  # It is already sorted by the previous algorithm
        microsat_list = [(int(t[0]), t[1], t[2], t[3], int(t[4]), int(t[5]), int(t[6]), int(t[7]), int(t[8])) for t in
                         microsat_list]
        if start is not None:
            microsat_list = [t for t in microsat_list if t[6]>=start]
        if end is not None:
            microsat_list = [t for t in microsat_list if t[7]<=end]

    output = {}
    for l in microsat_list:
        if l[2] not in output.keys():
            output[l[2]] = l[8]
        else:
            output[l[2]]+= l[8]
    output = [(k,output[k]) for k  in output.keys()]
    output.sort(key=lambda x:x[1], reverse=True)
    output = output[:min(ssr_number, len(output))]

    # Normalize
    total_ssr_len = sum([x[1] for x in output])
    output = [(x[0],x[1]/total_ssr_len*100) for x in output]
    return output


def microsat_occ_report(chromosome, ssr_number=math.inf, start=None, end=None):
    with open(f"../data/tmp/{chromosome}.tsv") as f:
        microsat_list = f.read().split("\n")
        microsat_list = [t.split("\t") for t in microsat_list]
        microsat_list = [t for t in microsat_list if len(t) > 1]  # It is already sorted by the previous algorithm
        microsat_list = [(int(t[0]), t[1], t[2], t[3], int(t[4]), int(t[5]), int(t[6]), int(t[7]), int(t[8])) for t in
                         microsat_list]
        if start is not None:
            microsat_list = [t for t in microsat_list if t[6]>=start]
        if end is not None:
            microsat_list = [t for t in microsat_list if t[7]<=end]

    output = {}
    for l in microsat_list:
        if l[2] not in output.keys():
            output[l[2]] = 1
        else:
            output[l[2]]+= 1
    output = [(k,output[k]) for k  in output.keys()]
    output.sort(key=lambda x:x[1], reverse=True)
    output = output[:min(ssr_number, len(output))]

    # Normalize
    total_microsat = len(microsat_list)
    output = [(x[0],x[1]/total_microsat *100) for x in output]
    return output