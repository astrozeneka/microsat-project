from Bio import SeqIO
import math
import random

"""
    This code has to demonstrate the Sand lizard Z chromosome
    has homology to chicken chromosome 9 and 6

    ====
"""

SSR_NUMBER = 20

def microsat_report(chromosome, ssr_number, start=None, end=None):
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
    output = output[:ssr_number]

    # Normalize
    total_ssr_len = sum([x[1] for x in output])
    output = [(x[0],x[1]/total_ssr_len) for x in output]
    return output

def partial_distance(partial1, partial2):
    output = 0
    for ssr_1 in partial1:
        hom = [d for d in partial2 if d[0] == ssr_1[0]]
        if len(hom) == 0:
            output+= ssr_1[1]**2
        else:
            output+= (ssr_1[1] - hom[0][1])**2
    return math.sqrt(output)

if __name__ == '__main__':
    # Find the most common microsat on GGA6 (all length)
    # GGA06, about 36 Mbp
    chicken_chrom_seq = SeqIO.read("../data/fasta/GGA/GGA06.fasta", "fasta")
    lizard_chrom_seq = SeqIO.read("../data/fasta/LAG/LAGZ.fasta", "fasta")
    #chicken_report = microsat_report("GGA06", ssr_number=SSR_NUMBER)

    #res = partial_similarity(chicken_report, lizard_report)
    #lizard_1 = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=0, end=chicken_chrom_len-STEP0)
    #lizard_2 = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=STEP0, end=chicken_chrom_len)
    #now = partial_similarity(chicken_report, lizard_report)
    #trunk_tail = partial_distance(chicken_report, lizard_2)
    #trunk_head = partial_distance(chicken_report, lizard_1)
    #lizard_report = microsat_report("LAGZ", ssr_number=SSR_NUMBER)
    #actual_distance = partial_distance(chicken_report, lizard_report)
    #head_margin = 0
    #tail_margin = 0

    # Gradient descent (4 dimension)
    # dim1: lizard head margin
    # dim2: lizard tail margin
    # dim3: chicken head margin
    # dim4: chicken tail margin


    chicken_chrom_len = len(chicken_chrom_seq)
    lizard_chrom_len = len(lizard_chrom_seq)

    STEP = 1000000

    lizard_head = random.randint(0,int(lizard_chrom_len/2/STEP))
    lizard_tail = random.randint(int(lizard_chrom_len/2/STEP), int(lizard_chrom_len/STEP))
    chicken_head = random.randint(0,int(chicken_chrom_len/2/STEP))
    chicken_tail = random.randint(int(chicken_chrom_len/2/STEP), int(chicken_chrom_len/STEP))

    lz = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=lizard_head*STEP, end=lizard_tail*STEP)
    ch = microsat_report("GGA06", ssr_number=SSR_NUMBER, start=chicken_head*STEP, end=chicken_tail*STEP)
    euclidian_distance = partial_distance(ch, lz)
    while True:
        print((lizard_head, lizard_tail, chicken_head, chicken_tail), euclidian_distance)
        lh_n = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=(lizard_head-1)*STEP, end=lizard_tail*STEP)
        lh_p = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=(lizard_head+1)*STEP, end=lizard_tail*STEP)
        lt_n = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=lizard_head*STEP, end=(lizard_tail-1)*STEP)
        lt_p = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=lizard_head*STEP, end=(lizard_tail+1)*STEP)
        ch_n = microsat_report("GGA06", ssr_number=SSR_NUMBER, start=(chicken_head-1)*STEP, end=chicken_tail*STEP)
        ch_p = microsat_report("GGA06", ssr_number=SSR_NUMBER, start=(chicken_head+1)*STEP, end=chicken_tail*STEP)
        ct_n = microsat_report("GGA06", ssr_number=SSR_NUMBER, start=chicken_head*STEP, end=(chicken_tail-1)*STEP)
        ct_p = microsat_report("GGA06", ssr_number=SSR_NUMBER, start=chicken_head*STEP, end=(chicken_tail+1)*STEP)

        # Gradient descent using Taxicab geometry (the simpliest)
        # https://en.wikipedia.org/wiki/Taxicab_geometry
        pgd_lh_n = partial_distance(ch, lh_n)
        pgd_lh_p = partial_distance(ch, lh_p)
        pgd_lt_n = partial_distance(ch, lt_n)
        pgd_lt_p = partial_distance(ch, lt_p)
        pgd_ch_n = partial_distance(ch_n, lz)
        pgd_ch_p = partial_distance(ch_p, lz)
        pgd_ct_n = partial_distance(ct_n, lz)
        pgd_ct_p = partial_distance(ct_p, lz)
        arr = [pgd_lh_n, pgd_lh_p, pgd_lt_n, pgd_lt_p, pgd_ch_n, pgd_ch_p, pgd_ct_n, pgd_ct_p]
        if min(arr) == pgd_lh_n:
            lizard_head = lizard_head-1
            euclidian_distance = pgd_lh_n
            continue
        if min(arr) == pgd_lh_p:
            lizard_head = lizard_head+1
            euclidian_distance = pgd_lh_p
            continue
        if min(arr) == pgd_lt_n:
            lizard_tail = lizard_tail-1
            euclidian_distance = pgd_lt_n
            continue
        if min(arr) == pgd_lt_p:
            lizard_tail = lizard_tail+1
            euclidian_distance = pgd_lt_p
            continue
        if min(arr) == pgd_ch_n:
            chicken_head = chicken_head-1
            euclidian_distance = pgd_ch_n
            continue
        if min(arr) == pgd_ch_p:
            chicken_head = chicken_head+1
            euclidian_distance = pgd_ch_p
            continue
        if min(arr) == pgd_ct_n:
            chicken_tail = chicken_tail-1
            euclidian_distance = pgd_ct_n
            continue
        if min(arr) == pgd_ct_p:
            chicken_tail = chicken_tail+1
            euclidian_distance = pgd_ct_p
            continue
        print("ERROR")
    print()

    """
    while True:
        _lizard_report = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=head_margin+STEP0, end=chicken_chrom_len)
        distance_after_head_trunk = partial_distance(chicken_report, _lizard_report)
        _lizard_report = microsat_report("LAGZ", ssr_number=SSR_NUMBER, start=head_margin, end=chicken_chrom_len-tail_margin-STEP0)
        distance_after_tail_trunk = partial_distance(chicken_report, _lizard_report)

        if min([actual_distance, distance_after_head_trunk, distance_after_tail_trunk]) == actual_distance:
            print()
            break
        if min([actual_distance, distance_after_head_trunk, distance_after_tail_trunk]) == distance_after_head_trunk:
            actual_distance = distance_after_head_trunk
            head_margin+=STEP0
            continue
        if min([actual_distance, distance_after_head_trunk, distance_after_tail_trunk]) == distance_after_tail_trunk:
            actual_distance = distance_after_tail_trunk
            tail_margin+=STEP0
            continue
    """

    print()
