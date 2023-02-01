from glob import glob

GENOME = "GGA"

SSRS = [
    "AAGG",
    "ATAG",
    "AAAC",
    "ACAG",
    "AATC",
    "AAAT",
    "AATCCC"
]

chromosome_orders = ['01', '02', '03', '04', '05', '06', '07', '08', '09',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99',
                     'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'W', 'X', 'Y', 'Z']

def analyze_ratio(filename):
    with open(f"../data/tmp/{filename}") as f:
        data = f.read().split("\n")
        data = [t.split("\t") for t in data]
        data = [t for t in data if len(t) > 1]
        total_len = len(data)
        output = {k: len([l for l in data if l[2]==k])/total_len * 100 for k in SSRS}
    return output


if __name__ == '__main__':

    # FOR EACH GENOME

    # Get tsv datalist
    fasta_list = glob(f"../data/tmp/{GENOME}_*.tsv")
    chroms = [fasta_file[fasta_file.index(GENOME + "_"):-4] for fasta_file in fasta_list]
    chroms.sort(key=lambda x:chromosome_orders.index(x[x.index("_")+1:]), reverse=False)
    print()