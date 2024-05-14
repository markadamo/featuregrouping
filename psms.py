import csv
from sys import argv
from matplotlib import pyplot as plt
from collections import defaultdict
import statistics
import numpy as np

DELIMITER = '\t'
TOP_N = 100

label_symbols = ['!', '@', '$']

mtraq_base = 140.0949630177
mtraq_delta = 4.0070994
mtraq_deltas = {None:-mtraq_base,
                '!':0,
                '@':mtraq_delta,
                '$':2*mtraq_delta
                }

def get_psm_groups(csvfile, label_symbols=label_symbols):
    reader = csv.DictReader(csvfile, delimiter=DELIMITER)

    zseq2psms = defaultdict(list)
    for row in reader:
        label = None
        num_tags = 0
        for symbol in label_symbols:
            if symbol in row['sequence']:
                label = symbol
                num_tags = row['sequence'].count(symbol)
        row['label'] = label
        row['num_tags'] = num_tags
        zseq2psms[(int(row['charge']), row['bare_sequence'])].append(row)

    #remove PSM groups that do not contain all labeled + unlabeled versions of the peptide
    for (z, seq) in list(zseq2psms.keys()):
        label_set = frozenset([psm['label'] for psm in zseq2psms[(z, seq)]])
        if label_set != frozenset(label_symbols + [None]):
            zseq2psms.pop((z, seq))

    return list(zseq2psms.values())

def psm_group_metrics(psm_groups):
    rt_ranges = []
    unlabeled_rt_diffs = []
    mz_errors = []

    for psm_group in psm_groups:
        labeled_rts = []
        unlabeled_rts = []
        tag_mzs = {None:[], '!':[], '@':[], '$':[]}

        for psm in psm_group:
            if psm['label'] is not None:
                labeled_rts.append(float(psm['rt']))
                charge = float(psm['charge'])
                num_tags = psm['num_tags']
            else:
                unlabeled_rts.append(float(psm['rt']))
            tag_mzs[psm['label']].append(float(psm['precursor_mz']))

        for tag in [None, '@', '$']:
            for mz0 in tag_mzs['!']:
                for mz in tag_mzs[tag]:
                    mz_error = (mz0 + mtraq_deltas[tag]*num_tags/charge) - mz
                    mz_errors.append(mz_error / statistics.mean(tag_mzs[tag]) * 1000000)

        rt_ranges.append(max(labeled_rts) - min(labeled_rts))
        unlabeled_rt_diffs.append(statistics.median(unlabeled_rts) - min(labeled_rts))

    rt_bins = np.linspace(0, 20, 100)

    plt.figure()
    plt.title('RT range for labeled peptides (s)')
    plt.xlabel('RT')
    plt.ylabel('count')
    plt.hist(rt_ranges, rt_bins)
    plt.figure()
    plt.title('Unlabeled RT vs. labeled RT diff (s)')
    plt.xlabel('RT')
    plt.ylabel('count')
    plt.hist(unlabeled_rt_diffs, 100)
    plt.figure()
    plt.title('M/z errors vs. nominal tag mass delta')
    plt.xlabel('M/z error (ppm)')
    plt.ylabel('count')
    plt.hist(mz_errors, 100)

    plt.show()

csvfile = open(argv[1])
psm_groups = get_psm_groups(csvfile)
psm_group_metrics(psm_groups)

if argv[2:]:
    outfile = open(argv[2], 'w')
    fieldnames = ['scan_num', 'rt', 'precursor_mz', 'precursor_scan_num', 'charge', 'sequence', 'bare_sequence', 'ppm', 'xcorr', 'dcn', 'label', 'num_tags', 'nterm_aa', 'sample_filename', 'refs', 'search', 'cterm_aa']
    writer = csv.DictWriter(outfile, fieldnames, delimiter=DELIMITER)
    writer.writeheader()
    for group in psm_groups:
        group.sort(key=lambda psm: str(psm['label']))
        for psm in group:
            writer.writerow(psm)
    outfile.close()
