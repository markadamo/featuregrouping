import xml.etree.ElementTree as ET
import bisect
import csv
import statistics
import numpy as np
import argparse
from pyopenms import FeatureXMLFile, FeatureMap, Feature
from collections import defaultdict
from matplotlib import pyplot as plt

mtraq_base = 140.0949630177
mtraq_delta = 4.0070994
mtraq_deltas = {4:mtraq_delta,
                8:2*mtraq_delta
                }

class Group:
    def __init__(self, base_feature, n_tags):
        self.n_tags = n_tags
        self.charge = base_feature.getCharge()
        self.features = {None:None, 0:base_feature, 4:None, 8:None}
        self.psms = {None:[], 0:[], 4:[], 8:[]}

    @property
    def num_features(self):
        return sum(f is not None for f in self.features.values())
    @property    
    def score(self):
        return sum(f.getOverallQuality() for f in self.features.values() if f is not None)

def load_featureXML(path):
    map = FeatureMap()
    FeatureXMLFile().load(path, map)
    return map

def load_biosaur2(path):
    tsvfile = open(path)
    reader = csv.DictReader(tsvfile, delimiter='\t')
    features = []
    for row in reader:
        f = Feature()
        f.setRT(float(row['rtApex'])*60)
        f.setMZ(float(row['mz']))
        f.setCharge(int(row['charge']))
        #quality and width are unsupported for now
        f.setOverallQuality(0)
        f.setWidth(0)
        features.append(f)
    return features

parsers = {'featureXML': load_featureXML,
           'biosaur2': load_biosaur2}
        
def group_features(featuremap, mz_tol_ppm=4, rt_tol=8, max_tags_per_peptide=3):
    z_features = defaultdict(list)

    #partition features into mz-sorted lists by charge
    for f in featuremap:
        bisect.insort(z_features[f.getCharge()], f, key=lambda f: f.getMZ())

    groups = []

    #for each charge state...
    for z, flist in z_features.items():
        #iterate over features, anchoring each as a potential mTraq0 tagged peptide
        for (i, f) in enumerate(flist[:-1]):
            mz = f.getMZ()
            
            #consider that this peptide may be labeled with up to [max_tags_per_peptide]
            best_group = None
            for n_tags in range(1, max_tags_per_peptide+1):

                #calculate mz deltas between isotopes, based on tag count and charge
                tag_mz_deltas = {tag:(m*n_tags)/z for (tag, m) in mtraq_deltas.items()}

                #calculate ppm-corrected mz tolerance ranges
                mz_tol_ranges = {tag:((mz+mtd)-(mz+mtd)*mz_tol_ppm/1000000,
                                      (mz+mtd)+(mz+mtd)*mz_tol_ppm/1000000) for (tag, mtd) in tag_mz_deltas.items()}
                mz_upper_bound = max(max(mz_tol_ranges.values()))
                
                group = Group(f, n_tags)
                #iterate over subsequent M/z-sorted features to find higher isotopes
                for fn in flist[i+1:]:
                    #validate that this feature is within RT range and local M/z limit
                    if fn.getMZ() > mz_upper_bound:
                        break 
                    rt_diff = abs(f.getRT() - fn.getRT())
                    if rt_diff > rt_tol:
                        continue

                    #check if this feature is within mz range for any tag mass delta
                    for tag, (min_mz, max_mz)  in mz_tol_ranges.items():
                        if min_mz <= fn.getMZ() <= max_mz:
                            fg = group.features[tag]
                            if fg is None or fn.getOverallQuality() > fg.getOverallQuality():
                                #store as a tagged variant of the base feature
                                #(if another feature is already stored, take the one with the highest OverallQuality) 
                                group.features[tag] = fn
                                
                #discard the group if it does not contain all 3 tagged variants
                if sum(f is not None for f in group.features.values()) < 3:
                    continue

                #calculate unlabeled mz and search backward in the mz-sorted feature list
                unlabeled_mz = f.getMZ() - (mtraq_base*n_tags)/f.getCharge()
                mz_lower_bound = unlabeled_mz - unlabeled_mz*mz_tol_ppm/1000000
                for fn in flist[:i][::-1]:
                    if fn.getMZ() < mz_lower_bound:
                        break
                    if fn.getRT() > f.getRT(): continue
                    if abs(unlabeled_mz - fn.getMZ()) < fn.getMZ()*mz_tol_ppm/1000000:
                        fg = group.features[None]
                        #assume that features closer in mz will be more likely matches
                        if fg is None or abs(f.getMZ()-fn.getMZ()) < abs(f.getMZ()-fg.getMZ()):
                            group.features[None] = fn

                if best_group is None or group.num_features > best_group.num_features or group.score > best_group.score:
                    best_group = group

            if best_group is not None and best_group.num_features >= 3:
                groups.append(best_group)

    return sorted(groups, key=lambda g: g.score, reverse=True)

def write_tsv(feature_groups, outfile, write_psms=False):
    fieldnames = ['charge',
                  'mz_unlabeled', 'mz_mtraq0', 'mz_mtraq4', 'mz_mtraq8',
                  'rt_unlabeled', 'rt_mtraq0', 'rt_mtraq4', 'rt_mtraq8',
                  'quality_unlabeled', 'quality_mtraq0', 'quality_mtraq4', 'quality_mtraq8',
                  'width_unlabeled', 'width_mtraq0', 'width_mtraq4', 'width_mtraq8']
    if write_psms:
                  fieldnames += ['sequences_unlabeled', 'sequences_mtraq0', 'sequences_mtraq4', 'sequences_mtraq8',]
                  
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(fieldnames)
    for group in feature_groups:
            row = [group.charge,
                   *[f.getMZ() if f is not None else '' for f in group.features.values()],
                   *[f.getRT() if f is not None else '' for f in group.features.values()],
                   *[f.getOverallQuality() if f is not None else '' for f in group.features.values()],
                   *[f.getWidth() if f is not None else '' for f in group.features.values()],
                   *[','.join(set([psm['sequence'] for psm in psms])) for psms in group.psms.values()],
                         
                   ]
            if write_psms:
                row += [','.join(set([psm['sequence'] for psm in psms])) for psms in group.psms.values()],
            writer.writerow(row)
            
    outfile.close()
                      
def annotate_psms(feature_groups, psm_tsv_path, mz_tol_ppm=4, rt_tol=8, show_plots=False):
    #read PSMs and group by charge state
    psm_tsv = open(psm_tsv_path)
    reader = csv.DictReader(psm_tsv, delimiter='\t')
    charge2psms = defaultdict(list)
    for psm in reader:
        psm['label'] = None
        if '!' in psm['sequence']:
            psm['label'] = 0
        elif '@' in psm['sequence']:
            psm['label'] = 4
        elif '$' in psm['sequence']:
            psm['label'] = 8
        charge2psms[int(psm['charge'])].append(psm)

    unmatched = {None:0, 0:0, 4:0, 8:0}
    matched = {None:0, 0:0, 4:0, 8:0}
    mismatched = {None:0, 0:0, 4:0, 8:0}
    label2category = {'unmatched': unmatched,
                      'matched (correct tag)': matched,
                      'matched (incorrect tag)': mismatched}
    match3 = 0
    match4 = 0
    
    for group in feature_groups:
        #match PSMs to features by mz and RT tolerance
        for (tag, f) in group.features.items():
            if f is None: continue
            matched_psms = [psm for psm in charge2psms[f.getCharge()] if abs(float(psm['precursor_mz']) - f.getMZ()) < f.getMZ()*mz_tol_ppm/1000000 and abs(f.getRT() - float(psm['rt'])) < rt_tol]
            
            group.psms[tag] = matched_psms

            #record whether each match is correct with respect to the feature's putative labe
            if not matched_psms:
                unmatched[tag] += 1
            elif all(psm['label'] == tag for psm in matched_psms):
                matched[tag] += 1
            else:
                mismatched[tag] += 1

        #additionally check whether PSM-matched features are in concordance WRT sequence
        def check_psms(tags):
            bare_sequences = set()
            for tag in tags:
                psms = group.psms[tag]
                if not psms:
                    return False
                for psm in psms:
                    if psm['label'] != tag:
                        return False
                    bare_sequences.add(psm['bare_sequence'])
            return len(bare_sequences) == 1

        if check_psms([None, 0, 4, 8]):
            match4 += 1
        elif check_psms([0, 4, 8]):
            match3 += 1

    n_groups = len(feature_groups)
    print(f'total groups = {n_groups}')
    for (label, category) in label2category.items():
        print(f'{label}:')
        for (tag, count) in category.items():
            print(f'    tag={tag}: {count} {(count/len(feature_groups)*100):5.2f}%')

    print(f'3 concordant matches: {match3} {(match3/n_groups*100):5.2f}%')
    print(f'4 concordant matches: {match4} {(match4/n_groups*100):5.2f}%')

    if show_plots:
        fig, ax = plt.subplots()
        ax.set_title('PSM match distribution')
        ax.set_ylabel('features')
        ax.set_xlabel('tag')
        bottoms = np.zeros(4)
        width = 0.5
        for (category, color) in zip((unmatched, matched, mismatched), ('gray', 'dodgerblue', 'red')):
            labels = list(str(x) for x in category.keys())
            counts = list(category.values())
            ax.bar(labels, counts, width, bottom=bottoms, label=labels, color=color)
            bottoms += counts

        plt.show()
    
    return feature_groups

if __name__ == '__main__':
    import sys
    from logging import StreamHandler
    parser = argparse.ArgumentParser(
                    prog='features.py',
                    description='mTRAQ feature grouping')

    parser.add_argument('feature_file')
    parser.add_argument('-t', '--input_type', choices=('featureXML', 'biosaur2'), default='featureXML', help="Specify 'biosaur2' if the input file is a biosaur2 tsv")
    parser.add_argument('-p', '--psm_file', help='Optionally provide a PSM tsv file to perform QC. Adds matched peptide sequences to output.')
    parser.add_argument('-o', '--output_file', help='Output path for tsv file. Omit to output to stdout.')
    parser.add_argument('-ppm', '--mz_tol_ppm', type=float, default=4, help='M/z tolerance in ppm.')
    parser.add_argument('-rt', '--rt_tol', type=float, default=8, help='Retention time tolerance in seconds.')
    parser.add_argument('-nt', '--max_tags_per_peptide', type=int, default=3, help='Maximum number of mTRAQ tags to consider per peptide.')

    ns = parser.parse_args(sys.argv[1:])

    parser = parsers[ns.input_type]
    print('Reading features from input file...')
    features = parser(ns.feature_file)

    print('Grouping features...')
    feature_groups = group_features(features, mz_tol_ppm=ns.mz_tol_ppm, rt_tol=ns.rt_tol,
                                    max_tags_per_peptide=ns.max_tags_per_peptide)

    if ns.psm_file:
        print('Annotating features with PSM matches...')
        annotate_psms(feature_groups, ns.psm_file, mz_tol_ppm=ns.mz_tol_ppm, rt_tol=ns.rt_tol)
    
    if ns.output_file:
        outfile = open(ns.output_file, 'w')
    else:
        outfile = sys.stdout

    write_tsv(feature_groups, outfile, write_psms=bool(ns.psm_file))

    

    
    
