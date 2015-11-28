import pandas as pd
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("original")
    parser.add_argument("converted")
    args = parser.parse_args()

    original = pd.read_csv(args.original, delimiter="\t")
    original.columns = ['chromosome', 'start', 'end', 'gene', 'strand',
                        'annot1', 'annot2', 'alu?', 'non_alu_repetitive?',
                        'conservation_chimp', 'conservation_rhesus',
                        'conservation_mouse']
    converted = pd.read_csv(args.converted, delimiter="\t", header=None)
    converted.columns = ['orig_chr', 'orig_start', 'orig_end', 'skip',
                         'new_chr', 'new_start', 'new_end']
    merged = original.merge(converted, left_on=['chromosome', 'start'],
                            right_on=['orig_chr', 'orig_start'])
    merged = merged[['new_chr', 'new_start', 'new_end', 'gene', 'strand',
                     'annot1', 'annot2', 'alu?', 'non_alu_repetitive?',
                     'conservation_chimp', 'conservation_rhesus',
                     'conservation_mouse']]
    merged.columns = ['chromosome', 'start', 'end', 'gene', 'strand',
                        'annot1', 'annot2', 'alu?', 'non_alu_repetitive?',
                        'conservation_chimp', 'conservation_rhesus',
                        'conservation_mouse']
    import ipdb
    ipdb.set_trace()
