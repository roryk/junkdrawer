import pandas as pd
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("original")
    parser.add_argument("converted")
    args = parser.parse_args()

    original = pd.read_csv(args.original, delimiter="\t", dtype=str)
    original.columns = ['chromosome', 'start', 'end', 'gene', 'strand',
                        'annot1', 'annot2', 'alu?', 'non_alu_repetitive?',
                        'conservation_chimp', 'conservation_rhesus',
                        'conservation_mouse']
    converted = pd.read_csv(args.converted, delimiter="\t", header=None, dtype=str)
    converted.columns = ['orig_chr', 'orig_start', 'orig_end', 'skip',
                         'new_chr', 'new_start', 'new_end']
    merged = original.merge(converted, left_on=['chromosome', 'start'],
                            right_on=['orig_chr', 'orig_start'], how="right")
    merged = merged[['new_chr', 'new_start', 'new_end', 'gene', 'strand',
                     'annot1', 'annot2', 'alu?', 'non_alu_repetitive?',
                     'conservation_chimp', 'conservation_rhesus',
                     'conservation_mouse']]
    merged.columns = ['chromosome', 'start', 'end', 'gene', 'strand',
                        'annot1', 'annot2', 'alu?', 'non_alu_repetitive?',
                        'conservation_chimp', 'conservation_rhesus',
                        'conservation_mouse']
    merged = merged.dropna()
    merged['chromosome'] = merged['chromosome'].str[3:]
    merged['start'] = merged['start'].astype(int)
    merged = merged.sort(["chromosome", "start"])
    merged['chromosome'] = 'chr' + merged['chromosome']
    merged.to_csv("out.bed", sep="\t", index=False)
