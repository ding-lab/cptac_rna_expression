import argparse
import csv
import gzip
from pathlib import Path


def main(args):
    sample_name = Path(args.tsv).stem
    with open(args.tsv) as f:
        comment = f.readline()
        reader = csv.DictReader(f, dialect='excel-tab')
        sample_col = reader.fieldnames[-1]

        print(comment.rstrip())
        cols = ['Geneid', 'Length', sample_col]
        print('\t'.join([*cols[:2], sample_name]))
        for row in reader:
            vals = [row[col] for col in cols]
            print('\t'.join(vals))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Shrink featureCounts output"
    )
    parser.add_argument('tsv', help="Path to featureCounts tsv")
    args = parser.parse_args()
    main(args)


