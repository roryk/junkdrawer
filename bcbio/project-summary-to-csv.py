from __future__ import print_function
from argparse import ArgumentParser
import yaml

if __name__ == "__main__":
    parser = ArgumentParser(description="Convert metrics from a bcbio run to CSV.")
    parser.add_argument("yamlfile",
                        help="project-summary.yaml from a bcbio project")
    args = parser.parse_args()

    with open(args.yamlfile) as in_handle:
        dat = yaml.load(in_handle)
    summaries = [x["summary"] for x in dat["samples"]]
    metrics = [x["metrics"] for x in summaries]
    samples = [x["description"] for x in dat["samples"]]
    header = ["sample"] + metrics[0].keys()

    print(",".join(header))
    for i, sample in enumerate(samples):
        print(",".join(map(str, [sample] + metrics[i].values())))
