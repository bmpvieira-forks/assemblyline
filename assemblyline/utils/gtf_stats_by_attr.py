import argparse
import logging
import collections

from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--group-by", dest="group_by", default="gene_id")
    parser.add_argument("gtf_file")
    parser.add_argument("gtf_attr")
    args = parser.parse_args()
    d = collections.defaultdict(lambda: set())
    for f in GTFFeature.parse(open(args.gtf_file)):
        if args.gtf_attr in f.attrs:
            v = f.attrs[args.gtf_attr]
            f_id = f.attrs[args.group_by]
            d[v].add(f_id)
    for k in sorted(d):
        print k, len(d[k])

if __name__ == '__main__':
    main()
