import argparse
import logging

from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--attr', dest='attrs', action='append') 
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    gtf_file = args.gtf_file
    comparisons = []
    for attr in args.attrs:
        key, op, value = attr.split()
        if (op == '=') or (op == '=='):
            func = lambda a,b: (a == b)
        elif (op == '!') or (op == '!='):
            func = lambda a,b: (a != b)
        else:
            assert False
        comparisons.append((key, value, func))
    for f in GTFFeature.parse(open(gtf_file)):
        match = True
        for k,v,func in comparisons:
            if not func(v, f.attrs.get(k, None)):
                match = False
                break
        if match:
            print str(f)

if __name__ == '__main__':
    main()