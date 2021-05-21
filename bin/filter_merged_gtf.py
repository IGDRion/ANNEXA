#! /usr/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('fai', type=str, help='path to fasta index file')
    parser.add_argument('gtf', type=str, help='path to gtf file')

    args = parser.parse_args()

    fai = {}
    with open(args.fai) as fd:
        for line in fd:
            line = line.split("\t")
            fai[line[0]] = int(line[1])

    with open(args.gtf) as fd:
        for line in fd:
            splitted = line.split("\t")
            if not(int(splitted[4]) > fai[splitted[0]] or int(splitted[3]) > fai[splitted[0]]):
                print(line.rstrip())

if __name__ == "__main__":
    main()
