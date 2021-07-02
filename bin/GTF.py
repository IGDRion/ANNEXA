#!/usr/bin/env python3


class GTFRecord:
    def __init__(self, line):
        # Remove comment
        line_splitted = line.split("#")[0].rstrip().split("\t")
        if len(line_splitted) != 9:
            exit(f"{line}\nUnable to parse this line, maybe badly formatted")
        (
            self.seqname,
            self.source,
            self.feature,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
        ) = line_splitted[:-1]
        self.start = int(self.start)
        self.end = int(self.end)
        self._parse_attributes_as_dict(line_splitted[-1])

    def __contains__(self, attribute):
        return attribute in self.attributes

    def __str__(self):
        begin = "\t".join(
            [
                self.seqname,
                self.source,
                self.feature,
                str(self.start),
                str(self.end),
                self.score,
                self.strand,
                self.frame,
            ]
        )
        end = " ".join(
            [f'{attribute} "{value}";' for attribute, value in self.attributes.items()]
        )
        return begin + "\t" + end

    def __len__(self):
        return abs(self.end - self.start)

    def __getitem__(self, item):
        return self.attributes[item]

    def __setitem__(self, key, value):
        self.attributes[key] = value

    def __delitem__(self, key):
        del self.attributes[key]

    def _parse_attributes_as_dict(self, att):
        self.attributes = dict(
            i.split("=")
            for i in att[:-2].replace('"; ', "|").replace(' "', "=").split("|")
        )

    def remove_attributes(self, attributes):
        """Remove attributes passed in arg (as list)"""
        for attribute in list(self.attributes.keys()):
            if attribute in attributes:
                del self.attributes[attribute]

    def filter_attributes(self, attributes):
        """Only keep attributes passed as arg (list)"""
        for attribute in list(self.attributes.keys()):
            if attribute not in attributes:
                del self.attributes[attribute]


class GTFRecordWithChildren(GTFRecord):
    def __init__(self, line):
        super().__init__(line)
        self.children = []

    def add_child(self, child, check_position=False):
        if check_position:
            if child.start < self.start:
                self.start = child.start
            if child.end > self.end:
                self.end = child.end
        self.children.append(child)

    def format_to_gtf(self):
        gtf_seq = str(self) + "\n"
        for child in self.children:
            if isinstance(child, GTFRecordWithChildren):
                gtf_seq += child.format_to_gtf() + "\n"
            elif isinstance(child, GTFRecord):
                gtf_seq += str(child) + "\n"
        return gtf_seq.rstrip()


class Gene(GTFRecordWithChildren):
    def __init__(self, line, check_attributes=False):
        super().__init__(line)
        self.feature = "gene"
        if check_attributes:
            to_filter = [
                att for att in self.attributes if ("exon" in att or "transcript" in att)
            ]
            self.remove_attributes(to_filter)

    @property
    def transcripts(self):
        return self.children

    @property
    def exons(self):
        return [exon for transcript in self.transcripts for exon in transcript.exons]


class Transcript(GTFRecordWithChildren):
    def __init__(self, line, check_attributes=False):
        super().__init__(line)
        self.feature = "transcript"
        if check_attributes:
            to_filter = [att for att in self.attributes if "exon" in att]
            self.remove_attributes(to_filter)

    @property
    def exons(self):
        return self.children


class GTF:
    @staticmethod
    def parse(fd, feature=None, by_line=False):
        if by_line:
            out = []
            for line in fd:
                if line.startswith("#"):
                    continue
                line = line.split("#")[0]
                record = GTFRecord(line)
                if record.feature == feature or feature is None:
                    out.append(record)
            return out
        else:
            genes = {}
            transcripts = {}
            other = []

            for line in fd:
                if line.startswith("#"):
                    continue

                # Remove comments
                line = line.split("#")[0]

                record = GTFRecord(line)
                if record.feature == "gene" and record["gene_id"] not in genes:
                    genes[record["gene_id"]] = Gene(line)
                elif record.feature == "transcript":
                    transcripts[record["transcript_id"]] = Transcript(line)
                elif record.feature == "exon":
                    other.append(record)

            for child in other:
                transcripts[child["transcript_id"]].add_child(child)

            for transcript in transcripts.values():
                genes[transcript["gene_id"]].add_child(transcript)

            return genes

    @staticmethod
    def reconstruct_full_gtf(file):
        genes = {}
        transcripts = {}

        # Have to do 2 loops if exons not sorted
        for record in GTF.parse(file, by_line=True):
            if record["gene_id"] not in genes:
                gene = Gene(str(record), check_attributes=True)
                genes[record["gene_id"]] = gene

            if record["transcript_id"] not in transcripts:
                transcript = Transcript(str(record), check_attributes=True)
                transcripts[record["transcript_id"]] = transcript

            transcripts[record["transcript_id"]].add_child(record, check_position=True)

        for transcript in transcripts.values():
            genes[transcript["gene_id"]].add_child(transcript, check_position=True)

        for gene in genes.values():
            yield gene

    @staticmethod
    def stats(file):
        exons = 0
        transcripts = set()
        genes = set()
        for exon in GTF.parse(file, feature="exon", by_line=True):
            exons += 1
            genes.add(exon["gene_id"])
            transcripts.add(exon["transcript_id"])
        return len(genes), len(transcripts), exons


##################################################
if __name__ == "__main__":
    import sys
    import os
    import argparse

    parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
    parser.add_argument(
        "mode",
        choices=["stats", "format"],
        type=str,
        help="Basic stats about your file | Format a gtf to with exon lines to gene and transcript levels",
    )
    parser.add_argument(
        "-i",
        "--input-file",
        help="Path to your GTF file. Use stdin by default",
        type=argparse.FileType("r"),
        default=(None if sys.stdin.isatty() else sys.stdin),
    )

    args = parser.parse_args()

    if args.input_file is None:
        print(
            "\033[91mPlease specify your GTF file or use stdin... See below for usage:\n\x1b[0m"
        )

        sys.exit(parser.print_help())

    if args.mode == "format":
        for gene in GTF.reconstruct_full_gtf(args.input_file):
            print(gene.format_to_gtf())

    elif args.mode == "stats":
        genes, transcripts, exons = GTF.stats(args.input_file)
        if args.input_file.name != "<stdin>":
            print(f"FILE: {os.path.abspath(args.input_file.name)}")
        print(f"# genes:\t{genes}\n# transcripts:\t{transcripts}\n# exons:\t{exons}")
