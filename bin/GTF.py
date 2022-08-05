#!/usr/bin/env python3
from typing import Generator, Tuple, List, Union
import re


FIELDS_IDX = {
    "seqname": 0,
    "source": 1,
    "feature": 2,
    "start": 3,
    "end": 4,
    "score": 5,
    "strand": 6,
    "frame": 7,
}


class Attributes(dict):
    motif = r"(\w+) (\"){0,1}((?(2)[^;]*|\d+))(?(2)\"|);"  # 'key "str";' or 'key num;'
    regex = re.compile(motif)  # extract 1 key/value

    @classmethod
    def from_str(cls, attr_str: str):
        matches = cls.regex.findall(attr_str)
        if len(matches) != attr_str.count(";"):
            raise Exception(f"Unable to parse maybe misformatted attributes:\n{attr_str}")

        attr = {}
        for match in matches:
            attr[match[0]] = match[2]
        return cls(attr)

    def __str__(self) -> str:
        return " ".join([f'{attr} "{val}";' for attr, val in self.items()])

    def remove(self, attributes: List[str]):
        """Remove attributes passed in arg (as list)"""
        for attribute in list(self.keys()):
            if attribute in attributes:
                del self[attribute]

    def filter(self, attributes: List[str]):
        """Only keep attributes passed as arg (list)"""
        for attribute in list(self.keys()):
            if attribute not in attributes:
                del self[attribute]


class GtfObject(object):
    def __len__(self):
        return abs(self.end - self.start)

    def __getitem__(self, key: str):
        return self.attributes[key]

    def __setitem__(self, key: str, value):
        self.attributes[key] = value

    def __contains__(self, name: str):
        return name in self.attributes

    def __getattr__(self, name):
        if name in FIELDS_IDX:
            return self.fields[FIELDS_IDX[name]]
        return super().__getattribute__(name)

    def __setattr__(self, name, value):
        if name not in FIELDS_IDX:
            super().__setattr__(name, value)
        else:
            self.fields[FIELDS_IDX[name]] = value

    @property
    def attributes(self):
        if isinstance(self.fields[8], str):
            self.fields[8] = Attributes.from_str(self.fields[8])
        return self.fields[8]

    @attributes.setter
    def attributes(self, attributes: Attributes):
        self.fields[8] = attributes


class GtfRecord(GtfObject):
    def __init__(self, fields) -> None:
        object.__setattr__(self, "fields", fields)
        self.start = int(self.start)
        self.end = int(self.end)

    def __str__(self):
        return "\t".join(map(str, self.fields))

    @classmethod
    def from_line(cls, line: str):
        fields = line.rstrip().split("#", 7)[0].split("\t")
        if len(fields) != 9:
            raise Exception(f"Unable to parse line:\n{line}")
        return cls(fields)

    @classmethod
    def from_record(cls, record: "GtfRecord"):
        return cls(record.fields)


class GtfParent(GtfObject):
    def __init__(self) -> None:
        self.children = []

    def __getattr__(self, name: str):
        if name == "start":
            return min(child.start for child in self.children)
        if name == "end":
            return max(child.end for child in self.children)
        return self.first_child.__getattr__(name)

    @property
    def first_child(self):
        return next(iter(self.children))

    def add_child(self, child: Union[GtfRecord, "GtfParent"], id=None) -> None:
        if isinstance(self.children, list):
            self.children.append(child)
        elif isinstance(self.children, GtfChildren):
            if id is None:
                raise Exception("id should be given when add_child is used on Gene")
            self.children[id] = child
        else:
            raise Exception("Error with children")

    def get_attributes(self, filters=[]) -> Attributes:
        attributes = Attributes()
        for k, v in self.first_child.attributes.items():
            if all([filter not in k for filter in filters]):
                attributes[k] = v
        return attributes

    def to_record(self) -> GtfRecord:
        fields = {}
        for field in FIELDS_IDX:
            fields[field] = self.__getattr__(field)
        fields["attributes"] = self.attributes
        return GtfRecord(list(fields.values()))

    def format_to_gtf(self, filters=[]) -> str:
        lines = str(self.to_record()) + "\n"
        for child in self.children:
            if isinstance(child, GtfParent):
                lines += str(child.format_to_gtf(filters)) + "\n"
            elif isinstance(child, GtfRecord):
                lines += str(child) + "\n"
        return lines.rstrip()


class GtfTranscript(GtfParent):
    def __init__(self) -> None:
        self.children = []

    def __getattr__(self, name: str):
        if name == "feature":
            return "transcript"
        return super().__getattr__(name)

    @property
    def exons(self) -> List[GtfObject]:
        return [child for child in self.children if child.feature == "exon"]

    @property
    def attributes(self):
        return super().get_attributes(["exon"])


class GtfChildren(dict):
    def __iter__(self):
        return iter(self.values())


class GtfGene(GtfParent):
    def __init__(self) -> None:
        self.children = GtfChildren()

    def __getattr__(self, name):
        if name == "feature":
            return "gene"
        return super().__getattr__(name)

    @property
    def transcripts(self) -> GtfChildren:
        return self.children

    @property
    def exons(self) -> List[GtfObject]:
        return [exon for transcript in self.transcripts for exon in transcript.exons]

    @property
    def attributes(self):
        return super().get_attributes(["exon", "transcript"])


class GTF(dict):
    def __iter__(self):
        return iter(self.values())

    def add_record(self, record: GtfRecord):
        try:
            tx_id = record["transcript_id"]
            g_id = record["gene_id"]
        except:
            raise Exception(f"'transcript_id' or 'gene_id' not found in:\n{str(record)}")

        if g_id not in self:
            self[g_id] = GtfGene()
        if tx_id not in self[g_id].transcripts:
            self[g_id].add_child(GtfTranscript(), tx_id)
        self[g_id].transcripts[tx_id].add_child(record)

    @staticmethod
    def parse_by_line(fd, feature=None) -> Generator[GtfRecord, None, None]:
        for line in fd:
            if line.startswith("#"):
                continue
            line_wo_comment = line.split("#", 1)[0].rstrip()
            record = GtfRecord.from_line(line_wo_comment)
            if feature is None or record.feature == feature:
                yield record

    @classmethod
    def parse(cls, fd) -> "GTF":
        gtf = cls()
        for record in cls.parse_by_line(fd):
            if record.feature == "gene" or record.feature == "transcript":
                continue
            gtf.add_record(record)
        return gtf

    @staticmethod
    def stats(file) -> Tuple[int, int, int]:
        exons = 0
        transcripts = set()
        genes = set()
        for exon in GTF.parse_by_line(file, feature="exon"):
            exons += 1
            genes.add(exon["gene_id"])
            transcripts.add(exon["transcript_id"])
        return len(genes), len(transcripts), exons

    def write(self, out, level="gene") -> None:
        for gene in self:
            if level == "gene":
                out.write(gene.format_to_gtf() + "\n")
                continue

            for transcript in gene.transcripts:
                if level == "transcript":
                    out.write(transcript.format_to_gtf() + "\n")
                    continue

                for child in transcript.children:
                    out.write(str(child) + "\n")


##################################################
if __name__ == "__main__":
    import sys
    import os
    import argparse

    parser = argparse.ArgumentParser(description="Utility tools for GTF files.")
    parser.add_argument(
        "mode",
        choices=["stats", "format"],
        type=str,
        help="Basic stats about your file | Format a gtf to with exon lines to gene and transcript levels",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Path to your GTF file. Use stdin by default.",
        type=argparse.FileType("r"),
        default=(None if sys.stdin.isatty() else sys.stdin),
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to your output file. Use stdout by default.",
        type=argparse.FileType("w"),
        default=sys.stdout,
    )
    parser.add_argument(
        "--level",
        help="Level to format gtf. Values: gene, transcript or exon",
        choices=["gene", "transcript", "exon"],
        default="gene",
        type=str,
    )
    args = parser.parse_args()

    if args.input is None:
        print("\033[91mPlease specify your GTF file or use stdin... See below for usage:\n\x1b[0m")
        sys.exit(parser.print_help())

    if args.mode == "format":
        GTF.parse(args.input).write(args.output, level=args.level)

    elif args.mode == "stats":
        genes, transcripts, exons = GTF.stats(args.input)
        if args.input.name != "<stdin>":
            args.output.write(f"FILE: {os.path.abspath(args.input.name)}\n")
        args.output.write(f"# genes:\t{genes}\n# transcripts:\t{transcripts}\n# exons:\t{exons}\n")
