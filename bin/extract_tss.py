#!/usr/bin/env python3

from GTF import GTF

SEPARATOR = "\t"


def get_tss_interval(transcript, length):
    if transcript.strand == "+" or transcript.strand == ".":
        end = transcript.start + (length // 2)
        start = transcript.start - (length // 2) + 1
    else:
        start = transcript.end - (length // 2)
        end = transcript.end + (length // 2) - 1
    return start, end


def get_intervals(transcripts, length):
    intervals = []
    for transcript in transcripts:
        start, end = get_tss_interval(transcript, length)
        if start > 0:
            intervals.append(
                (
                    transcript.seqname,
                    start - 1,
                    end,
                    f"{transcript['gene_id']}::{transcript['transcript_id']}",
                    transcript.score,
                    transcript.strand,
                )
            )
        else:
            print(
                f"{transcript['transcript_id']} skipped: Start Coordinate detected < 0.",
                file=sys.stderr,
            )
    return intervals


if __name__ == "__main__":
    import argparse
    import sys

    # CLI
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--input",
        help="Path to your GTF file. Use stdin by default",
        type=argparse.FileType("r"),
        default=(None if sys.stdin.isatty() else sys.stdin),
    )
    parser.add_argument(
        "-l",
        "--length",
        help="",
        type=int,
        default=512,
    )
    args = parser.parse_args()

    # Extract all transcripts from input GTF
    transcripts = []
    for gene in GTF.parse(args.input):
        transcripts += gene.transcripts

    # Compute TSS for all intervals
    intervals = get_intervals(transcripts, args.length)

    # Print TSS interval as bed6 format for all individual transcript
    for interval in intervals:
        print(SEPARATOR.join(map(str, interval)))
