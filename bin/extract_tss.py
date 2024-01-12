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


def get_interval_record(tx, length):
    start, end = get_tss_interval(tx, length)
    if start > 0:
        return (
            tx.seqname,
            start - 1,
            end,
            f"{tx['gene_id']}::{tx['transcript_id']}::{tx.strand}",
            tx.score,
            tx.strand,
        )
    else:
        print(
            f"{tx['transcript_id']} skipped: Start Coordinate detected < 0.",
            file=sys.stderr,
        )
    return None


def get_intervals(transcripts, length):
    intervals = []
    for transcript in transcripts:
        # If defined strand, test only TSS
        if transcript.strand in ("+", "-"):
            interval_record = get_interval_record(transcript, length)
            if interval_record is not None:
                intervals.append(interval_record)

        # If undefined strand, simulate and test both extremities
        else:
            for strand in ("+", "-"):
                transcript.strand = strand
                interval_record = get_interval_record(transcript, length)
                if interval_record is not None:
                    intervals.append(interval_record)

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
