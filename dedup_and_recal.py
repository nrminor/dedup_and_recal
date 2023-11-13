#!/usr/bin/env python3

""" 
The module `dedup_and_recal.py` is meant to combine redundant
assembly contigs with identical consensus sequences and combine them.
It also has the ability to parse a read support/depth value at the end
of the FASTA defline and include a combined read support in any sequences
it deduplicates. 

To run the module, use this as your template:
```
python3 dedup_and_recal.py \
--fasta /path/to/fasta/file.fasta \
--output "desired_output_prefix"
```
"""

import argparse
from typing import Tuple, Dict, TextIO, Optional, List
from dataclasses import dataclass


@dataclass
class SeqWithSupport:
    """
    The `SeqWithSupport` dataclass stores a) each FASTA
    sequence, and b) the read-support (depth of coverage)
    for the contig represented in that FASTA.
    """

    sequence: str
    read_support: Optional[int]


def parse_command_line_args() -> Tuple[str, str]:
    """
    Parse two possible arguments provided with keyword flags in the command line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        required=True,
        help="FASTA file of assembly contig consensus sequences to deduplicate.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        default="deduplicated_contigs",
        help="Prefix for output FASTA file.",
    )
    args = parser.parse_args()

    return args.fasta, args.output


def _try_parse_int(value: str) -> Optional[int]:
    """
    Helper function that handles the possibility that a read support
    cannot be parsed as an integer from the FASTA defline and returns
    `None` instead of raising an unrecoverable error.
    """
    try:
        return int(value)
    except ValueError:
        return None


def generate_seq_dict(input_fasta: list[str]) -> Dict[str, SeqWithSupport]:
    """
    Generate a dictionary that will structure FASTA deflines as its keys and use the
    `SeqWithSupport` dataclass to store sequences and read supports as its values
    """

    deflines = [line for line in input_fasta if line.startswith(">")]
    supports = [
        _try_parse_int(line.split("_")[-1])
        for line in input_fasta
        if line.startswith(">")
    ]
    sequences: List[str] = []
    sequence_parts: List[str] = []
    for line in input_fasta:
        line = line.strip()  # Remove newline and any trailing whitespace
        if line.startswith(">"):
            if (
                sequence_parts
            ):  # If there's any sequence collected, join it and add to sequences
                sequences.append("".join(sequence_parts))
                sequence_parts = []  # Reset for the next sequence
        else:
            sequence_parts.append(line)
    if sequence_parts:
        sequences.append("".join(sequence_parts))

    seqs_and_support = [
        SeqWithSupport(sequence=seq, read_support=support)
        for seq, support in zip(sequences, supports)
    ]

    assert len(deflines) == len(
        seqs_and_support
    ), "Mismatch between the number of deflines and number of sequences"

    seq_dict = dict(zip(deflines, seqs_and_support))

    return seq_dict


def _make_new_defline(old_defline: str, support: Optional[int]) -> str:
    """
    Helper function that makes a new defline with updated read support.
    """

    if support is None:
        return old_defline

    new_defline_list = old_defline.split("_")[:-1] + ["combined"] + [str(support)]
    new_defline = "_".join(new_defline_list)

    return new_defline


def _write_wrapped(output_file: TextIO, sequence: str, wrap_at: int = 80) -> None:
    """
    Basic bioinformatic wheel I have reinvented that wraps FASTA sequences every
    80 characters.
    """

    buffer: List[str] = []
    for nucleotide in sequence:
        buffer.append(nucleotide)
        if len(buffer) == wrap_at:
            line = "".join(buffer)
            output_file.write(f"{line}\n")
            buffer.clear()
    if buffer:
        output_file.write("".join(buffer) + "\n")


def deduplicate_and_recalibrate(
    seq_struct: Dict[str, SeqWithSupport], output_fasta: TextIO
) -> None:
    """
    Function `deduplicate_and_recalibrate()` combines any contigs with identical sequences
    and recalibrates their previously separate read support by adding them.
    """
    # establish a bound for looping
    bound = len(seq_struct)

    # unpack the dictionary into indexable lists
    deflines, seqs = zip(*seq_struct.items())

    i = 0
    while i < bound - 1:
        sequence = seqs[i].sequence
        next_sequence = seqs[i + 1].sequence

        if sequence == next_sequence:
            support = seqs[i].read_support
            next_support = seqs[i + 1].read_support
            new_support = support + next_support
            new_defline = _make_new_defline(deflines[i], new_support)
            output_fasta.write(f"{new_defline}\n")
            _write_wrapped(output_fasta, sequence)
            i += 2
        else:
            output_fasta.write(f"{deflines[i]}")
            _write_wrapped(output_fasta, sequence)
            i += 1


def main() -> None:
    """
    Main coordinates the flow of data through the above-defined functions.
    """

    input_fasta, output_prefix = parse_command_line_args()

    with open(input_fasta, "r", encoding="utf-8") as old_file, open(
        f"{output_prefix}.fasta", "w", encoding="utf-8"
    ) as new_fasta:
        old_fasta = old_file.readlines()
        seq_dict = generate_seq_dict(old_fasta)
        deduplicate_and_recalibrate(seq_dict, new_fasta)


if __name__ == "__main__":
    main()
