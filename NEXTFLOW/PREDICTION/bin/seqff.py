#!/usr/bin/env python

# import subprocess
# import argparse
# import os
# import pandas as pd

# r_exec = 'Rscript'

# current_dir = os.path.dirname(os.path.abspath(__file__))
# assets_dir = os.path.join(current_dir, "../assets")
# r_script = os.path.join(assets_dir, "seqff.r") 
# seqffdir = assets_dir

# def execute_seqff(sample_path):
#     """
#     Call SeqFF Rscript to compute ff for one BAM file.

#     :param sample_path: Path to a single BAM file
#     :return: Floating-point FF value (as string)
#     """
#     # Call Rscript with BAM path and seqffdir
#     cmd = [r_exec, r_script, sample_path, seqffdir]

#     result = subprocess.run(
#         cmd,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE,
#         text=True,
#     )

#     # If R failed, raise a clear error instead of crashing later
#     if result.returncode != 0:
#         raise RuntimeError(
#             f"SeqFF failed for {sample_path} (exit code {result.returncode}).\n"
#             f"Command: {' '.join(cmd)}\n"
#             f"STDERR:\n{result.stderr}"
#         )

#     # R prints a named numeric vector like:
#     #   seqff     Enet     WRSC 
#     #  12.3456  11.2345  13.4567
#     lines = [ln for ln in result.stdout.splitlines() if ln.strip()]

#     if len(lines) < 2:
#         raise RuntimeError(
#             f"Unexpected SeqFF output for {sample_path}:\n{result.stdout}"
#         )

#     value_line = lines[1]              # second line: "12.3456 11.2345 13.4567"
#     ff_value = value_line.split()[0]   # first number = seqff

#     value_line = lines[1]              # second line: "12.3456 11.2345 13.4567"
#     Enet = value_line.split()[1]       # second number = Enet

#     value_line = lines[1]              # second line: "12.3456 11.2345 13.4567"
#     WRSC = value_line.split()[2]       # third number = WRSC

#     return ff_value, Enet, WRSC

# def write_result(sample, ff_value, Enet, WRSC, output_path):
#     """
#     Writes a single-sample FF result to a TSV file.
#     """
#     df = pd.DataFrame({
#         "sample": [sample],
#         "ff (SeqFF)": [ff_value],
#         "Enet": [Enet],
#         "WRSC": [WRSC]
#     })
#     df.to_csv(output_path, sep='\t', index=False, na_rep='NA')

# def main():
#     parser = argparse.ArgumentParser(description='Run SeqFF for a single BAM file')
#     parser.add_argument('sample', type=str, help='Absolute path to a single BAM file')
#     parser.add_argument('output', type=str, help='Output TSV file')

#     args = parser.parse_args()
    
#     sample_name = os.path.splitext(os.path.basename(args.sample))[0]
    
#     ff_value, Enet, WRSC = execute_seqff(args.sample)
    
#     write_result(sample_name, ff_value, args.output)

# if __name__ == '__main__':
#     main()
#####################################################################################

import subprocess
import argparse
import os
import pandas as pd
import sys

r_exec = 'Rscript'

current_dir = os.path.dirname(os.path.abspath(__file__))
assets_dir = os.path.join(current_dir, "../assets")
r_script = os.path.join(assets_dir, "seqff.r") 
seqffdir = assets_dir

def execute_seqff(sample_path):
    """Call SeqFF Rscript to compute ff for one BAM file."""
    if not os.path.isfile(sample_path):
        raise FileNotFoundError(f"BAM file does not exist: {sample_path}")

    cmd = [r_exec, r_script, sample_path, seqffdir]

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True  # Raises CalledProcessError if returncode != 0
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"SeqFF failed for {sample_path} (exit code {e.returncode}).\n"
            f"Command: {' '.join(cmd)}\nSTDERR:\n{e.stderr}"
        )

    lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise RuntimeError(
            f"Unexpected SeqFF output for {sample_path}:\n{result.stdout}"
        )

    value_line = lines[1].split()
    if len(value_line) < 3:
        raise RuntimeError(
            f"SeqFF output line does not have 3 values:\n{lines[1]}"
        )

    ff_value, Enet, WRSC = value_line[:3]
    return ff_value, Enet, WRSC

def write_result(sample, ff_value, Enet, WRSC, output_path):
    """Writes a single-sample FF result to a TSV file."""
    df = pd.DataFrame({
        "sample": [sample],
        "ff (SeqFF)": [ff_value],
        "Enet": [Enet],
        "WRSC": [WRSC]
    })
    try:
        df.to_csv(output_path, sep='\t', index=False, na_rep='NA')
    except Exception as e:
        raise IOError(f"Failed to write output file {output_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Run SeqFF for a single BAM file')
    parser.add_argument('sample', type=str, help='Absolute path to a single BAM file')
    parser.add_argument('output', type=str, help='Output TSV file')
    args = parser.parse_args()
    
    sample_name = os.path.splitext(os.path.basename(args.sample))[0]
    
    try:
        ff_value, Enet, WRSC = execute_seqff(args.sample)
        write_result(sample_name, ff_value, Enet, WRSC, args.output)
        print(f"SeqFF results written to {args.output}")
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()