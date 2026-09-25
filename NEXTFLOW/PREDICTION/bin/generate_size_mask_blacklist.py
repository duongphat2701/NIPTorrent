#!/usr/bin/env python3
"""
TMB Size Mask Blacklist Generator
==================================
Generate size mask blacklist bed files from input bin files based on z-score thresholds.

This script processes bed files with z-score columns and creates blacklist files by:
1. Identifying regions with absolute z-score > 5
2. Filtering for 1-2 consecutive bins only (1-2 Mb)
3. Excluding X chromosome regions
4. Combining with predefined basic blacklist
5. Sorting chromosomes in numerical order (1,2,3,...,22,X)

Usage:
    python generate_size_mask_blacklist.py --input-folder <path> --output-folder <path> [options]

Example:
    python generate_size_mask_blacklist.py \
        --input-folder data/input_bins \
        --output-folder data/size_mask_regions \
        --basic-blacklist data/basic_blacklist/basic_blacklist.bed \
        --zscore-threshold 5.0 \
        --max-consecutive 2

Author: TMB Pipeline
Date: October 2025
"""

import pandas as pd
import os
import glob
import argparse
import sys
from pathlib import Path


def filter_zscore_regions(df, zscore_threshold=5.0, max_consecutive=2):
    """
    Filter regions based on z-score threshold and consecutive bins constraint.
    Excludes X chromosome regions as per requirements.
    
    Parameters:
    - df: DataFrame with columns ['chr', 'start', 'end', 'id', 'ratio', 'zscore']
    - zscore_threshold: Absolute z-score threshold (default: 5.0)
    - max_consecutive: Maximum consecutive bins allowed (default: 2)
    
    Returns:
    - DataFrame with filtered regions (chr, start, end columns only)
    """
    # Remove rows with NaN z-scores
    df_clean = df.dropna(subset=['zscore']).copy()
    
    if df_clean.empty:
        return pd.DataFrame(columns=['chr', 'start', 'end'])
    
    # Exclude X chromosome regions
    df_clean = df_clean[df_clean['chr'] != 'X'].copy()
    
    if df_clean.empty:
        return pd.DataFrame(columns=['chr', 'start', 'end'])
    
    # Find regions with absolute z-score > threshold
    high_zscore_mask = abs(df_clean['zscore']) > zscore_threshold
    df_clean['high_zscore'] = high_zscore_mask
    
    # Group consecutive regions
    df_clean['group'] = (df_clean['high_zscore'] != df_clean['high_zscore'].shift()).cumsum()
    
    # Filter for valid regions (high z-score groups with <= max_consecutive bins)
    valid_regions = []
    
    for group_id, group_df in df_clean.groupby('group'):
        if group_df['high_zscore'].iloc[0]:  # Only consider high z-score groups
            group_size = len(group_df)
            if 1 <= group_size <= max_consecutive:
                # Add all regions in this valid group
                regions = group_df[['chr', 'start', 'end']].copy()
                valid_regions.append(regions)
    
    if valid_regions:
        result_df = pd.concat(valid_regions, ignore_index=True)
        
        # Create a chromosome sorting key for proper numerical order (1,2,3,...,22,X)
        def chr_sort_key(chr_name):
            """Create sorting key for chromosomes in numerical order"""
            chr_str = str(chr_name)
            if chr_str == 'X':
                return (23, 0)
            elif chr_str == 'Y':
                return (24, 0)
            else:
                try:
                    return (int(chr_str), 0)
                except ValueError:
                    return (25, chr_str)
        
        # Add sorting key column
        result_df['_sort_key'] = result_df['chr'].apply(chr_sort_key)
        
        # Sort by chromosome (numerical order) and start position
        result_df = result_df.sort_values(['_sort_key', 'start']).reset_index(drop=True)
        
        # Remove the sorting key column
        result_df = result_df.drop('_sort_key', axis=1)
        
        return result_df
    else:
        return pd.DataFrame(columns=['chr', 'start', 'end'])


def load_basic_blacklist(blacklist_file):
    """
    Load the predefined basic blacklist bed file.
    
    Parameters:
    - blacklist_file: Path to the basic blacklist bed file
    
    Returns:
    - DataFrame with basic blacklist regions
    """
    try:
        if os.path.exists(blacklist_file):
            df = pd.read_csv(blacklist_file, sep='\t', 
                           names=['chr', 'start', 'end'],
                           header=None)
            # Convert chromosome to string to handle numeric chromosomes
            df['chr'] = df['chr'].astype(str)
            print(f"✓ Loaded basic blacklist: {len(df)} regions from {blacklist_file}")
            return df
        else:
            print(f"⚠ Warning: Basic blacklist file not found: {blacklist_file}")
            return pd.DataFrame(columns=['chr', 'start', 'end'])
    except Exception as e:
        print(f"✗ Error loading basic blacklist: {str(e)}")
        return pd.DataFrame(columns=['chr', 'start', 'end'])


def combine_with_blacklist(filtered_regions, basic_blacklist_df):
    """
    Combine filtered z-score regions with basic blacklist.
    
    Parameters:
    - filtered_regions: DataFrame with filtered z-score regions
    - basic_blacklist_df: DataFrame with basic blacklist regions
    
    Returns:
    - Combined DataFrame sorted by chromosome and start position in numerical order (1,2,3,...,22,X)
    """
    if basic_blacklist_df.empty:
        return filtered_regions
    
    if filtered_regions.empty:
        # Still need to sort basic blacklist
        combined_df = basic_blacklist_df.copy()
    else:
        # Combine both DataFrames
        combined_df = pd.concat([basic_blacklist_df, filtered_regions], ignore_index=True)
    
    # Create a chromosome sorting key for proper numerical order (1,2,3,...,22,X)
    def chr_sort_key(chr_name):
        """Create sorting key for chromosomes in numerical order"""
        chr_str = str(chr_name)
        if chr_str == 'X':
            return (23, 0)
        elif chr_str == 'Y':
            return (24, 0)
        else:
            try:
                return (int(chr_str), 0)
            except ValueError:
                return (25, chr_str)
    
    # Add sorting key column
    combined_df['_sort_key'] = combined_df['chr'].apply(chr_sort_key)
    
    # Sort by chromosome (numerical order) and start position
    combined_df = combined_df.sort_values(['_sort_key', 'start']).reset_index(drop=True)
    
    # Remove the sorting key column
    combined_df = combined_df.drop('_sort_key', axis=1)
    
    # Remove duplicates (if any)
    combined_df = combined_df.drop_duplicates().reset_index(drop=True)
    
    return combined_df


def process_single_bed_file(bed_file, zscore_threshold, max_consecutive):
    """
    Process a single bed file and extract valid z-score regions.
    
    Parameters:
    - bed_file: Path to input bed file
    - zscore_threshold: Absolute z-score threshold
    - max_consecutive: Maximum consecutive bins allowed
    
    Returns:
    - DataFrame with filtered regions
    - Sample name
    """
    try:
        # Extract sample name from filename
        filename = os.path.basename(bed_file)
        sample_name = filename.replace("_bins.bed", "")
        
        # Read bed file
        df = pd.read_csv(bed_file, sep='\t', 
                       names=['chr', 'start', 'end', 'id', 'ratio', 'zscore'],
                       skiprows=1)  # Skip header row
        
        # Convert chromosome to string
        df['chr'] = df['chr'].astype(str)
        
        # Filter regions based on z-score
        filtered_df = filter_zscore_regions(df, zscore_threshold, max_consecutive)
        
        return filtered_df, sample_name
        
    except Exception as e:
        print(f"✗ Error processing {bed_file}: {str(e)}")
        return pd.DataFrame(columns=['chr', 'start', 'end']), None


def generate_blacklist_files(input_folder, output_folder, basic_blacklist_file=None,
                            zscore_threshold=5.0, max_consecutive=2, 
                            include_basic_blacklist=True, verbose=True):
    """
    Main function to generate blacklist files from input bed files.
    
    Parameters:
    - input_folder: Path to folder containing input bed files
    - output_folder: Path to output directory for blacklist files
    - basic_blacklist_file: Path to basic blacklist bed file (optional)
    - zscore_threshold: Absolute z-score threshold (default: 5.0)
    - max_consecutive: Maximum consecutive bins allowed (default: 2)
    - include_basic_blacklist: Whether to include predefined blacklist (default: True)
    - verbose: Print detailed progress (default: True)
    
    Returns:
    - Number of successfully generated files
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Find all bed files in the input folder
    bed_pattern = os.path.join(input_folder, "*_bins.bed")
    bed_files = glob.glob(bed_pattern)
    
    if len(bed_files) == 0:
        print(f"✗ Error: No bed files found in {input_folder}")
        print(f"  Expected pattern: *_bins.bed")
        return 0
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"TMB SIZE MASK BLACKLIST GENERATOR")
        print(f"{'='*70}")
        print(f"Configuration:")
        print(f"  • Input folder: {input_folder}")
        print(f"  • Output folder: {output_folder}")
        print(f"  • Z-score threshold: ±{zscore_threshold}")
        print(f"  • Max consecutive bins: {max_consecutive}")
        print(f"  • Include basic blacklist: {include_basic_blacklist}")
        if include_basic_blacklist and basic_blacklist_file:
            print(f"  • Basic blacklist file: {basic_blacklist_file}")
        print(f"\n✓ Found {len(bed_files)} bed files to process")
        print(f"{'='*70}\n")
    
    # Load basic blacklist if enabled
    if include_basic_blacklist and basic_blacklist_file:
        basic_blacklist_df = load_basic_blacklist(basic_blacklist_file)
    else:
        basic_blacklist_df = pd.DataFrame(columns=['chr', 'start', 'end'])
        if verbose:
            print("  (Basic blacklist not included)\n")
    
    # Process each bed file
    successful_count = 0
    failed_count = 0
    
    for i, bed_file in enumerate(bed_files, 1):
        if verbose and i % 50 == 0:
            print(f"  Processing: {i}/{len(bed_files)} files...")
        
        # Process bed file
        filtered_df, sample_name = process_single_bed_file(bed_file, zscore_threshold, max_consecutive)
        
        if sample_name is None:
            failed_count += 1
            continue
        
        # Combine with basic blacklist
        final_blacklist = combine_with_blacklist(filtered_df, basic_blacklist_df)
        
        # Generate output filename
        output_filename = f"{sample_name}_size_mask_blacklist.bed"
        output_path = os.path.join(output_folder, output_filename)
        
        # Save to file (no header, tab-separated)
        try:
            final_blacklist.to_csv(output_path, sep='\t', header=False, index=False)
            successful_count += 1
            
            if verbose and i <= 3:  # Show details for first 3 files
                print(f"  ✓ {output_filename}: {len(final_blacklist)} regions", end="")
                if len(filtered_df) == 0 and not basic_blacklist_df.empty:
                    print(" (basic blacklist only)")
                elif include_basic_blacklist and not basic_blacklist_df.empty:
                    print(f" ({len(basic_blacklist_df)} basic + {len(filtered_df)} z-score)")
                else:
                    print(f" ({len(filtered_df)} z-score only)")
                    
        except Exception as e:
            print(f"  ✗ Error writing {output_filename}: {str(e)}")
            failed_count += 1
    
    # Summary
    if verbose:
        print(f"\n{'='*70}")
        print(f"PROCESSING COMPLETE")
        print(f"{'='*70}")
        print(f"✓ Successfully generated: {successful_count} files")
        if failed_count > 0:
            print(f"✗ Failed: {failed_count} files")
        print(f"📁 Output location: {output_folder}")
        print(f"{'='*70}\n")
    
    return successful_count


def main():
    """Main entry point for command-line usage"""
    
    parser = argparse.ArgumentParser(
        description='Generate TMB size mask blacklist files from bin bed files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python %(prog)s --input-folder data/input_bins --output-folder data/size_mask_regions

  # With custom basic blacklist
  python %(prog)s --input-folder data/input_bins \\
                  --output-folder data/size_mask_regions \\
                  --basic-blacklist data/basic_blacklist/basic_blacklist.bed

  # Custom thresholds
  python %(prog)s --input-folder data/input_bins \\
                  --output-folder data/size_mask_regions \\
                  --zscore-threshold 6.0 \\
                  --max-consecutive 3

  # Without basic blacklist
  python %(prog)s --input-folder data/input_bins \\
                  --output-folder data/size_mask_regions \\
                  --no-basic-blacklist
        """
    )
    
    # Required arguments
    parser.add_argument('--input-folder', '-i', required=True,
                       help='Path to folder containing input bin bed files (*_bins.bed)')
    
    parser.add_argument('--output-folder', '-o', required=True,
                       help='Path to output folder for blacklist files')
    
    # Optional arguments
    parser.add_argument('--basic-blacklist', '-b', default=None,
                       help='Path to basic blacklist bed file (default: None)')
    
    parser.add_argument('--zscore-threshold', '-z', type=float, default=5.0,
                       help='Absolute z-score threshold (default: 5.0)')
    
    parser.add_argument('--max-consecutive', '-m', type=int, default=2,
                       help='Maximum consecutive bins allowed (default: 2)')
    
    parser.add_argument('--no-basic-blacklist', action='store_true',
                       help='Do not include basic blacklist in output')
    
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress progress output')
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Validate input folder
    if not os.path.exists(args.input_folder):
        print(f"✗ Error: Input folder does not exist: {args.input_folder}")
        sys.exit(1)
    
    # Validate basic blacklist file if provided
    if args.basic_blacklist and not os.path.exists(args.basic_blacklist):
        print(f"✗ Error: Basic blacklist file does not exist: {args.basic_blacklist}")
        sys.exit(1)
    
    # Run the generator
    include_blacklist = not args.no_basic_blacklist
    
    num_generated = generate_blacklist_files(
        input_folder=args.input_folder,
        output_folder=args.output_folder,
        basic_blacklist_file=args.basic_blacklist,
        zscore_threshold=args.zscore_threshold,
        max_consecutive=args.max_consecutive,
        include_basic_blacklist=include_blacklist,
        verbose=not args.quiet
    )
    
    # Exit with appropriate status code
    if num_generated > 0:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
