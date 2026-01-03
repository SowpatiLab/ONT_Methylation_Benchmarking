import pysam
import numpy as np
from tqdm import tqdm
import multiprocessing
import os
import tempfile
import shutil
import argparse

# Precompute error probabilities
ERR_TAB = np.power(10.0, -np.arange(0, 129) / 10.0)

def ave_qual(quals, qround=False, tab=ERR_TAB):
    """Calculate average basecall quality (vectorized)."""
    if quals:
        quals = np.fromiter(quals, dtype=np.uint8)
        err_probs = tab[quals]
        mean_err = np.mean(err_probs)
        mq = -10 * np.log10(mean_err)
        return round(mq) if qround else mq
    else:
        return None

def process_chunk(args):
    input_bam, output_bam, region, quality_threshold, length_threshold = args
    total_reads = 0
    passed_reads = 0
    failed_read_ids = []
    passed_quals = []
    passed_lens = []
    pre_filter_qualities = []
    pre_filter_lengths = []

    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:
            for read in bam_in.fetch(region=region):
                if read.is_secondary or read.is_supplementary:
                    continue
                total_reads += 1

                if read.query_qualities and len(read.query_qualities) > 0:
                    avg_quality = ave_qual(read.query_qualities)
                    read_length = read.query_length

                    pre_filter_qualities.append(avg_quality)
                    pre_filter_lengths.append(read_length)

                    if avg_quality >= quality_threshold and read_length >= length_threshold:
                        bam_out.write(read)
                        passed_reads += 1
                        passed_quals.append(avg_quality)
                        passed_lens.append(read_length)
                    else:
                        failed_read_ids.append(read.query_name)
                else:
                    failed_read_ids.append(read.query_name)

    return (output_bam, total_reads, passed_reads, failed_read_ids, passed_quals, passed_lens, pre_filter_qualities, pre_filter_lengths)

def filter_bam_by_average_quality_parallel(input_bam, output_bam, quality_threshold, n_processes=4, length_threshold=0, failed_output_path=None):
    temp_dir = tempfile.mkdtemp(prefix="bam_filter_")
    temp_files = []

    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        references = bam_in.references

    regions = [f"{ref}" for ref in references]

    print(f"Splitting work into {len(regions)} regions across {n_processes} processes...")

    args_list = []
    for idx, region in enumerate(regions):
        temp_out_bam = os.path.join(temp_dir, f"temp_part_{idx}.bam")
        temp_files.append(temp_out_bam)
        args_list.append((input_bam, temp_out_bam, region, quality_threshold, length_threshold))

    total_reads = 0
    passed_reads = 0
    all_failed_read_ids = []
    all_passed_quals = []
    all_passed_lens = []
    all_pre_filter_quals = []
    all_pre_filter_lens = []

    with multiprocessing.Pool(processes=n_processes) as pool:
        for result in tqdm(pool.imap_unordered(process_chunk, args_list), total=len(args_list), desc="Processing chunks"):
            out_file, reads_in_chunk, reads_passed_in_chunk, failed_ids, passed_quals, passed_lens, chunk_quals, chunk_lens = result
            total_reads += reads_in_chunk
            passed_reads += reads_passed_in_chunk
            all_failed_read_ids.extend(failed_ids)
            all_passed_quals.extend(passed_quals)
            all_passed_lens.extend(passed_lens)
            all_pre_filter_quals.extend(chunk_quals)
            all_pre_filter_lens.extend(chunk_lens)

    print("Merging filtered parts...")
    pysam.merge("-f", output_bam, *temp_files)
    pysam.index(output_bam)

    if failed_output_path:
        with open(failed_output_path, 'w') as f:
            for read_id in all_failed_read_ids:
                f.write(f"{read_id}\n")

    shutil.rmtree(temp_dir)

    print(f"\nTotal reads processed: {total_reads}")
    print(f"Reads passed quality and length filter: {passed_reads}")
    print(f"Reads failed quality/length filter: {len(all_failed_read_ids)}")
    print(f"\n--- Before Filtering ---")
    print(f"Minimum average quality seen: {min(all_pre_filter_quals):.2f}")
    print(f"Minimum read length seen: {min(all_pre_filter_lens)}")
    print(f"Median average quality seen: {np.median(all_pre_filter_quals):.2f}")
    print(f"Median read length seen: {int(np.median(all_pre_filter_lens))}")
    print(f"\n--- After Filtering ---")
    print(f"Minimum average quality (passed reads): {min(all_passed_quals):.2f}")
    print(f"Minimum read length (passed reads): {min(all_passed_lens)}")
    print(f"Median average quality (passed reads): {np.median(all_passed_quals):.2f}")
    print(f"Median read length (passed reads): {int(np.median(all_passed_lens))}")

    print(f"\nFiltered BAM created at: {output_bam}")
    if failed_output_path:
        print(f"Failed read IDs written to: {failed_output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a BAM file based on average read quality and read length.")
    parser.add_argument("-i", "--input_bam", required=True, help="Path to input BAM file")
    parser.add_argument("-o", "--output_bam", required=True, help="Path to output BAM file")
    parser.add_argument("-q", "--quality_threshold", type=float, default=10, help="Minimum average base quality (default: 10)")
    parser.add_argument("-l", "--length_threshold", type=int, default=0, help="Minimum read length (default: 0)")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of parallel processes (default: 4)")
    parser.add_argument("-f", "--failed_output", type=str, help="Optional path to output failed read IDs")

    args = parser.parse_args()

    filter_bam_by_average_quality_parallel(
        args.input_bam, 
        args.output_bam, 
        args.quality_threshold, 
        args.processes,
        args.length_threshold,
        args.failed_output
    )