import os
import subprocess
import sys
import time
from pathlib import Path

# Paths to directories
example_dirs = ["example_files1", "example_files2"]

# Path to the mapper script (can be overridden: python3 run_tests.py other_mapper.py)
mapper_script = "/Users/alicja/Documents/26:27Z/algo/project/debug_memory/mapper_memory.py"
if len(sys.argv) > 1:
    mapper_script = sys.argv[1]

# Output directory
output_dir = "test_outputs2"
os.makedirs(output_dir, exist_ok=True)

# A read is correctly mapped if both coordinates are within this distance
TOLERANCE = 20


def read_ids(fasta_file):
    """Identifiers of all reads in a fasta file (first word after '>')."""
    ids = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def read_mapping(path):
    """Read a file with lines: read_id <tab> start <tab> end.
    Returns dict read_id -> (start, end). Malformed lines are skipped."""
    mapping = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                mapping[parts[0]] = (int(parts[1]), int(parts[2]))
            except ValueError:
                continue
    return mapping


def compare(reads_file, output_file, report_file, elapsed):
    """Compare mapper output with the answers (same name as the reads file,
    ending .txt, in the same directory). Writes a detailed report and
    returns a one-line summary."""
    answers_file = Path(reads_file).with_suffix(".txt")
    if not answers_file.exists():
        return f"no answers file {answers_file}"

    total = len(read_ids(reads_file))
    answers = read_mapping(answers_file)      # mappable reads and true coordinates
    mapped = read_mapping(output_file)        # what the mapper reported

    correct = 0
    wrong = []            # mapped, mappable, but coordinates too far
    unmappable = []       # mapped, but the read is not in the answers
    for rid, (s, e) in mapped.items():
        if rid not in answers:
            unmappable.append(rid)
            continue
        ts, te = answers[rid]
        if abs(s - ts) <= TOLERANCE and abs(e - te) <= TOLERANCE:
            correct += 1
        else:
            wrong.append((rid, s, e, ts, te))
    missed = [rid for rid in answers if rid not in mapped]
    incorrect = len(wrong) + len(unmappable)

    def pct(a, b):
        return 100.0 * a / b if b else 0.0

    lines = [
        f"reads total:                 {total}",
        f"mappable reads (in answers): {len(answers)}",
        f"reads mapped by mapper:      {len(mapped)}  ({pct(len(mapped), total):.2f}% of all reads)",
        f"correctly mapped:            {correct}  ({pct(correct, len(answers)):.2f}% of mappable)",
        f"incorrectly mapped:          {incorrect}  ({pct(incorrect, len(mapped)):.2f}% of mapped)",
        f"   wrong coordinates:        {len(wrong)}",
        f"   unmappable read mapped:   {len(unmappable)}",
        f"mappable but not mapped:     {len(missed)}",
        f"time:                        {elapsed:.1f} s  ({elapsed / total if total else 0:.3f} s per read)",
        "",
        "== wrong coordinates: read, reported start/end, true start/end ==",
    ]
    for rid, s, e, ts, te in wrong:
        lines.append(f"{rid}\t{s}\t{e}\t{ts}\t{te}")
    lines.append("")
    lines.append("== unmappable reads that were mapped ==")
    lines.extend(unmappable)
    lines.append("")
    lines.append("== mappable reads not mapped ==")
    lines.extend(missed)

    with open(report_file, "w") as f:
        f.write("\n".join(lines) + "\n")

    return (f"mapped {pct(len(mapped), total):.1f}% of all, "
            f"correct {pct(correct, len(answers)):.1f}% of mappable, "
            f"incorrect {incorrect} ({pct(incorrect, len(mapped)):.2f}% of mapped), "
            f"{elapsed / total if total else 0:.3f} s/read")


# Short summary of all runs
summary_file = os.path.join(output_dir, "summary.log")

with open(summary_file, "w") as summary:
    for example_dir in example_dirs:
        # Get all reads*.fasta files
        reads_files = sorted(Path(example_dir).glob("reads*.fasta"))

        # Get reference file(s)
        if example_dir == "example_files1":
            reference_files = [Path(example_dir) / "reference.fasta"]
        else:
            reference_files = sorted(Path(example_dir).glob("reference*.fasta"))

        # Run mapper for each reads file and reference file
        for reads_file in reads_files:
            for reference_file in reference_files:
                # Common name of this run (example dir included to avoid name clashes)
                name = f"{example_dir}_{reads_file.stem}_to_{reference_file.stem}"

                # 1) aligned reads: written by the mapper itself (3rd argument)
                output_file = os.path.join(output_dir, f"{name}.out")
                # 2) console text (memory usage prints, errors): stdout + stderr
                console_file = os.path.join(output_dir, f"{name}.console.log")
                # 3) comparison with the answers
                report_file = os.path.join(output_dir, f"{name}.compare.txt")

                command = [
                    "python3", mapper_script,
                    str(reference_file), str(reads_file), output_file
                ]

                print(f"Running command: {' '.join(command)}")

                start = time.time()
                with open(console_file, "w") as console:
                    result = subprocess.run(
                        command, stdout=console, stderr=subprocess.STDOUT
                    )
                elapsed = time.time() - start

                if result.returncode != 0 or not os.path.exists(output_file):
                    comparison = "mapper failed, see console log"
                else:
                    comparison = compare(reads_file, output_file, report_file, elapsed)

                line = (f"{name}: return code {result.returncode}, "
                        f"{elapsed:.1f} s | {comparison}\n")
                print(line, end="")
                summary.write(line)
                summary.flush()