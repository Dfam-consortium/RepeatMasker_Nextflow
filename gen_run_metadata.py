import argparse
import json
import sys
import subprocess
from datetime import datetime
import os
import hashlib


def main(*args):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-a", "--assembly_accession", required=True)
    parser.add_argument("-r", "--rep_base", required=True)
    parser.add_argument("-g", "--algorithm", required=True)
    parser.add_argument("-c", "--rmsk_commands", required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--species")
    group.add_argument("-l", "--lib")

    args = parser.parse_args()

    rep_masker_dir = "/lustre/work/daray/software/RepeatMasker-4.1.2-p1"
    famdb_py = os.path.join(rep_masker_dir, "famdb.py")
    famdb_lib = os.path.join(rep_masker_dir, "Libraries/famdb")
    rep_masker = os.path.join(rep_masker_dir, "RepeatMasker")

    metadata = {
        "assembly": args.assembly_accession,
        "repbase_ver": args.rep_base if args.rep_base else None,
        "run_date": datetime.now().strftime("%Y-%m-%d"),
        "algorithm": args.algorithm,
        "rmsk_commands": args.rmsk_commands,
    }

    if args.species:
        metadata["lib_sequence_count"] = None
        metadata["lib_md5"] = None
        metadata["lib_base_ct"] = None

        # FamDB Names Command
        names_result = subprocess.run(
            [
                famdb_py,
                "-i",
                famdb_lib,
                "names",
                args.species,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout.decode("utf-8")
        lines = names_result.split("\n")
        if "Exact Matches" not in lines:
            print("Error: Exact name match not found")
            exit(1)
        for line in lines:
            if line.startswith("Taxon: "):
                names = line.split(",")
                for name in names:
                    if name.startswith("Taxon:"):
                        metadata["tax_id"] = name.split(" ")[-1]
                    if "(sanitized scientific name)" in name:
                        metadata["species"] = name.strip().split(" ")[0]
                break

        # FamDB Lineage Command
        lineage_result = subprocess.run(
            [
                famdb_py,
                "-i",
                famdb_lib,
                "lineage",
                "-a",
                "--format",
                "totals",
                metadata["species"],
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout.decode("utf-8")
        lines = lineage_result.split(";")
        for line in lines:
            if "entries in ancestors" in line:
                metadata["ancestral_fams"] = int(line.strip().split(" ")[0])
            if "lineage-specific entries" in line:
                metadata["specific_fams"] = int(line.strip().split(" ")[0])

        # FamDB Info Command
        info_result = subprocess.run(
            [
                famdb_py,
                "-i",
                famdb_lib,
                "info",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout.decode("utf-8")
        lines = info_result.split("\n")
        for line in lines:
            if line.startswith("FamDB Creation Format Version : ") or line.startswith(
                "FamDB Format Version:"
            ):
                metadata["famdb_ver"] = line.strip().split(" ")[-1]
            if line.startswith("Version : "):
                metadata["dfam_ver"] = line.strip().split(" ")[-1]

    elif args.lib:
        metadata["tax_id"] = None
        metadata["species"] = None
        metadata["ancestral_fams"] = None
        metadata["specific_fams"] = None
        metadata["famdb_ver"] = None
        metadata["dfam_ver"] = None

        try:
            with open(args.lib, "r") as file:
                libfile = file.read()
            metadata["lib_sequence_count"] = libfile.count(">")

            metadata["lib_md5"] = hashlib.md5(open(args.lib, "rb").read()).hexdigest()

            lines = libfile.split("\n")
            nucleotides = {"A", "T", "C", "G"}
            base_count = 0
            for line in lines:
                sline = line.strip()
                if set(sline).issubset(nucleotides):
                    base_count += len(sline)
            metadata["lib_base_ct"] = base_count

        except Exception as e:
            print(f"Couldn't Open {args.lib}, {e}")

    # RepeatMasker Version Command
    rm_result = subprocess.run(
        [rep_masker, "-v"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout.decode("utf-8")
    metadata["rmsk_ver"] = rm_result.strip().split(" ")[-1]

    # check if the metadata is complete
    element_check = [
        "assembly",
        "tax_id",
        "species",
        "famdb_ver",
        "repbase_ver",
        "run_date",
        "algorithm",
        "rmsk_commands",
        "ancestral_fams",
        "specific_fams",
        "rmsk_ver",
        "lib_sequence_count",
        "lib_md5",
        "lib_base_ct",
    ]
    missing = []
    for element in element_check:
        if element not in metadata:
            missing.append(element)

    if missing:
        print(f"Missing metadata elements: {missing}")
        exit(1)

    metadata["public"] = True
    metadata["path"] = None
    with open(
        f"{metadata['assembly']}_{metadata['algorithm']}-run_data.json", "w"
    ) as file:
        json.dump(metadata, file)


if __name__ == "__main__":
    main(*sys.argv)
