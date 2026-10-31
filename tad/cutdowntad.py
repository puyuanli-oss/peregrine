import argparse
import os
import tempfile


def cutdowntad(eqtl_file, heirarchical_tad_file, tad_file, output_file, chia_file=None):
    """Keep TAD pairs supported by eQTL, hierarchical TAD, or enabled ChIA."""
    supported = set()
    sources = [eqtl_file, heirarchical_tad_file]
    if chia_file is not None:
        sources.append(chia_file)

    for source in sources:
        with open(source, "r") as input_file:
            for number, line in enumerate(input_file, 1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 2 or not all(fields[:2]):
                    raise ValueError(f"{source}:{number}: expected enhancer and PANTHER gene")
                if fields[:2] == ["enhancer", "gene"]:
                    continue
                if source == chia_file and (len(fields) != 4 or not all(fields) or fields[3] != "1"):
                    raise ValueError(f"{source}:{number}: expected four fields with ChIA assay 1")
                supported.add((fields[0], fields[1]))

    retained = 0
    output_dir = os.path.dirname(os.path.abspath(output_file))
    os.makedirs(output_dir, exist_ok=True)
    fd, temporary_path = tempfile.mkstemp(prefix=".selectTAD-", dir=output_dir)
    try:
        with os.fdopen(fd, "w") as out, open(tad_file, "r") as tad:
            for number, line in enumerate(tad, 1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 4 or not all(fields):
                    raise ValueError(f"{tad_file}:{number}: expected four nonempty fields")
                enhancer, panthid, tissue, assay = fields
                if (enhancer, panthid) in supported:
                    retained += 1
                    out.write(f"{enhancer}\t{panthid}\t{tissue}\t{assay}\n")
        os.replace(temporary_path, output_file)
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
    print(f"Number of TAD links found in TADintxn, eQTL, or ChIA: {retained}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_path")
    parser.add_argument("--chia-links", help="Explicit linksDBchia from this run; omitted means no ChIA support")
    args = parser.parse_args()
    target_path = os.path.abspath(args.target_path)
    os.makedirs(target_path, exist_ok=True)
    base = os.path.dirname(target_path)
    cutdowntad(
        os.path.join(base, "eQTL", "linksDBnumeqtl"),
        os.path.join(base, "heirarchicalTAD", "PSYCHIClinksDB"),
        os.path.join(target_path, "linksTADtissues"),
        os.path.join(target_path, "selectTAD"),
        chia_file=args.chia_links,
    )


if __name__ == "__main__":
    main()
