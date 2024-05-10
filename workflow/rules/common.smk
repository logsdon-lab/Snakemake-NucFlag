from collections import defaultdict
from typing import Iterator, DefaultDict


def read_fofn_file(path: str) -> Iterator[str]:
    with open(path, "rt") as fh:
        for line in fh.readlines():
            yield line.strip()


def get_dir_files(dirname: str, ext: str) -> Iterator[str]:
    escaped_ext = re.escape(f".{ext}")
    path_pattern = re.compile(r"([^/]+)(" + escaped_ext + ")")
    for root, read_dirs, fnames in os.walk(dirname):
        for file in fnames:
            read_dir_path = os.path.join(root, file)
            try:
                flowcell_id, _ = re.search(path_pattern, file).groups()
            except (ValueError, AttributeError):
                continue

            yield read_dir_path


def get_sample_assemblies_and_reads() -> tuple[DefaultDict[str, list[str]], DefaultDict[str, list[str]]]:
    SAMPLE_ASSEMBLIES = defaultdict(list)
    SAMPLE_READS = defaultdict(list)

    for sm in config["samples"]:
        sm_name = sm["name"]
        asm_fofn = sm.get("asm_fofn")
        read_fofn = sm.get("read_fofn")

        if asm_fofn:
            for file in read_fofn_file(asm_fofn):
                SAMPLE_ASSEMBLIES[sm_name].append(file)
        elif sm.get("asm_dir") and sm.get("asm_ext"):
            for file in get_dir_files(sm["asm_dir"], sm["asm_ext"]):
                SAMPLE_ASSEMBLIES[sm_name].append(file)
        else:
            raise ValueError("Must provide either asm_fofn or asm_dir and asm_ext.")

        if read_fofn:
            for file in read_fofn_file(read_fofn):
                SAMPLE_READS[sm_name].append(file)
        elif sm.get("read_dir") and sm.get("read_ext"):
            for file in get_dir_files(sm["read_dir"], sm["read_ext"]):
                SAMPLE_READS[sm_name].append(file)
        else:
            raise ValueError("Must provide either read_fofn or read_dir and read_ext.")

    return SAMPLE_ASSEMBLIES, SAMPLE_READS
