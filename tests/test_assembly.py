"""拼接核心逻辑的单元测试。"""
from collections import Counter

from app.assembly import (
    assemble,
    canonical_and_alignment,
    min_rotation,
    reverse_complement,
)

UNIQUE = ["AAA", "AAC", "ACC", "CCC", "CCA", "CAA"]
DUPLICATES = ["AAC", "AAA", "ACA", "AAA", "CAA", "AAA"]
AMBIGUOUS = ["ACT", "CTG", "ACA", "CAT", "ATG", "TGA", "GAC", "TGT", "GTA", "TAC"]
IMBALANCED = ["ACG"] * 6
DISCONNECTED = ["AAA"] * 3 + ["CCC"] * 3


def check_witness(witness: dict, sequences: list[str]) -> None:
    """复核一份见证：使用次序、重叠证据与条码三者必须自洽。"""
    n = len(sequences)
    k = len(sequences[0])
    usage = witness["usage_order"]
    assert sorted(usage) == list(range(n)), "每条序列必须恰好使用一次"

    overlaps = witness["overlaps"]
    assert len(overlaps) == n, "重叠证据须覆盖每一对相邻序列（含环绕）"
    for i, ov in enumerate(overlaps):
        a = sequences[ov["from_index"]]
        b = sequences[ov["to_index"]]
        assert ov["from_index"] == usage[i]
        assert ov["to_index"] == usage[(i + 1) % n]
        assert ov["from_sequence"] == a and ov["to_sequence"] == b
        assert ov["overlap"] == a[1:] == b[:-1]
        assert ov["overlap_length"] == k - 1

    assembled = witness["assembled"]
    assert len(assembled) == n
    # 环状串的 k-mer 多重集必须等于输入多重集
    kmers = [(assembled[i:] + assembled[:i])[:k] if i + k > n else assembled[i : i + k] for i in range(n)]
    assert Counter(kmers) == Counter(sequences)
    # 条码与正向拼接串的关系由 strand 决定
    if witness["strand"] == "+":
        assert witness["barcode"] == min_rotation(assembled)
    else:
        assert witness["barcode"] == min_rotation(reverse_complement(assembled))
    # 条码必须是规范代表
    canonical, _, _ = canonical_and_alignment(assembled)
    assert witness["barcode"] == canonical


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAACCC") == "GGGTTT"
    assert reverse_complement("") == ""


def test_min_rotation():
    assert min_rotation("CAAACC") == "AAACCC"
    assert min_rotation("AAAAAA") == "AAAAAA"


def test_unique_barcode():
    r = assemble(UNIQUE)
    assert r["status"] == "unique"
    assert r["barcode"] == "AAACCC"
    assert r["strand"] == "+"
    assert r["usage_order"] == [0, 1, 2, 3, 4, 5]
    check_witness(r, UNIQUE)


def test_duplicate_sequences_stable_assignment():
    r = assemble(DUPLICATES)
    assert r["status"] == "unique"
    assert r["barcode"] == "AAAAAC"
    # 三条 AAA 按输入序号 1,3,5 稳定分配，不得因互换产生别的结果
    assert r["usage_order"] == [1, 3, 5, 0, 2, 4]
    check_witness(r, DUPLICATES)


def test_input_order_does_not_change_barcode():
    rotated = UNIQUE[2:] + UNIQUE[:2]
    assert assemble(rotated)["barcode"] == assemble(UNIQUE)["barcode"]


def test_reverse_complement_reads_give_same_barcode():
    rc = [reverse_complement(s) for s in UNIQUE]
    r = assemble(rc)
    assert r["status"] == "unique"
    assert r["barcode"] == assemble(UNIQUE)["barcode"]
    assert r["strand"] == "-"
    check_witness(r, rc)


def test_ambiguous_two_distinct_witnesses():
    r = assemble(AMBIGUOUS)
    assert r["status"] == "ambiguous"
    witnesses = r["witnesses"]
    assert len(witnesses) == 2
    assert witnesses[0]["barcode"] != witnesses[1]["barcode"]
    for w in witnesses:
        check_witness(w, AMBIGUOUS)


def test_no_solution_degree_imbalance():
    r = assemble(IMBALANCED)
    assert r["status"] == "no_solution"
    assert r["reason"] == "degree_imbalance"
    frags = {f["fragment"]: (f["in_degree"], f["out_degree"]) for f in r["details"]["imbalanced_fragments"]}
    assert frags == {"AC": (0, 6), "CG": (6, 0)}


def test_no_solution_disconnected():
    r = assemble(DISCONNECTED)
    assert r["status"] == "no_solution"
    assert r["reason"] == "disconnected"
    assert r["details"]["components"] == [["AA"], ["CC"]]


def test_longer_kmers_unique():
    # 由环状串 ACGTTCGA（长度 8）取 4-mer，恰有 8 条序列
    s = "ACGTTCGA"
    seqs = [(s[i:] + s[:i])[:4] for i in range(8)]
    r = assemble(seqs)
    assert r["status"] == "unique"
    assert r["kmer_length"] == 4
    check_witness(r, seqs)
