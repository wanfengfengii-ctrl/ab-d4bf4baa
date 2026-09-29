"""一次性验收服务：等待 API 就绪 → 代码测试 → 构建检查 → 三类拼接冒烟。

退出码：
  0  全部通过
  2  API 未在限定时间内就绪
  3  代码测试（pytest）失败
  4  构建检查失败
  5  拼接冒烟失败
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import Counter

EXIT_OK = 0
EXIT_NOT_READY = 2
EXIT_TESTS_FAILED = 3
EXIT_BUILD_FAILED = 4
EXIT_SMOKE_FAILED = 5

BASE_URL = os.environ.get("API_BASE_URL", "http://api:8000").rstrip("/")
READY_TIMEOUT = float(os.environ.get("VERIFY_READY_TIMEOUT", "90"))

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app.assembly import canonical_and_alignment, reverse_complement  # noqa: E402


def log(msg: str) -> None:
    print(f"[verify] {msg}", flush=True)


def http(method: str, path: str, payload: dict | None = None, timeout: float = 10.0):
    req = urllib.request.Request(
        BASE_URL + path,
        method=method,
        data=None if payload is None else json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.status, json.loads(resp.read().decode())


def wait_ready() -> bool:
    deadline = time.monotonic() + READY_TIMEOUT
    while time.monotonic() < deadline:
        try:
            status, body = http("GET", "/health", timeout=3.0)
            if status == 200 and body.get("status") == "ok":
                return True
        except (urllib.error.URLError, OSError, ValueError):
            pass
        time.sleep(1.0)
    return False


def run_unit_tests() -> bool:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/", "-q"],
        cwd=os.path.dirname(os.path.abspath(__file__)),
    )
    return proc.returncode == 0


def build_check() -> bool:
    # 源码可编译、应用包可导入、服务能给出 OpenAPI 文档
    if subprocess.run([sys.executable, "-m", "compileall", "-q", "app"]).returncode != 0:
        log("compileall 失败")
        return False
    try:
        import app.main  # noqa: F401
    except Exception as exc:  # pragma: no cover
        log(f"应用包导入失败: {exc}")
        return False
    try:
        status, _ = http("GET", "/openapi.json")
    except Exception as exc:
        log(f"openapi.json 不可达: {exc}")
        return False
    if status != 200:
        log(f"openapi.json 返回 {status}")
        return False
    return True


def check_witness(witness: dict, sequences: list[str], label: str) -> list[str]:
    """独立复核一份见证的自洽性，返回错误列表（空表示通过）。"""
    errors = []
    n, k = len(sequences), len(sequences[0])
    usage = witness.get("usage_order", [])
    if sorted(usage) != list(range(n)):
        errors.append(f"{label}: usage_order 不是 0..{n-1} 的排列")
        return errors
    overlaps = witness.get("overlaps", [])
    if len(overlaps) != n:
        errors.append(f"{label}: 重叠证据条数 {len(overlaps)} != {n}")
    for i, ov in enumerate(overlaps):
        a, b = sequences[ov["from_index"]], sequences[ov["to_index"]]
        if ov["from_index"] != usage[i] or ov["to_index"] != usage[(i + 1) % n]:
            errors.append(f"{label}: 第 {i} 条重叠证据与使用次序不符")
        if ov["overlap"] != a[1:] or ov["overlap"] != b[:-1] or ov["overlap_length"] != k - 1:
            errors.append(f"{label}: 第 {i} 条重叠证据重叠区不正确")
    assembled = witness.get("assembled", "")
    kmers = [(assembled[i:] + assembled[:i])[:k] for i in range(n)]
    if Counter(kmers) != Counter(sequences):
        errors.append(f"{label}: 拼接串的 k-mer 多重集与输入不一致")
    canonical, strand, _ = canonical_and_alignment(assembled)
    if witness.get("barcode") != canonical or witness.get("strand") != strand:
        errors.append(f"{label}: barcode/strand 与拼接串的规范代表不符")
    expect = assembled if strand == "+" else reverse_complement(assembled)
    if witness.get("barcode") != canonical_and_alignment(expect)[0]:
        errors.append(f"{label}: barcode 与 strand 声明的方向关系不成立")
    return errors


def smoke_unique() -> list[str]:
    seqs = ["AAA", "AAC", "ACC", "CCC", "CCA", "CAA"]
    _, body = http("POST", "/assemble", {"sequences": seqs})
    errors = []
    if body.get("status") != "unique":
        return [f"unique 用例状态异常: {body.get('status')}"]
    if body.get("barcode") != "AAACCC":
        errors.append(f"unique 用例条码异常: {body.get('barcode')}")
    errors += check_witness(body, seqs, "unique")
    # 循环移位与反向互补必须给出同一规范条码
    rc = [reverse_complement(s) for s in seqs]
    _, body2 = http("POST", "/assemble", {"sequences": rc})
    if body2.get("barcode") != body.get("barcode"):
        errors.append("反向互补输入未给出同一规范条码")
    return errors


def smoke_ambiguous() -> list[str]:
    seqs = ["ACT", "CTG", "ACA", "CAT", "ATG", "TGA", "GAC", "TGT", "GTA", "TAC"]
    _, body = http("POST", "/assemble", {"sequences": seqs})
    if body.get("status") != "ambiguous":
        return [f"ambiguous 用例状态异常: {body.get('status')}"]
    errors = []
    witnesses = body.get("witnesses", [])
    if len(witnesses) != 2:
        errors.append(f"见证条数 {len(witnesses)} != 2")
    if len(witnesses) == 2 and witnesses[0].get("barcode") == witnesses[1].get("barcode"):
        errors.append("两条见证的规范条码相同，不构成歧义证据")
    for i, w in enumerate(witnesses):
        errors += check_witness(w, seqs, f"witness[{i}]")
    return errors


def smoke_no_solution() -> list[str]:
    errors = []
    _, body = http("POST", "/assemble", {"sequences": ["ACG"] * 6})
    if body.get("status") != "no_solution" or body.get("reason") != "degree_imbalance":
        errors.append(f"度数失衡用例异常: {body.get('status')}/{body.get('reason')}")
    elif not body.get("details", {}).get("imbalanced_fragments"):
        errors.append("度数失衡用例缺少可复核的失衡片段明细")
    _, body = http("POST", "/assemble", {"sequences": ["AAA"] * 3 + ["CCC"] * 3})
    if body.get("status") != "no_solution" or body.get("reason") != "disconnected":
        errors.append(f"不连通用例异常: {body.get('status')}/{body.get('reason')}")
    elif len(body.get("details", {}).get("components", [])) < 2:
        errors.append("不连通用例缺少可复核的连通分量明细")
    return errors


def main() -> int:
    log(f"等待 API 就绪: {BASE_URL}")
    if not wait_ready():
        log("API 未就绪，退出")
        return EXIT_NOT_READY
    log("API 已就绪")

    log("运行代码测试（pytest）…")
    if not run_unit_tests():
        log("代码测试失败")
        return EXIT_TESTS_FAILED
    log("代码测试通过")

    log("运行构建检查…")
    if not build_check():
        log("构建检查失败")
        return EXIT_BUILD_FAILED
    log("构建检查通过")

    log("运行拼接冒烟：唯一 / 歧义 / 无解 …")
    failures = smoke_unique() + smoke_ambiguous() + smoke_no_solution()
    if failures:
        for f in failures:
            log(f"冒烟失败: {f}")
        return EXIT_SMOKE_FAILED
    log("三类拼接冒烟全部通过")
    log("验收通过")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
