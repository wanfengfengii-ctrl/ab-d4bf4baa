"""环状 DNA 条码拼接核心逻辑。

模型：把每条输入序列（等长 k-mer）看作 de Bruijn 多重图中的一条有向边，
前/后缀（长度 k-1）为顶点。一次合法拼接 = 恰好使用每条边一次且相邻重叠
k-1 的闭环 ⟺ 图中的一条欧拉回路。回路经过各边的首字符即组成长度为 n 的
环状条码。

结果等价：循环移位与整条反向互补视为同一结果（双链环状分子）；重复序列
互换只影响输入序号的分配，不影响条码本身，因此按字符串等价类判重，
天然不会把重复片段的排列差异误报为不同条码。
"""
from __future__ import annotations

from collections import defaultdict

DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    """整条反向互补。"""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def min_rotation(s: str) -> str:
    """字典序最小的循环移位（条码长度 ≤ 24，直接枚举即可）。"""
    return min(s[i:] + s[:i] for i in range(len(s)))


def canonical_and_alignment(circular: str) -> tuple[str, str, str]:
    """计算环状串的规范代表及与正向拼接串的对齐关系。

    返回 (canonical, strand, aligned)：
    - canonical：所有循环移位与整条反向互补的所有循环移位中字典序最小者，
      即“循环移位 + 反向互补”等价类的唯一规范代表；
    - strand：'+' 表示 canonical 取自正向拼接串，'-' 表示取自反向互补；
    - aligned：circular 的某个循环移位，与 canonical 位置对齐，满足
      aligned == canonical（'+'）或 reverse_complement(aligned) == canonical（'-'）。
    """
    fwd = min_rotation(circular)
    rev = min_rotation(reverse_complement(circular))
    if fwd <= rev:
        return fwd, "+", fwd
    return rev, "-", reverse_complement(rev)


def _rotation_offset(circular: str, target: str) -> int:
    """target 作为 circular 循环移位的最小偏移。"""
    n = len(circular)
    for i in range(n):
        if circular[i:] + circular[:i] == target:
            return i
    raise ValueError("target 不是 circular 的循环移位")  # 不应发生


def _build_witness(labels: list[str], label_indices: dict[str, list[int]]) -> dict:
    """由一条欧拉回路（边标签序列）构建一份可复核的拼接见证。

    - 使用次序：重复序列按输入序号升序稳定分配（第 j 次使用对应该序列
      第 j 小的输入序号，序号从 0 开始）；
    - 相邻重叠证据：包含最后一回到第一的环绕边，重叠长度恒为 k-1。
    """
    n = len(labels)
    circular = "".join(label[0] for label in labels)
    barcode, strand, aligned = canonical_and_alignment(circular)
    offset = _rotation_offset(circular, aligned)
    rotated = labels[offset:] + labels[:offset]

    pools = {label: list(idxs) for label, idxs in label_indices.items()}
    usage_order = [pools[label].pop(0) for label in rotated]

    overlaps = []
    for i in range(n):
        a = rotated[i]
        b = rotated[(i + 1) % n]
        overlaps.append(
            {
                "from_index": usage_order[i],
                "to_index": usage_order[(i + 1) % n],
                "from_sequence": a,
                "to_sequence": b,
                "overlap": a[1:],
                "overlap_length": len(a) - 1,
            }
        )

    return {
        "barcode": barcode,
        "strand": strand,
        "assembled": aligned,
        "circular_length": n,
        "usage_order": usage_order,
        "overlaps": overlaps,
    }


def _enumerate_circuits(adj: dict[str, list[list]], start: str, total: int, limit: int = 2):
    """枚举欧拉回路产生的不同规范等价类，取到 limit 个即停。

    在同一顶点按边标签分组展开（平行重边不重复展开），因此重复序列互换
    不会被重复计数；每走一步检查剩余边（忽略零度顶点）是否仍弱连通——
    这是能从当前顶点完成欧拉路径的必要条件，可安全剪枝。
    """
    radj: dict[str, list[list]] = defaultdict(list)
    for edges in adj.values():
        for edge in edges:
            radj[edge[1]].append(edge)

    def remaining_connected() -> bool:
        nonzero = set()
        for u, edges in adj.items():
            for label, v, cnt in edges:
                if cnt:
                    nonzero.add(u)
                    nonzero.add(v)
        if not nonzero:
            return True
        seed = next(iter(nonzero))
        seen = {seed}
        stack = [seed]
        while stack:
            x = stack.pop()
            for edge in adj.get(x, ()):
                if edge[2] and edge[1] not in seen:
                    seen.add(edge[1])
                    stack.append(edge[1])
            for edge in radj.get(x, ()):
                u = edge[0][:-1]
                if edge[2] and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen == nonzero

    found: list[list[str]] = []
    seen_canonical: set[str] = set()
    path: list[str] = []

    def dfs(v: str) -> None:
        if len(found) >= limit:
            return
        if len(path) == total:
            circular = "".join(label[0] for label in path)
            canonical, _, _ = canonical_and_alignment(circular)
            if canonical not in seen_canonical:
                seen_canonical.add(canonical)
                found.append(list(path))
            return
        for edge in adj[v]:
            if edge[2] == 0:
                continue
            edge[2] -= 1
            path.append(edge[0])
            if remaining_connected():
                dfs(edge[1])
            path.pop()
            edge[2] += 1

    dfs(start)
    return found


def assemble(sequences: list[str]) -> dict:
    """拼接主入口。sequences 已校验：6–24 条、等长、长度 3–8、仅 ACGT。"""
    n = len(sequences)
    k = len(sequences[0])
    base = {"sequence_count": n, "kmer_length": k}

    label_indices: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(sequences):
        label_indices[s].append(i)

    indeg: dict[str, int] = defaultdict(int)
    outdeg: dict[str, int] = defaultdict(int)
    undirected: dict[str, set] = defaultdict(set)
    for label, idxs in label_indices.items():
        u, v = label[:-1], label[1:]
        outdeg[u] += len(idxs)
        indeg[v] += len(idxs)
        undirected[u].add(v)
        undirected[v].add(u)

    nodes = set(outdeg) | set(indeg)

    # 1) 度数平衡：任一 (k-1)-mer 的入度必须等于出度
    imbalanced = [
        {"fragment": x, "in_degree": indeg[x], "out_degree": outdeg[x]}
        for x in sorted(nodes)
        if indeg[x] != outdeg[x]
    ]
    if imbalanced:
        return {
            "status": "no_solution",
            "reason": "degree_imbalance",
            **base,
            "details": {"imbalanced_fragments": imbalanced},
            "message": "无可行闭环：以下 (k-1)-mer 片段的入度与出度不平衡，"
            "任何闭合拼接都无法经过它们恰好守恒。",
        }

    # 2) 连通性：所有非零度片段必须落在同一弱连通分量
    start = min(nodes)
    seen = {start}
    stack = [start]
    while stack:
        x = stack.pop()
        for y in undirected[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    if len(seen) != len(nodes):
        remaining = nodes - seen
        components = [sorted(seen)]
        while remaining:
            seed = min(remaining)
            comp = {seed}
            stack = [seed]
            while stack:
                x = stack.pop()
                for y in undirected[x]:
                    if y not in comp:
                        comp.add(y)
                        stack.append(y)
            components.append(sorted(comp))
            remaining -= comp
        return {
            "status": "no_solution",
            "reason": "disconnected",
            **base,
            "details": {"components": components},
            "message": "无可行闭环：非零度 (k-1)-mer 片段不连通，"
            "各连通分量之间没有任何重叠证据可以衔接。",
        }

    # 3) 枚举欧拉回路对应的规范等价类（最多取 2 个见证）
    adj: dict[str, list[list]] = defaultdict(list)
    for label, idxs in label_indices.items():
        adj[label[:-1]].append([label, label[1:], len(idxs)])
    for edges in adj.values():
        edges.sort(key=lambda e: e[0])

    circuits = _enumerate_circuits(adj, start, n, limit=2)
    witnesses = [_build_witness(labels, label_indices) for labels in circuits]

    if len(witnesses) == 1:
        return {
            "status": "unique",
            **base,
            **witnesses[0],
            "message": "仅存在一个规范等价类：所有合法拼接在循环移位与"
            "反向互补意义下给出同一条码。",
        }
    return {
        "status": "ambiguous",
        **base,
        "witnesses": witnesses,
        "message": "存在多个规范等价类：以下两条规范见证互不为循环移位或"
        "反向互补，仅凭本批序列无法判定唯一条码。",
    }
