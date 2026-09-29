# 环状条码拼接服务

病原监测场景：一批来自同一环状条码的等长短序列（每条视为一份证据），判断能否
拼成可信条码，避免把重复片段的排列差异误报为不同病原株。

## 模型

把每条 k-mer 序列看作 de Bruijn 多重图中的一条有向边（前/后缀 k-1 为顶点）：

- **可行闭环** ⟺ 图中存在欧拉回路：每条序列恰好使用一次，相邻重叠长度恒为 k-1；
- **同一结果**：环状条码的循环移位与整条反向互补（双链分子）视为同一等价类，
  规范代表取该类中字典序最小者；重复序列互换只影响输入序号分配，不产生新条码；
- **判定**：
  - 仅一个规范等价类 → `unique`，返回唯一条码；
  - 多个等价类 → `ambiguous`，返回两条不同的规范见证；
  - 欧拉回路不存在 → `no_solution`，给出可复核原因：
    `degree_imbalance`（(k-1)-mer 入/出度失衡，附失衡片段明细）或
    `disconnected`（非零度片段不连通，附连通分量明细）。

## 运行

```bash
# 构建并启动，verify 一次性验收后自动退出，退出码即验收结果
docker compose up --build --exit-code-from verify

# 自定义端口（默认 8000）
API_PORT=9000 docker compose up --build --exit-code-from verify
```

`verify` 服务流程：等待 API 就绪（`depends_on: service_healthy` + 主动轮询
`/health`）→ 代码测试（pytest）→ 构建检查（源码编译、应用导入、OpenAPI 可达）
→ 唯一 / 歧义 / 无解三类拼接冒烟 → 自行退出。

退出码：`0` 全部通过；`2` API 未就绪；`3` 代码测试失败；`4` 构建检查失败；
`5` 拼接冒烟失败。

## API

服务端口由环境变量 `API_PORT` 控制（Dockerfile 的 `HEALTHCHECK` 同步跟随）。

### `GET /health`

```json
{"status": "ok"}
```

### `POST /assemble`

请求：6–24 条等长 DNA 序列，长度 3–8，仅含 A/C/G/T（大小写不敏感）。

```json
{"sequences": ["AAA", "AAC", "ACC", "CCC", "CCA", "CAA"]}
```

唯一响应（`status: "unique"`）：

```json
{
  "status": "unique",
  "sequence_count": 6,
  "kmer_length": 3,
  "barcode": "AAACCC",
  "strand": "+",
  "assembled": "AAACCC",
  "circular_length": 6,
  "usage_order": [0, 1, 2, 3, 4, 5],
  "overlaps": [
    {"from_index": 0, "to_index": 1, "from_sequence": "AAA",
     "to_sequence": "AAC", "overlap": "AA", "overlap_length": 2}
  ]
}
```

- `barcode`：规范条码（等价类字典序最小代表）；`strand` 为 `-` 时表示
  `barcode == reverse_complement(assembled)`；
- `usage_order`：按输入序号（从 0 开始）的使用次序，重复序列按序号升序稳定分配；
- `overlaps`：相邻重叠证据，含最后一回到第一的环绕边，重叠长度恒为 k-1。

歧义响应（`status: "ambiguous"`）：`witnesses` 为两份上述结构的见证，
规范条码互不相同。无解响应（`status: "no_solution"`）：`reason` +
`details` 给出失衡片段或连通分量明细。输入不合法返回 422。

## 本地开发

```bash
pip install -r requirements.txt
pytest tests/
uvicorn app.main:app --port ${API_PORT:-8000}
```

## 目录

```
app/assembly.py   拼接核心（欧拉回路、等价类枚举、见证构建）
app/main.py       FastAPI 入口与输入校验
tests/            单元与 API 测试
verify.py         一次性验收脚本
Dockerfile        API 镜像（含 HEALTHCHECK，API_PORT 可配置）
docker-compose.yml  api + verify 双服务编排
```
