"""环状条码拼接 API 服务。"""
from __future__ import annotations

from fastapi import FastAPI
from pydantic import BaseModel, Field, field_validator

from .assembly import assemble

app = FastAPI(
    title="环状条码拼接服务",
    version="1.0.0",
    description=(
        "将同一环状条码的一批等长短序列（每条恰好使用一次、相邻重叠 k-1）"
        "拼接为规范条码；循环移位与整条反向互补视为同一结果。"
    ),
)


class AssembleRequest(BaseModel):
    sequences: list[str] = Field(..., description="6–24 条等长 DNA 序列，长度 3–8，仅含 A/C/G/T")

    @field_validator("sequences")
    @classmethod
    def _validate_and_normalize(cls, value: list[str]) -> list[str]:
        if not 6 <= len(value) <= 24:
            raise ValueError("序列条数须在 6 至 24 之间")
        normalized = []
        for s in value:
            if not isinstance(s, str):
                raise ValueError("每条序列必须为字符串")
            t = s.strip().upper()
            if not 3 <= len(t) <= 8:
                raise ValueError("每条序列长度须在 3 至 8 之间")
            if any(c not in "ACGT" for c in t):
                raise ValueError("序列仅允许包含 A/C/G/T")
            normalized.append(t)
        if len({len(t) for t in normalized}) != 1:
            raise ValueError("所有序列必须等长")
        return normalized


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/assemble")
def assemble_route(req: AssembleRequest) -> dict:
    return assemble(req.sequences)
