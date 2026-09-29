"""API 层测试：健康检查、正常拼接与输入校验。"""
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

UNIQUE = ["AAA", "AAC", "ACC", "CCC", "CCA", "CAA"]
AMBIGUOUS = ["ACT", "CTG", "ACA", "CAT", "ATG", "TGA", "GAC", "TGT", "GTA", "TAC"]


def test_health():
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_assemble_unique():
    r = client.post("/assemble", json={"sequences": UNIQUE})
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "unique"
    assert body["barcode"] == "AAACCC"
    assert body["kmer_length"] == 3
    assert body["sequence_count"] == 6


def test_assemble_ambiguous():
    r = client.post("/assemble", json={"sequences": AMBIGUOUS})
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "ambiguous"
    assert len(body["witnesses"]) == 2
    assert body["witnesses"][0]["barcode"] != body["witnesses"][1]["barcode"]


def test_assemble_no_solution():
    r = client.post("/assemble", json={"sequences": ["ACG"] * 6})
    assert r.status_code == 200
    assert r.json()["reason"] == "degree_imbalance"


def test_lowercase_normalized():
    r = client.post("/assemble", json={"sequences": [s.lower() for s in UNIQUE]})
    assert r.status_code == 200
    assert r.json()["barcode"] == "AAACCC"


def test_reject_too_few_sequences():
    assert client.post("/assemble", json={"sequences": UNIQUE[:5]}).status_code == 422


def test_reject_too_many_sequences():
    assert client.post("/assemble", json={"sequences": ["AAA"] * 25}).status_code == 422


def test_reject_bad_length():
    assert client.post("/assemble", json={"sequences": ["AC"] * 6}).status_code == 422
    assert client.post("/assemble", json={"sequences": ["ACGTACGTA"] * 6}).status_code == 422


def test_reject_unequal_lengths():
    seqs = ["AAA", "AAC", "ACC", "CCC", "CCA", "CAAT"]
    assert client.post("/assemble", json={"sequences": seqs}).status_code == 422


def test_reject_invalid_alphabet():
    seqs = ["AAA", "AAC", "ACC", "CCC", "CCA", "CAN"]
    assert client.post("/assemble", json={"sequences": seqs}).status_code == 422
