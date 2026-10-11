"""Utils components of the 3DPrintase pipeline."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import urllib.request
from .config import VERSION


class PipelineError(RuntimeError):
    pass


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def request(url, *, data=None, timeout=60):
    headers = {"User-Agent": f"3DPrintase/{VERSION}", "Accept": "application/json"}
    if data is not None:
        data = json.dumps(data).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            return response.read()
    except (urllib.error.URLError, TimeoutError) as exc:
        raise PipelineError(f"Network request failed: {url}: {exc}") from exc


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()
