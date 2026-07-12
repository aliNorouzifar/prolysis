"""Shared key/value state for prolysis, with an automatic in-memory fallback.

prolysis stashes small pieces of intermediate analysis state (JSON strings, keyed
by name) so that separate steps — and, in the web tools, separate Dash background
callbacks — can hand results to one another. When a Redis server is reachable
(configured via ``REDIS_HOST`` / ``REDIS_PORT``, default ``localhost:6379``) that
state lives in Redis and can be shared across processes.

When Redis is not installed or not reachable, prolysis transparently falls back to a
process-local :class:`InMemoryStore` implementing the same ``get``/``set``/``ping``
surface. This keeps the package importable and usable for single-process runs, tests,
and demos without requiring a running Redis server.

``get``, ``set`` and ``ping`` are used by prolysis; the web tools additionally call
``flushall`` on the shared ``redis_client``. The fallback implements exactly that
surface. Import ``redis_client`` from this module — it is either a live
``redis.StrictRedis`` client or an :class:`InMemoryStore`, and both expose the same API.
"""
from __future__ import annotations

import os
from typing import Any, Optional


class InMemoryStore:
    """Process-local stand-in for a Redis client (``get``/``set``/``ping``/``flushall``).

    Mimics ``redis.StrictRedis(decode_responses=True)`` closely enough for prolysis:
    values are coerced to ``str`` on write and missing keys return ``None`` on read.
    State is not shared across processes and is lost when the process exits.
    """

    def __init__(self) -> None:
        self._store: dict[str, str] = {}

    def set(self, key: Any, value: Any, *args: Any, **kwargs: Any) -> bool:
        if isinstance(value, (bytes, bytearray)):
            value = value.decode()
        elif not isinstance(value, str):
            value = str(value)
        self._store[str(key)] = value
        return True

    def get(self, key: Any) -> Optional[str]:
        return self._store.get(str(key))

    def ping(self) -> bool:
        return True

    def flushall(self, *args: Any, **kwargs: Any) -> bool:
        self._store.clear()
        return True


def _connect() -> Any:
    """Return a live Redis client, or an :class:`InMemoryStore` if unavailable."""
    try:
        import redis
    except ImportError:
        print("prolysis: 'redis' not installed - using in-memory state fallback.")
        return InMemoryStore()

    client = redis.StrictRedis(
        host=os.getenv("REDIS_HOST", "localhost"),
        port=int(os.getenv("REDIS_PORT", "6379")),
        decode_responses=True,
    )
    try:
        client.ping()
    except Exception as exc:  # redis.ConnectionError and related transport errors
        print(f"prolysis: Redis unreachable ({exc}) - using in-memory state fallback.")
        return InMemoryStore()
    print("prolysis: connected to Redis.")
    return client


# Either a live redis.StrictRedis or an InMemoryStore; both expose get/set/ping.
redis_client = _connect()
