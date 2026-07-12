"""Tests for the Redis-with-in-memory-fallback state layer.

These are fully environment-independent: they do not require a running Redis server
(nor the ``redis`` package) and exercise the in-memory fallback directly.
"""
from prolysis.util.redis_connection import InMemoryStore, redis_client


def test_inmemory_get_set_ping_roundtrip():
    store = InMemoryStore()
    assert store.ping() is True
    assert store.get("absent") is None
    store.set("k", "v")
    assert store.get("k") == "v"


def test_inmemory_coerces_values_to_str_like_redis():
    store = InMemoryStore()
    store.set("n", 5)
    assert store.get("n") == "5"           # ints -> str, as redis stores them
    store.set("b", b"bytes")
    assert store.get("b") == "bytes"       # bytes decoded under decode_responses
    store.set("d", {"a": 1})
    assert store.get("d") == "{'a': 1}"


def test_module_client_exposes_the_common_surface():
    # redis_client is a live client or an InMemoryStore; both expose get/set/ping.
    for name in ("get", "set", "ping"):
        assert callable(getattr(redis_client, name))
