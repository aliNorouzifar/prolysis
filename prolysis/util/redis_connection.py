
import redis
import os

# Initialize Redis client
redis_client = redis.StrictRedis(
    host=os.getenv('REDIS_HOST', 'localhost'),
    port=6379,
    decode_responses=True
)

# Test connection (optional)
try:
    redis_client.ping()
    print("Connected to Redis!")
except redis.ConnectionError as e:
    print(f"Failed to connect to Redis: {e}")
