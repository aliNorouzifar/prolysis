import os
from datetime import datetime

# Custom log function
def log_command(message, log_file="log.log"):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"{timestamp} - {message}\n"

    # Append log entry to the log file
    with open(log_file, "a") as f:
        f.write(log_entry)

    print(log_entry.strip())  # Optional: print to the console
