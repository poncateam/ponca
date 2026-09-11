"""
    This files serves as a proxy to specify parameters to python unittest classes.

    For now, it only supports a single arguments . 
"""

import os
import sys
import argparse
import subprocess

if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument("--env-name", required=True)
    parser.add_argument("--env-value", required=True)
    parser.add_argument("cmd", nargs=argparse.REMAINDER) 
    args = parser.parse_args()

    env = os.environ.copy()
    env[args.env_name] = args.env_value

    result = subprocess.run([sys.executable, *args.cmd], env=env)
    sys.exit(result.returncode)