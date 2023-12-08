"""
Mixing subprocess and concurrent futures can make it difficult to catch 
exceptions in a timely manner. This small example is a trial on how to do that
properly in a production code
"""
import subprocess
from concurrent.futures import ProcessPoolExecutor


def task(x, break_on=5):
    if x == break_on:
        raise Exception
    return x**5


def run():
    for x_ in range(10):
        subprocess.run(
