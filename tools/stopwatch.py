"""
Command line stop watch
"""
import os
import time

if __name__ == "__main__":
    input("Enter to start")
    start = time.time()
    seconds = 0
    while True:
        end = time.time()
        elapsed = end - start
        if elapsed > seconds:
            minutes = elapsed // 60
            seconds = int(elapsed % 60)
            hours = int(minutes // 60)
            minutes = int(minutes % 60)
            os.system("clear")
            print(f"{hours:0>2}:{minutes:0>2}:{seconds:0>2}")
            seconds += 1
