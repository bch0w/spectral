#!/usr/bin/env python3
import sys

if __name__ == "__main__":
    x = int(sys.argv[1])
    if x == 5:
        print("BREAK")
        sys.exit(-1)
    print(x**5)
    
