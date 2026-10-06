from math import log10, floor

def find_exp(number) -> int:
    base10 = log10(abs(number))
    return abs(floor(base10))