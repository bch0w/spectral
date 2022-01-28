"""
Create a stacked plot of waveforms that spell out a letter or word or series
of words, in the same vein as the Joy Division Untold Pleasure album art, which
itself is a stacked plot of pulsar radio signals.
"""


class SixteenSegment:
    """
    A class to define the look of a sixteen-segment display
    """
    def __init__(a, b, c, d, e, f, g):
        """
        Segments are labelled a-g, starting with a at the top, moving CCW around
        the segment, ending with g in the middle. e.g.
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.g = g

    def ab(self, choice):
        """
        Defines the alphabet in sixteen segment display bits.
        Bit values were defined: https://github.com/dmadison/LED-Segment-ASCII
        """
        return {"A": "1000100011001111",
			    "B": "0010101000111111",
			    "C": "0000000011110011",
			    "D": "0010001000111111",
			    "E": "1000000011110011",
			    "F": "1000000011000011",
			    "G": "0000100011111011",
			    "H": "1000100011001100",
			    "I": "0010001000110011",
			    "J": "0000000001111100",
			    "K": "1001010011000000",
			    "L": "0000000011110000",
			    "M": "0000010111001100",
			    "N": "0001000111001100",
			    "O": "0000000011111111",
			    "P": "1000100011000111",
			    "Q": "0001000011111111",
			    "R": "1001100011000111",
			    "S": "1000100010111011",
			    "T": "0010001000000011",
			    "U": "0000000011111100",
			    "V": "0100010011000000",
			    "W": "0101000011001100",
			    "X": "0101010100000000",
			    "Y": "1000100010111100",
			    "Z": "0100010000110011",
			    "a": "1010000001110000",
			    "b": "1010000011100000",
			    "c": "1000000001100000",
			    "d": "0010100000011100",
			    "e": "1100000001100000",
			    "f": "1010101000000010",
			    "g": "1010001010100001",
			    "h": "1010000011000000",
			    "i": "0010000000000000",
			    "j": "0010001001100000",
			    "k": "0011011000000000",
			    "l": "0000000011000000",
			    "m": "1010100001001000",
			    "n": "1010000001000000",
			    "o": "1010000001100000",
			    "p": "1000001011000001",
			    "q": "1010001010000001",
			    "r": "1000000001000000",
			    "s": "1010000010100001",
			    "t": "1000000011100000",
			    "u": "0010000001100000",
			    "v": "0100000001000000",
			    "w": "0101000001001000",
			    "x": "0101010100000000",
			    "y": "0000101000011100",
			    "z": "1100000000100000",}

