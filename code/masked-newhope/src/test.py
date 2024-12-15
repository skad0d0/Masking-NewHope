from math import floor

def flipabs(r):
    NEWHOPE_Q = 12289  # Replace with the actual value of NEWHOPE_Q
    r = r - (NEWHOPE_Q // 2)
    m = r >> 15  # In Python, right-shifting a negative number behaves as expected
    return (r + m) ^ m

def to_unsigned(t, bits):
    """Convert t to an unsigned integer with the specified number of bits."""
    mask = (1 << bits) - 1  # Create a mask for the desired bit size
    return t & mask

q = 12289
a1 = 3*q//4+1
a2 = 3*q//4

a1 = flipabs(a1)
# print (a1)
a2 = flipabs(a2)

# t = 0x2686
# print(t)
t = a1 + a2
t = int(t - q//2)
t = to_unsigned(t, 16)
t = int(t) >> 15

print(t)
print("--------------------")

kq = 3329
r1 = kq - 1
r1 = ((r1 << 1)+kq//2) // kq
r1 = int(r1)
print(r1)
print(r1&1)