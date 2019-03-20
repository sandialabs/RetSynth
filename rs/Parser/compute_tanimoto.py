# compute_tanimoto.py
import sys
import array

# Create a lookup table from byte value to 'on' bit count
# chr(0) has 0 'on' bits
# chr(1) has 1 'on' bit
# chr(2) has 1 'on' bit
# chr(3) has 2 'on' bits
#  ...
# chr(255) has 8 'on' bits
popcount_in_byte = (
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
    )


def tanimoto(query, target):
    a = b = c = 0
    for query_byte, target_byte in zip(query, target):
        a += popcount_in_byte[query_byte]
        b += popcount_in_byte[target_byte]
        c += popcount_in_byte[query_byte & target_byte]
    # I'll define that if neither fingerprint has any
    # set bits (so a == b == 0) then the Tanimoto is 0.
    # Otherwise the equation will try to compute 0/0
    # and Python will throw an exception.
    if not c:
        return 0
    return float(c)/(a+b-c)

# if len(sys.argv) > 1:
#     filename = sys.argv[1]
# else:
#     filename = "drugs.binary_fp"

# # Use the first fingerprint as the query
# infile = open(filename, "rb")
# s = infile.read(1024//8)  # 128 bytes in a fingerprint
# query = array.array("b", s)

# # Reopen and compute the Tanimoto against all fingerprints
# # including itself.
# infile = open(filename, "rb")

# while 1:
#     s = infile.read(1024//8)
#     if not s:
#         # End of file
#         break

#     target = array.array("b", s)
#     print tanimoto(query, target)