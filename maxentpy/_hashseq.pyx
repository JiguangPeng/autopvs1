cdef int base_to_int(char base):
    if base == 65: # ord('A')
        return 0
    elif base == 67: # ord('C')
        return 1
    elif base == 71: # ord('G')
        return 2
    elif base == 84: # ord('T')
        return 3
    else:
        raise KeyError()

cpdef int hashseq(str fa):
    cdef int length = len(fa)
    cdef int hashseq = 0
    cdef int i
    cdef int val
    for i in range(length):
        val = base_to_int(ord(fa[i]))
        hashseq += val * 4**(length - i - 1)
    return hashseq
