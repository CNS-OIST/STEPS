
def tolerable(a, b, tolerance):
    return a == b == 0 or (a + b != 0 and abs(2 * (a - b) / (a + b)) <= tolerance)

