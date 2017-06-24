
def tolerable(first_point, second_point, tolerance):
    if (abs(2*(first_point-second_point)/(first_point+second_point)) <= tolerance): return True
    else: return False