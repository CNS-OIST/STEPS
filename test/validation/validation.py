
failed_tests=[]

# # # # # # # # # # # # # # # REACTIONS # # # # # # # # # # # # # # # # #

print "\n1: Reaction - First order, irreversible:"
import first_order_irev
if first_order_irev.passed: print "PASSED\n"
else: 
    print "FAILED\n"
    failed_tests.append("1: Reaction - First order, irreversible")

print "2: Reaction - First order, reversible:"
import first_order_rev
if first_order_rev.passed: print "PASSED\n"
else: 
    print "FAILED\n"
    failed_tests.append("2: Reaction - First order, reversible")

print "3: Reaction - Second order, irreversible, A0=B0:"
import second_order_irev_AA
if second_order_irev_AA.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("3: Reaction - Second order, irreversible, A0=B0")

print "4: Reaction - Second order, irreversible, A0!=B0:"
import second_order_irev_AB
if second_order_irev_AB.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("4: Reaction - Second order, irreversible, A0!=B0")

print "5: Reaction - Second-order, irreversible, 2D:"
import second_order_irev_2D
if second_order_irev_2D.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("5: Reaction - Second order, irreversible, 2D")

print "6: Reaction - Production and degradation:"
import masteq
if masteq.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("6: Reaction - Production and degradation")
    num_failed = 1

# # # # # # # # # # # # # # # DIFFUSION # # # # # # # # # # # # # # # # #

print "7: Diffusion - Unbounded:"
import unbdiff
if unbdiff.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("7: Diffusion - Unbounded")
    num_failed = 1

print "8: Diffusion - Bounded:"
import bounddiff
if bounddiff.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("8: Diffusion - Bounded")
    num_failed = 1

print "9: Diffusion - Clamped:"
import csd_clamp
if csd_clamp.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("9: Diffusion - Clamped")
    num_failed = 1

# # # # # # # # # # # # # #  REACTION-DIFFUSION # # # # # # # # # # # # # 

print "10: Reaction-diffusion - Degradation-diffusion:"
import kisilevich
if kisilevich.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("10: Reaction-diffusion - Degradation-diffusion")
    num_failed = 1

print "11: Reaction-diffusion - Production and second order degradation:"
import masteq_diff
if masteq_diff.passed: print "PASSED\n"
else:
    print "FAILED\n"
    failed_tests.append("11: Reaction-diffusion - Production and second order degradation")
    num_failed = 1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


print ""
print len(failed_tests), "tests failed:"
for t in failed_tests:
    print t, 
