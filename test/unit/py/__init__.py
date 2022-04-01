import unittest

from . import API_1
from . import API_2
from . import dist

def suite():
    all_tests = []
    all_tests.append(API_1.suite())
    all_tests.append(API_2.suite())
    all_tests.append(dist.suite())
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

