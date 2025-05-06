import pygel3d as gel
import unittest

class TestClass(unittest.TestCase):

    def test_succeed(self):
        print(gel.__all__)
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()