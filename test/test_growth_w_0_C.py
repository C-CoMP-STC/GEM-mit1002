import unittest

class TestGrowthWithNoCarbon(unittest.TestCase):
    def test_growth_w_0_C(self):
        self.assertEqual('foo'.upper(), 'FOO')

if __name__ == '__main__':
    unittest.main()
