#!/usr/bin/env python3

import unittest
import numpy as np

from simulations import *


class TestEvents(unittest.TestCase):
    def setUp(self):
        np.random.seed(2024)

    def test_time_to_next(self):
        s = 1e-2
        theta = 1e-3
        r = 1e-1
        num_samples = 1000
        # Timesteps should get shorter with increasing k
        t0 = sum(time_to_next(0, s, theta, r) for _ in range(num_samples))
        t1 = sum(time_to_next(1, s, theta, r) for _ in range(num_samples))
        t2 = sum(time_to_next(2, s, theta, r) for _ in range(num_samples))
        self.assertGreater(t0, t1)
        self.assertGreater(t1, t2)

    def test_choose_event(self):
        s = 1e-2
        theta = 1e-2
        r = 1e-1
        num_samples = 1000
        # Only mutations when k = 0
        with self.subTest(k=0):
            events_0 = [choose_event(0, s, theta, r) for _ in range(num_samples)]
            self.assertEqual(set(events_0), set(["m"]))
        # When k > 0, should sample all events
        with self.subTest(k=1):
            events_1 = [choose_event(1, s, theta, r) for _ in range(num_samples)]
            self.assertEqual(set(events_1), set(["b", "d", "m", "s"]))

    def test_sampling(self):
        k = 2
        N = 100
        n = 10
        max_allele_count = 10
        sample_temp = sample_sfs(k,N,n,max_allele_count)
        self.assertAlmostEqual(np.sum(sample_temp), 1, delta = 1e-8)
        self.assertEqual(len(sample_temp),max_allele_count+1)


class TestLocations(unittest.TestCase):
    def setUp(self):
        self.L = 1.0
        self.locations = np.array(
            [[0.24, 0.81], [np.nan, np.nan], [0.01, 0.01], [0.99, 0.99]]
        )

    def test_get_alive(self):
        self.assertEqual(list(get_alive(self.locations)), [0, 2, 3])

    def test_get_nan(self):
        self.assertEqual(list(get_nan(self.locations)), [1])

    def test_extend_locations(self):
        orig_len = self.locations.shape[0]
        extended = extend_locations(self.locations)
        self.assertEqual(extended.shape[0], orig_len * 2)
        # The beginning of extended should be equal to the original
        self.assertTrue(
            np.array_equal(extended[:orig_len], self.locations, equal_nan=True)
        )
        self.assertTrue(np.all(np.isnan(extended[orig_len:])))

    def test_wrap_locations(self):
        # Wrapping shouldn't do anything if we're already in bounds
        with self.subTest("self"):
            self.assertTrue(
                np.array_equal(
                    wrap_locations(self.locations, self.L),
                    self.locations,
                    equal_nan=True,
                )
            )

        # If we transpose by L in both dimensions and wrap, should get the same locations
        with self.subTest("+L"):
            self.assertTrue(
                np.allclose(
                    wrap_locations(self.locations + self.L, self.L),
                    self.locations,
                    equal_nan=True,
                )
            )

        # If we transpose by -L in both dimensions and wrap, should get the same locations
        with self.subTest("-L"):
            self.assertTrue(
                np.allclose(
                    wrap_locations(self.locations - self.L, self.L),
                    self.locations,
                    equal_nan=True,
                )
            )

    def test_update_locations(self):
        # Choose a big sigma to check that we're wrapping
        sigma = 0.5
        t_next = 1.0
        num_tests = 100
        for i in range(num_tests):
            with self.subTest(i=i):
                updated = update_locations(self.locations, sigma, t_next, self.L)

                with self.subTest("preserves alive"):
                    self.assertEqual(
                        list(get_alive(updated)), list(get_alive(self.locations))
                    )

                with self.subTest("preserves dead"):
                    self.assertEqual(
                        list(get_nan(updated)), list(get_nan(self.locations))
                    )

                with self.subTest("in bounds"):
                    self.assertTrue(
                        np.array_equal(
                            wrap_locations(updated, self.L),
                            updated,
                            equal_nan=True,
                        )
                    )

if __name__ == "__main__":
    unittest.main()
