#!/usr/bin/env python3

import unittest
from copy import copy

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
            self.assertEqual(set(events_0), set([Event.MUTATION]))
        # When k > 0, should sample all events
        with self.subTest(k=1):
            events_1 = [choose_event(1, s, theta, r) for _ in range(num_samples)]
            self.assertEqual(set(events_1), set(Event))

    def test_sampling_prob(self):

        with self.subTest("in range"):
            L = 50
            w = 10
            rho = 2
            locations = np.array(
                [[0.24, 0.81], [np.nan, np.nan], [0.01, 0.01], [0.99, 0.99]]
            )
            centers = get_centers_grid(L,2)
            p_temp = sampling_probability_gaussian(locations, centers, w, L, rho)
            self.assertGreaterEqual(p_temp[0], 0)
            self.assertGreaterEqual(p_temp[1], 0)
            self.assertGreaterEqual(p_temp[2], 0)
            self.assertGreaterEqual(p_temp[3], 0)
            self.assertLessEqual(p_temp[0], 1)
            self.assertLessEqual(p_temp[1], 1)
            self.assertLessEqual(p_temp[2], 1)
            self.assertLessEqual(p_temp[3], 1)
            self.assertEqual(len(p_temp),4)

        with self.subTest("add carrier"):
            L = 1.0
            w = 0.1
            rho = 2
            locations = np.array(
                [[0.24, 0.81], [np.nan, np.nan], [0.01, 0.01], [0.99, 0.99]]
            )
            centers=[(L/2,L/2)]
            p_temp = sampling_probability_gaussian(locations, centers, w, L, rho)
            loc_new = np.append(locations, [[0, 0]], axis=0)
            p_new = sampling_probability_gaussian(loc_new, centers, w, L, rho)
            self.assertGreater(p_new[0], p_temp[0])

        with self.subTest("farther away"):
            L = 1.0
            w = 0.1
            rho = 2
            centers = [(L/2, L/2)]
            locations = np.array([[0.5, 0.5], [0.5, 0.5], [0.4, 0.4]])
            p_temp = sampling_probability_gaussian(locations, centers, w, L, rho)
            loc_new = locations - 0.2
            p_new = sampling_probability_gaussian(loc_new, centers, w, L, rho)
            self.assertLess(p_new[0], p_temp[0])

        with self.subTest("increase w"):
            L = 1.0
            w = 0.1
            rho = 2
            centers = [(L / 2, L / 2)]
            locations = np.array([[0.45, 0.45]])
            p_temp = sampling_probability_gaussian(locations, centers,w, L, rho)
            p_new = sampling_probability_gaussian(locations, centers,w * 2, L, rho)
            self.assertLess(p_new[0], p_temp[0])

        with self.subTest("one at boundary: 0"):
            L = 1.0
            w = 0.1
            rho = 2
            locations = np.array([[0, 0]])
            centers = [(L / 2, L / 2)]
            p_temp = sampling_probability_gaussian(locations, centers,w, L, rho)
            self.assertGreater(p_temp[0], 0)

        with self.subTest("one at boundary: 1"):
            L = 1.0
            w = 0.1
            rho = 2
            locations = np.array([[1, 1]])
            centers = [(L / 2, L / 2)]
            p_temp = sampling_probability_gaussian(locations, centers,w, L, rho)
            self.assertGreater(p_temp[0], 0)

        with self.subTest("one at boundary: epsilon away"):
            L = 1.0
            w = 0.1
            rho = 2
            epsilon = 1e-12
            locations = np.array([[1 - epsilon, 1 - epsilon]])
            centers = [(L / 2, L / 2)]
            p_temp = sampling_probability_gaussian(locations, centers,w, L, rho)
            self.assertGreater(p_temp[0], 0)

        with self.subTest("no carriers"):
            L = 1.0
            w = 0.1
            rho = 2
            locations = np.array([[np.nan, np.nan], [np.nan, np.nan]])
            centers = [(L / 2, L / 2)]
            p_temp = sampling_probability_gaussian(locations, centers,w, L, rho)
            self.assertEqual(p_temp[0], 0)

        with self.subTest("near edge"):
            L = 1.0
            w = 0.1
            rho = 2
            epsilon = 1e-12
            locations = np.array([[1 - epsilon, 1 - epsilon]])
            centers = [(L / 2, L / 2)]
            p_temp = sampling_probability_gaussian(locations, centers,w, L, rho)
            p_new = sampling_probability_gaussian(locations, centers,w * 2, L, rho)
            self.assertGreater(p_new[0], p_temp[0])

    def test_grid_centers(self):
        L = 50
        with self.subTest("return origin for single center"):
            self.assertEqual(
                get_centers_grid(L, 1), [(0, 0)]
            )
        with self.subTest("correct length"):
            self.assertEqual(
                len(get_centers_grid(L, 12)), 144
            )

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
        extended = copy(self.locations)
        extend_locations(extended)
        self.assertEqual(extended.shape[0], orig_len * 2)
        # The beginning of extended should be equal to the original
        self.assertTrue(
            np.array_equal(extended[:orig_len], self.locations, equal_nan=True)
        )
        self.assertTrue(np.all(np.isnan(extended[orig_len:])))

    def test_get_free_row(self):
        with self.subTest("free"):
            to_check = copy(self.locations)
            self.assertEqual(get_free_row(to_check), 1)
            # If there is a free row, should not change locations
            self.assertTrue(np.array_equal(to_check, self.locations, equal_nan=True))

        with self.subTest("no free"):
            full = np.array([[0.5, 0.5]])
            self.assertEqual(get_free_row(full), 1)
            # If there is not a free row, should double the length and fill with nan
            self.assertTrue(
                np.array_equal(full, [[0.5, 0.5], [np.nan, np.nan]], equal_nan=True)
            )

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
