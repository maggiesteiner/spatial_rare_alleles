#!/usr/bin/env python3

import unittest
import numpy as np

from simulations import *

np.random.seed(2024)


class TestEvents(unittest.TestCase):
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


# def get_alive(locations):
# def get_nan(locations):
# def wrap_locations(locations,L):
# def extend_locations(locations):
# def update_locations(locations,sigma,t_next,L):

if __name__ == "__main__":
    unittest.main()
