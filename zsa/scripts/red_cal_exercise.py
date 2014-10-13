#! /usr/bin/env python 
import numpy as n

coeffs = n.array([0.00  ,-8.62 ,-21.13, -0.96,-8.29 ,-11.23, -16.91, -18.56])

tau_ew_meas = n.array([-0.55, 5.74, 0.12, 2.7, -.15, -.82, 5.52])
tau_ew_true = 2.23
ants_layout_row1 =  n.array([49,41,47,19,29,28,34,51])

dly_diffs = coeffs[1:] - coeffs[:7] + tau_ew_meas

#0 is the first antenna coeffs when we mandate that tau_ew = 2.2 (true value)
true_coeffs = [0]

for i in range(len(dly_diffs)):
    true_c_prev = true_coeffs[-1]
    true_coeffs.append(dly_diffs[i] + true_c_prev - tau_ew_true)

print true_coeffs






