#!/usr/bin/env python
# Written by Magnus Lie Hetland, modified by Yunfan Zhang

def partition(array, start, end):
    pivot = array[end]                          # Partition around the last value
    pivotval = array[end][0]
    largest = start-1                           # Start outside the area to be partitioned
    smallest = end                                  # Ditto
    done = 0
    while not done:                            # Until all elements are partitioned...
        while not done:                        # Until we find an out of place element...
            largest = largest+1                  # ... move the largest up.
            if largest == smallest:                  # If we hit the smallest...
                done = 1                       # ... we are done.
                break
            if array[largest][0] < pivotval:           # Is the largest out of place?
                array[smallest] = array[largest]       # Then put it at the smallest...
                break                          # ... and start searching from the smallest.
        while not done:                        # Until we find an out of place element...
            smallest = smallest-1                        # ... move the smallest down.
            if smallest == largest:                  # If we hit the largest...
                done = 1                       # ... we are done.
                break
            if array[smallest][0] > pivotval:              # Is the smallest out of place?
                array[largest] = array[smallest]       # Then put it at the largest...
                break                          # ...and start searching from the largest.
    array[smallest] = pivot                          # Put the pivot in its place.
    return smallest                                 # Return the split point

def quick_sort(array, start, end):
    if start < end:                            # If there are two or more elements...
        split = partition(array, start, end)    # ... partition the subarray...
        quick_sort(array, start, split-1)        # ... and sort both halves.
        quick_sort(array, split+1, end)
    else:
        return
