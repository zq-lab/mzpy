
"""
Basic processors for m/z value(s)

This module provides utilities for processing mass-to-charge ratio (m/z) values 
or numpy arrays of m/z values. It includes functionality for determining whether 
two m/z values match, and for matching between two groups of m/z values.
It does not include intensity information processing.
"""

import numpy as np

def match(mz1, mz2, tol=0.003, is_abs_error=True):
    """Determine if two m/z values match within specified tolerance.

    Args:
        mz1 (float): First m/z value to compare.
        mz2 (float): Second m/z value to compare.
        tol (float): Tolerance threshold. Defaults to 0.003 Da (absolute error).
                     Absolute error (when is_abs_error is True), the error unit is Da (Dalton). 
                        Default value is 0.003; 
                     Relative error (when is_abs_error is False), the error unit is ppm, 
                        typically set to 5.
        is_abs_error (bool): If True, use absolute difference; 
                           if False, use relative ppm difference. Defaults to True.

    Returns:
        bool: True if the values match within tolerance, False otherwise.
    """
    if is_abs_error:
        return abs(mz1 - mz2) < tol
    else:
        return (1.0E6 * abs(mz1 - mz2) / mz2) < tol


def match_array(que, ref, tol=0.003):  
    """Match two m/z arrays within specified tolerance.
    
    Matches two m/z lists (que and ref) within a specified tolerance and   
    returns the unmatched and matched elements from each list.  
    
    Args:  
        que (np.ndarray): Array of m/z values.  
        ref (np.ndarray): Another array of m/z values to compare with.  
        tol (float): Allowed absolute difference for matching. Defaults to 0.003.  
    
    Returns:  
        tuple: Contains four numpy arrays in order:
               (que_unmatched, que_matched, ref_unmatched, ref_matched)  
    """
    # Sort input arrays in descending order
    que_sorted = np.sort(que)[::-1]
    ref_sorted = np.sort(ref)[::-1]
    
    len_que, len_ref = len(que_sorted), len(ref_sorted)
    
    # Pre-allocate result arrays with maximum possible length
    que_matched_list = np.zeros(min(len_que, len_ref), dtype=np.float64)
    que_unmatched_list = np.zeros(len_que, dtype=np.float64)
    ref_matched_list = np.zeros(min(len_que, len_ref), dtype=np.float64)
    ref_unmatched_list = np.zeros(len_ref, dtype=np.float64)
    
    # Counters
    que_matched_count = 0
    que_unmatched_count = 0
    ref_matched_count = 0
    ref_unmatched_count = 0
    
    i, j = 0, 0  
    while i < len_que and j < len_ref:  
        diff = que_sorted[i] - ref_sorted[j]  
        
        # Within tolerance -> consider as match
        if abs(diff) <= tol:  
            que_matched_list[que_matched_count] = que_sorted[i]
            ref_matched_list[ref_matched_count] = ref_sorted[j]
            que_matched_count += 1
            ref_matched_count += 1
            i += 1  
            j += 1  
        else:  
            if diff > 0:  
                # que_sorted[i] > ref_sorted[j] by too much, advance i pointer
                que_unmatched_list[que_unmatched_count] = que_sorted[i]
                que_unmatched_count += 1
                i += 1  
            else:  
                # ref_sorted[j] > que_sorted[i] by too much, advance j pointer
                ref_unmatched_list[ref_unmatched_count] = ref_sorted[j]
                ref_unmatched_count += 1
                j += 1  
    
    # Add any remaining elements as unmatched
    while i < len_que:  
        que_unmatched_list[que_unmatched_count] = que_sorted[i]
        que_unmatched_count += 1
        i += 1  
    
    while j < len_ref:  
        ref_unmatched_list[ref_unmatched_count] = ref_sorted[j]
        ref_unmatched_count += 1
        j += 1  
    
    # Trim arrays to actual size
    que_unmatched = que_unmatched_list[:que_unmatched_count]
    que_matched = que_matched_list[:que_matched_count]
    ref_unmatched = ref_unmatched_list[:ref_unmatched_count]
    ref_matched = ref_matched_list[:ref_matched_count]
    
    return que_unmatched, que_matched, ref_unmatched, ref_matched


def match_num(que, ref, tol=0.003):
    """Count matching elements between two m/z arrays.
    
    Calculates the number of elements that match between two m/z lists 
    within the specified tolerance.
    
    Args:
        que (np.ndarray): Array of m/z values.
        ref (np.ndarray): Another array of m/z values to compare with.
        tol (float): Allowed absolute difference for matching. Defaults to 0.003.
    
    Returns:
        int: Number of matching elements.
    """
    # Sort input arrays in descending order
    que_sorted = np.sort(que)[::-1]
    ref_sorted = np.sort(ref)[::-1]
    
    len_que, len_ref = len(que_sorted), len(ref_sorted)
    
    # Match counter
    match_count = 0
    
    i, j = 0, 0
    while i < len_que and j < len_ref:
        diff = que_sorted[i] - ref_sorted[j]
        
        # Within tolerance -> consider as match
        if abs(diff) <= tol:
            match_count += 1
            i += 1
            j += 1
        else:
            if diff > 0:
                # que_sorted[i] > ref_sorted[j] by too much, advance i pointer
                i += 1
            else:
                # ref_sorted[j] > que_sorted[i] by too much, advance j pointer
                j += 1
    
    return match_count


def unique(mz_values, tol=0.003):  
    """Deduplicate a list of m/z values.
    
    If the difference between any two m/z values is less than (or equal to) 
    the tolerance, they are considered duplicates and only one is kept.
    
    Args:
        mz_values (list or np.ndarray): List of m/z values to deduplicate.
        tol (float): Tolerance threshold for considering values as duplicates.
                    Defaults to 0.003.
    
    Returns:
        list: Deduplicated list of m/z values.
    """
    # Sort the input m/z list in ascending order
    mz_sorted = sorted(mz_values)  
    # Include the first item in the result list
    unique_list = [mz_sorted[0]]  

    # Add values only if they differ from the last added value by more than the tolerance
    for i in range(1, len(mz_sorted)):  
        if abs(mz_sorted[i] - unique_list[-1]) > tol:  
            unique_list.append(mz_sorted[i])  

    return unique_list


