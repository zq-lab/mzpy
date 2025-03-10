'''
basic transformer or process for MS spectra (ms)

In general ms is a 2D list or numpy matrix (2D array)
    1st column is mz
    2nd column is intensity
'''

def match(mz1, mz2, tol=0.003, tol_rel=5, mode='abs'):
    '''
    Determine if two mz values (mz1 and mz2) match.
    param:
        mz1, mz2, mz values to compare.
        tol, tolerance. If mode is set as 'both', both tol and tol_rel are required.
    return:
        True or False
    '''
    if mode == 'abs':
        return abs(mz1 - mz2) < tol
    elif mode == 'rel':
        return (1.0E6 * abs(mz1 - mz2) / mz2) < tol
    elif mode == 'both': 
        absolute = abs(mz1-mz2)
        rel = 1.0E6 * absolute / mz2
        return (absolute < tol) or (rel < tol_rel)


def match_list(que, ref, tol=0.003):  
    """  
    Matches two m/z lists (que and ref) within a specified tolerance and   
    returns the unmatched and matched elements from each list.  

    1. que and ref are both sorted in descending order at the start of the function.  
    2. The algorithm uses a two-pointer approach with indices i and j:  
       - If que_sorted[i] and ref_sorted[j] differ by at most 'tolerance', they are considered matched.  
         Both elements are then stored in que_matched and ref_matched, and i and j are incremented.  
       - Otherwise, if que_sorted[i] is significantly larger than ref_sorted[j], it cannot match that  
         (or any smaller) ref element. que_sorted[i] is marked unmatched and i is incremented.  
       - Conversely, if ref_sorted[j] is significantly larger, ref_sorted[j] is marked unmatched and j is incremented.  
    3. Any remaining elements in que_sorted or ref_sorted after one pointer reaches the end are all unmatched.  
    4. The function returns:  
       - que_unmatched: elements in que that find no match in ref  
       - que_matched: elements in que that successfully match with ref  
       - ref_unmatched: elements in ref that find no match in que  
       - ref_matched: elements in ref that successfully match with que  

    Parameters:  
    - que (list[float]): List of m/z values.  
    - ref (list[float]): Another list of m/z values to compare with.  
    - tolerance (float): Allowed absolute difference for matching. Defaults to 0.003.  

    Returns:  
    - Tuple[list[float], list[float], list[float], list[float]]:  
      (que_unmatched, que_matched, ref_unmatched, ref_matched)  

    Example Usage:  
        >>> que = [160.04518, 99.53137, 196.00356, 176.02272, 133.97444, 111.03919]  
        >>> ref = [154.98994, 201.97528, 111.03919, 112.01798, 115.96407, 156.99043, 158.98581]  
        >>> q_unmatched, q_matched, r_unmatched, r_matched = match_arrays_sorted_desc(que, ref)  
        >>> print("Que Unmatched:", q_unmatched)  
        >>> print("Que Matched:", q_matched)  
        >>> print("Ref Unmatched:", r_unmatched)  
        >>> print("Ref Matched:", r_matched)  
    """  
    que_sorted = sorted(que, reverse=True)  
    ref_sorted = sorted(ref, reverse=True)  

    i, j = 0, 0  
    que_matched, que_unmatched = [], []  
    ref_matched, ref_unmatched = [], []  

    len_que, len_ref = len(que_sorted), len(ref_sorted)  

    while i < len_que and j < len_ref:  
        diff = que_sorted[i] - ref_sorted[j]  

        # 在容差范围内 -> 视为匹配  
        if abs(diff) <= tol:  
            que_matched.append(que_sorted[i])  
            ref_matched.append(ref_sorted[j])  
            i += 1  
            j += 1  
        else:  
            if diff > 0:  
                # que_sorted[i] > ref_sorted[j] 太多，i 指针前进  
                que_unmatched.append(que_sorted[i])  
                i += 1  
            else:  
                # ref_sorted[j] > que_sorted[i] 太多，j 指针前进  
                ref_unmatched.append(ref_sorted[j])  
                j += 1  

    # 若还有剩余元素，则这些都无法匹配  
    while i < len_que:  
        que_unmatched.append(que_sorted[i])  
        i += 1  

    while j < len_ref:  
        ref_unmatched.append(ref_sorted[j])  
        j += 1  

    return que_unmatched, que_matched, ref_unmatched, ref_matched  


def unique(mz_values, tol=0.003):  
    """  
    Deduplicate the mass-to-charge ratio (m/z) list:
    If the difference between any two m/z values is less than (or equal to) the tolerance, 
        they are considered duplicates and only one is kept
    
    tol, tolerance
    """  

    # 先对输入 m/z 列表按升序排序  
    mz_sorted = sorted(mz_values)  
    # 将第一项纳入结果列表  
    unique_list = [mz_sorted[0]]  

    # 逐个遍历，若与 unique_list 最后一个值的差大于 tolerance，才加入列表  
    for i in range(1, len(mz_sorted)):  
        if abs(mz_sorted[i] - unique_list[-1]) > tol:  
            unique_list.append(mz_sorted[i])  

    return unique_list 