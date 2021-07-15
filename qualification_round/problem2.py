from unicodedata import lookup
from tqdm import tqdm
import concurrent.futures as cf
import pandas as pd
import numpy as np

# Problem 2: Metabolite Annotation

def build_lookup_list(masses, adducts):
    lookup = {}
    for i in tqdm(range(len(adducts))):
        for j in range(len(masses)):
            total = round(adducts[i] + masses[j],6)
            combo = f'{j+1} {i+1}'
            if total not in lookup:
                lookup[total] = [combo]
            # else:
            #     lookup[total].append(combo)
    lookup = np.array([
        [key for key in lookup.keys()], 
        [val for val in lookup.values()]
        ], dtype='object')
    lookup = lookup[:,lookup[0,:].argsort()]
    return lookup

def lookup_combo(lst, signal):
    masses = lst[0]
    combos = lst[1]
    if len(masses) <= 5:
        delta = float('inf')
        combo = ''
        for i in range(len(masses)):
            delt = abs(signal - masses[i])
            if delt < delta:
                delta = delt
                combo = combos[i]
        return combo
    midpoint = round(float(len(masses)/2))
    if signal >= lst[0,:][midpoint]:
        return lookup_combo(lst[:, midpoint:], signal)
    else:
        return lookup_combo(lst[:,:midpoint+1], signal)


# attempts ====================================================================

# zug zug exponential attempt ---------------------------------
# def annotate_metabolites(case):
#     results = []
#     for i in range(0, len(case['s'])):
#         result = {'combo': '', 'delta': float('inf')}
#         for j in range(0, len(case['m'])):
#             for k in range(0, len(case['a'])):
#                 mass = case['m'][j] + case['a'][k]
#                 if mass > 0:
#                     delta = abs(case['s'][i] - mass)
#                     if delta < result['delta']:
#                         result['delta'] = delta
#                         result['combo'] = f'{j+1} {k+1}'
#         results.append(f'{result["combo"]}')
#     return results
# ----------------------------------------------------------------

# flawed gradient descent attempt --------------------------------
# def test_combo(s, j, k):
#     if j and k and j+k > 0:
#         return abs(s-(j+k))
#     else:
#         return float('inf')
# 
# def decend_gradient(s, m, a):
#     j,k = 0,0
#     done = False
#     while not done:
#         options = {
#             'current_square': test_combo(s, m[j], a[k]),
#             'step_right': float('inf'),
#             'step_down': float('inf'),
#             'step_diag': float('inf'),
#             }
#         space_to_right = j < len(m)-1
#         space_to_down = k < len(a)-1
#         if space_to_right:
#             options['step_right'] = test_combo(s, m[j+1], a[k])
#         if space_to_down:
#             options['step_down'] = test_combo(s, m[j], a[k+1])
#         if space_to_right and space_to_down:
#             options['step_diag']=  test_combo(s, m[j+1], a[k+1])
# 
#         option = min(options, key=options.get)
#         if option == 'current_square':
#             return [option, m[j], a[k]]
#         elif option == 'step_right':
#             j += 1
#         elif option == 'step_down':
#             k += 1
#         elif option == 'step_diag':
#             j += 1
#             k += 1
# 
# def annotate_metabolites(case):
#     sorted_m = sorted(case['m'], reverse=True)
#     sorted_m_2 = sorted(case['m'])
#     sorted_a = sorted(case['a'], reverse=True)
# 
#     results = []
#     for i in range(0, len(case['s'])):
#         choice_1 = decend_gradient(case['s'][i], sorted_m, sorted_a,)
#         choice_2 = decend_gradient(case['s'][i], sorted_m_2, sorted_a,)
#         if choice_1[0] > choice_2[0]:
#             m_j = choice_1[1]
#             a_k = choice_1[2]
#         else:
#             m_j = choice_2[1]
#             a_k = choice_2[2]
#         j = case['m'].index(m_j) + 1
#         k = case['a'].index(a_k) + 1
#         results.append(f'{j} {k}')
#     return results
#----------------------------------------------------------------------

# multithreadded option --------------------------------------------
# def annotate_metabolites(case):
#     sorted_m = sorted(case['m'], reverse=True)
#     sorted_m_2 = sorted(case['m'])
#     sorted_a = sorted(case['a'], reverse=True)
# 
#     results = []
#     with cf.ThreadPoolExecutor() as executor:
#         futures = [executor.submit(decend_gradient, case['s'][i], sorted_m, sorted_a) for i in range(0, len(case['s']))]
#         for future in list(tqdm(cf.as_completed(futures), total = len(futures))):
#             res = future.result()
#             m_j = res[1]
#             a_k = res[2]
#             j = case['m'].index(m_j) + 1
#             k = case['a'].index(a_k) + 1
#             results.append(f'{j} {k}')
#     return results
# -------------------------------------------------------------------------

# greedy matrix traversal ------------------------------------------------
# def build_lookup_table(masses, adducts):
#     table = pd.DataFrame(index = list(range(len(adducts))), columns=list(range(len(masses))))
#     for i in range(len(adducts)):
#         for j in range(len(masses)):
#             table.iloc[i][j] = adducts[i] + masses[j]
#     table.index = [i+1 for i in range(len(table.index))]
#     table.columns = [i+1 for i in range(len(table.columns))]
#     table = table.sort_values(by=list(table.index), axis=1)
#     table = table.sort_values(by=list(table.columns)[-1], axis=0)
#     return table
#
# def lookup_combo(signal, table):
#     print(signal)
#     print(table)
#     print(table.shape)
#     i,j = 0,0
#     done = False
#     while not done:
#         print('loop')
#         options = {
#             'current_square': abs(signal - table.iloc[i,j])
#             }
#         space_to_right = j+1 < table.shape[1]
#         space_to_down = i+1 < table.shape[0]
#         options['step_right'] = abs(signal - table.iloc[i,j+1]) if space_to_right else float('inf')
#         options['step_down'] = abs(signal - table.iloc[i+1,j]) if space_to_down else float('inf')
#         options['step_diag']= abs(signal - table.iloc[i+1,j+1]) if space_to_right and space_to_down else float('inf')
#         print(options)
#         option = min(options, key=options.get)
#         print(option)
#         if option == 'current_square':
#             return f'{table.columns[j]} {table.index[i]}'
#         elif option == 'step_down':
#             i += 1
#         elif option == 'step_right':
#             j += 1
#         elif option == 'step_diag':
#             j += 1
#             i += 1
# -----------------------------------------------------------------

# numpy based lookup build -----------------------------------------
# def build_lookup_list(masses, adducts):
#     sums = np.full(len(masses)*len(adducts), float('inf'))
#     combos = np.full(len(masses)*len(adducts), '0 0')
#     n = 0
#     for i in range(len(masses)):
#         for j in range(len(adducts)):
#             total = round(masses[i] + adducts[j],6)
#             combo = f'{i+1} {j+1}'
#             if total not in sums:
#                 sums[n] = total
#                 combos[n] = combo
#                 n += 1
#     lookup_len = np.argmax(sums == float('inf'))
#     lookup = np.array([sums[:lookup_len], combos[:lookup_len]], dtype='object') if lookup_len > 0 else np.array([sums, combos], dtype='object')
#     lookup = lookup[:,lookup[0,:].argsort()]
#     return lookup
#------------------------------------------------------------------