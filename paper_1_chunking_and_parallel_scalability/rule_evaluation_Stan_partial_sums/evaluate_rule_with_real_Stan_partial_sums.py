##
## ======================================================================================================================
## evaluate_rule_with_real_Stan_partial_sums.py
##
## Paper 1, table:paper1_auto_n_chunks_rule_evaluation, Stan half: the same evaluation as
## rule_evaluation_v45/evaluate_rule_split_algorithm1_check.py (rule alone, and the check during burn-in of algorithm 1),
## but with every Stan N_chunks replaced by the number of partial sums which reduce_sum_static actually ran.
##
## reduce_sum_static (TBB simple partitioner on a blocked range) halves the range of N individuals until each piece holds
## at most chunk_size = ceiling(N / N_chunks) individuals, so a requested N_chunks runs as partial_sums(N, N_chunks)
## pieces (the next power of two for the tested values). Two requested values which ran as the same partial sums
## (e.g. WCP-only with 22 threads per chain and N_chunks = 25, both 32 partial sums) are the same configuration timed
## twice; their mean time is used.
##
## The rule's N_chunks for Stan (L3-cache per active thread, as before) runs as partial_sums(N, rule N_chunks); the check
## during burn-in doubles or halves the number of partial sums, never below the WCP-only partition (one requested chunk
## per thread). NicoStan+BayesMVP is unchanged (its chunks are its own loop).
##
## Run with:  cd <this folder> && nice -n 19 python3 evaluate_rule_with_real_Stan_partial_sums.py
##
import collections, csv, math

def partial_sums(N, N_chunks):
    grainsize = math.ceil(N / N_chunks)
    def leaves(n):
        if n <= grainsize: return 1
        return leaves(n // 2) + leaves(n - n // 2)
    return leaves(N)

rows = list(csv.DictReader(open('curves_extended_chunk_grid.csv')))
alloc = collections.defaultdict(list)
for r in rows:
    if r['case_set'] == 'chunking_only_tape_chunked_reference': continue
    alloc[r['allocation']].append(r)
HW = {'HPC': dict(L3=2**25, nL3=12, cores=8), 'Laptop': dict(L3=2**24, nL3=1, cores=8)}
B_row_Stan = 19152
keep = {'HPC': (96, 176, 180), 'Laptop': (8, 16)}

def curve_real(rs, N):
    times = collections.defaultdict(list)
    for r in rs:
        times[partial_sums(N, int(r['num_chunks']))].append(float(r['time_mean']))
    x = sorted(times); t = [sum(times[k]) / len(times[k]) for k in x]
    best = min(t); return x, [best / ti for ti in t]

def rel(x, y, c):
    if c <= x[0]: return y[0]
    if c >= x[-1]: return y[-1]
    if c in x: return y[x.index(c)]
    i = max(j for j in range(len(x)) if x[j] < c)
    w = (math.log(c) - math.log(x[i])) / (math.log(x[i+1]) - math.log(x[i]))
    return y[i] + w * (y[i+1] - y[i])

def k_of(d, n): return min(16, math.ceil(n / HW[d]['nL3']))

def rule_requested(d, n, t, W):
    c = math.ceil(W / (HW[d]['L3'] / k_of(d, n)))                     ## L3-cache per active thread
    return t * math.ceil(max(c, t) / t)                                ## WCP rounding (equation eq:paper1_auto_n_chunks_WCP)

def check_algorithm_1(x, y, s, lo):
    probed = {s: rel(x, y, s)}; current = s
    candidates = [2 * s] + ([s // 2] if s // 2 >= lo and s // 2 != s else [])
    for c in candidates: probed.setdefault(c, rel(x, y, c))
    best_first = max(candidates, key = lambda c: probed[c])
    if probed[best_first] > probed[current]:
        direction = 'up' if best_first > current else 'down'; current = best_first
        while True:
            nxt = 2 * current if direction == 'up' else current // 2
            if nxt == current or nxt < lo: break
            probed.setdefault(nxt, rel(x, y, nxt))
            if probed[nxt] > probed[current]: current = nxt
            else: break
    chosen = max(probed, key = lambda c: probed[c])
    return chosen

per_case = []
for key, rs in alloc.items():
    impl, cs, dev, N, nch, tpc = key.split('|')
    if impl != 'Stan': continue
    N, nch, tpc = int(N), int(nch), int(tpc)
    if nch * tpc not in keep[dev]: continue
    x, y = curve_real(rs, N)
    s_req = rule_requested(dev, nch * tpc, tpc, N * B_row_Stan)
    s = partial_sums(N, s_req)
    lo = x[0]
    chk = check_algorithm_1(x, y, s, lo)
    per_case.append(dict(case_set = cs, device = dev, N = N, n_chains = nch, threads_per_chain = tpc, n_threads = nch * tpc,
                         tested_partial_sums = ';'.join(map(str, x)), best_partial_sums = x[y.index(max(y))],
                         rule_requested_N_chunks = s_req, rule_partial_sums = s, rule_rel = round(rel(x, y, s), 4),
                         check_partial_sums = chk, check_rel = round(rel(x, y, chk), 4)))
per_case.sort(key = lambda r: (r['case_set'], r['device'], r['n_threads'], r['N'], r['n_chains']))
with open('rule_evaluation_Stan_real_partial_sums_per_case.csv', 'w', newline = '') as f:
    w = csv.DictWriter(f, fieldnames = list(per_case[0].keys())); w.writeheader(); w.writerows(per_case)

def pct(k, n): return f"{k} (100\\%)" if k == n else f"{k} ({100 * k / n:.1f}\\%)"
def cells(sel):
    n = len(sel); g = [r['rule_rel'] for r in sel]; c = [r['check_rel'] for r in sel]
    return (f"{n} & {pct(sum(v >= 0.95 for v in g), n)} & {pct(sum(v >= 0.90 for v in g), n)} & {min(g):.3f} & "
            f"{pct(sum(v >= 0.95 for v in c), n)} & {pct(sum(v >= 0.90 for v in c), n)} & {min(c):.3f} \\\\")
groups = {('HPC', 96): (96,), ('HPC', 180): (176, 180), ('Laptop', 8): (8,), ('Laptop', 16): (16,)}
for cs in ('chunking_only', 'chunking_WCP'):
    print('%%', cs)
    for (dev, nt), nts in groups.items():
        sel = [r for r in per_case if r['case_set'] == cs and r['device'] == dev and r['n_threads'] in nts]
        if sel: print(f"  & {dev} {nt:<4} & " + cells(sel))
print("  & Total & " + cells(per_case))
print('mean rule_rel', round(sum(r['rule_rel'] for r in per_case) / len(per_case), 3))
print('outside 5% (rule alone):', sum(r['rule_rel'] < 0.95 for r in per_case),
      '| of these, rule fewer partial sums than best:', sum(r['rule_rel'] < 0.95 and r['rule_partial_sums'] < r['best_partial_sums'] for r in per_case))
print('worked example (N = 50,000, chunking only):')
for r in per_case:
    if r['case_set'] == 'chunking_only' and r['N'] == 50000:
        print('  ', r['device'], r['n_threads'], 'rule', r['rule_requested_N_chunks'], '->', r['rule_partial_sums'], r['rule_rel'],
              'check', r['check_partial_sums'], r['check_rel'], 'best', r['best_partial_sums'])
