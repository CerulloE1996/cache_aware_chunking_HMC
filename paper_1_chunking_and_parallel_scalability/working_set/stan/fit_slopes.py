## Fits the per-individual working set of one reduce_sum_static partial sum from results.txt (the probe output of
## ws_harness): for every single-chunk case (chunk_size = N) at parameter point A, third repetition (steady state),
##   arena  = arena_fn + 144 * rows            (the chunk's tape, plus the deep copy of its u_raw slice: 6 vars x 24 B per row)
##   stacks = 8 * (var_stack_fn + var_nochain_fn + var_nochain_at_enter_minus_outer)   (pointers pushed for the chunk)
##   heap   = heap_peak_fn + heap_at_enter_minus_outer   (peak heap during the forward pass, plus the slice container)
##   data   = 52 * rows                        (the y row, 6 doubles, and the pop entry, 1 int)
## and regresses each on rows. Usage: python3 fit_slopes.py results.txt
import sys
rows = []
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'results.txt'):
    if line.startswith('CHUNK'):
        rows.append(dict(kv.split('=') for kv in line.split()[1:]))
single = [r for r in rows if r['rep'] == '2' and r['label'].startswith('A_N') and r['label'] == 'A_N' + r['rows'] + '_cs' + r['rows']]
single.sort(key = lambda r: int(r['rows']))
comp = {}
for r in single:
    n = int(r['rows'])
    comp.setdefault('arena', []).append((n, int(r['arena_fn']) + 144 * n))
    comp.setdefault('stacks', []).append((n, 8 * (int(r['var_stack_fn']) + int(r['var_nochain_fn']) + int(r['var_nochain_at_enter_minus_outer']))))
    comp.setdefault('heap', []).append((n, int(r['heap_peak_fn']) + int(r['heap_at_enter_minus_outer'])))
    comp.setdefault('data', []).append((n, 52 * n))
comp['total'] = [(n, sum(comp[k][i][1] for k in ('arena', 'stacks', 'heap', 'data'))) for i, (n, _) in enumerate(comp['arena'])]
print('n, arena, stacks, heap, data, total, total/n')
for i, (n, _) in enumerate(comp['arena']):
    vals = [comp[k][i][1] for k in ('arena', 'stacks', 'heap', 'data', 'total')]
    print(n, *vals, round(vals[-1] / n, 1))
def fit(points):
    n = len(points); sx = sum(p[0] for p in points); sy = sum(p[1] for p in points)
    sxx = sum(p[0] ** 2 for p in points); sxy = sum(p[0] * p[1] for p in points)
    slope = (n * sxy - sx * sy) / (n * sxx - sx ** 2); intercept = (sy - slope * sx) / n
    ss_res = sum((p[1] - slope * p[0] - intercept) ** 2 for p in points); mean = sy / n
    ss_tot = sum((p[1] - mean) ** 2 for p in points)
    return slope, intercept, 1 - ss_res / ss_tot
print()
for k, label in [('arena', 'arena (tape: forward-pass varis, callbacks, arena copies + slice deep copy)'),
                 ('stacks', 'AD stacks (var_stack_ + var_nochain_stack_ pointers)'),
                 ('heap', 'heap (Eigen locals, chainable_alloc, slice container)'),
                 ('data', 'data read (y row 48 B + pop 4 B)'), ('total', 'TOTAL')]:
    slope, intercept, r2 = fit(comp[k])
    print(label + ': slope=' + format(slope, '.1f') + ' B/row intercept=' + format(intercept, '.0f') + ' B R2=' + format(r2, '.8f'))
print('arena+stacks: slope=' + format(fit(comp['arena'])[0] + fit(comp['stacks'])[0], '.1f'))
## Several partial sums per gradient: arena per row of each block (should match the single-chunk slope).
for r in rows:
    if r['rep'] == '2' and r['label'] in ('A_N20000_cs5000', 'A_N20000_cs1000') and r['k'] in ('0', '1', '3'):
        n = int(r['rows']); print('multi', r['label'], 'k', r['k'], 'rows', n, 'arena/row', round((int(r['arena_fn']) + 144 * n) / n, 1))
print('ratio to BayesMVP 1608:', round(fit(comp['total'])[0] / 1608, 2))
