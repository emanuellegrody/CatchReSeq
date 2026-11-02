# Chemistry fingerprint for HIV-seq libraries. HIV-seq RT primers are HIV capture
# sequences appended to the same handle as the 10x poly(dT) RT primer, so in R2 the
# handle is followed by an HIV capture sequence only in Primer libraries; in a
# standard library it is followed by poly(T).
import gzip, sys, collections
HANDLE  = "CAACGCAGAGTAC"
PRIMERS = {"U5_PBS": "GGGCGCCACTGC", "gag": "GCACCCATCTCT", "RRE": "GTCTGGCCTGTA"}
label, path, nreads = sys.argv[1], sys.argv[2], int(sys.argv[3])
n = handle = polyT = 0
hits = collections.Counter(); other = collections.Counter()
with gzip.open(path, "rt") as f:
    for i, line in enumerate(f):
        if i % 4 != 1: continue
        n += 1
        if n > nreads: n -= 1; break
        j = line.find(HANDLE)
        if j < 0: continue
        handle += 1
        after = line[j + len(HANDLE):]
        if after.startswith("TTTTTTTT"): polyT += 1; continue
        name = next((p for p, s in PRIMERS.items() if after.startswith(s)), None)
        if name: hits[name] += 1
        elif len(after) >= 12 and after[:12].count("T") < 8: other[after[:12]] += 1
per_m = lambda x: 1e6 * x / n
top = ";".join(f"{k}:{v}" for k, v in other.most_common(3))
print("\t".join([label, path.split("/")[-1], str(n), str(handle), str(polyT),
                 *[str(hits[p]) for p in PRIMERS], f"{per_m(sum(hits.values())):.1f}", top]), flush=True)
