vir = open('codons_virus.txt').read().splitlines()
vax_tozinameran = open('codons_tozinameran.txt').read().splitlines()
vax_zorecimeran = open('codons_zorecimeran.txt').read().splitlines()

assert len(vir) == len(vax_tozinameran)
assert len(vax_tozinameran) == len(vax_zorecimeran)

diff = 0
same = 0
total = len(vir)

for i in range(total):
    a, b = vax_tozinameran[i], vax_zorecimeran[i]
    if a == b:
      d = "SAME"
      same += 1
    else:
      d = "DIFF"
      diff += 1
    print(i * 3, a, b, d)

print("DIFF: %d %0.2f%%" % (diff, 100.0 * diff / total))
print("SAME: %d %0.2f%%" % (same, 100.0 * same / total))
print("TOTAL: %d" % total)
