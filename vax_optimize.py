from dnachisel import *
from dnachisel import biotools
import python_codon_tables as pct

species = "9606"  # Homo sapiens
# species = "10090"  # Mus musculus
# species = "57486"  # Mus musculus molossinus
# species = "9544"   # Macaca mulatta


def sequence_codons(seq):
    assert len(seq) % 3 == 0
    for i in range(0, len(seq), 3):
        yield seq[i : i + 3]


def optimize_dnachisel(vir):
    vir = "".join(vir)
    problem = DnaOptimizationProblem(
        sequence=vir,
        constraints=[
            EnforceGCContent(mini=0.54, maxi=0.93, window=120),
            EnforceTranslation(),
            # UniquifyAllKmers(k=12),
            # AvoidHairpins(stem_size=20, hairpin_window=1000),
        ],
        objectives=[
            MaximizeCAI(species=species),
        ],
    )

    # problem.mutations_per_iterations = 10000
    problem.max_random_iters = 20000
    # problem.optimization_stagnation_tolerance = 20000

    problem.resolve_constraints()

    problem.optimize()

    return list(sequence_codons(problem.sequence))


def optimize_remap_pct(vir):
    ct = pct.get_codons_table(species)
    r = []
    for a in vir:
        b = biotools.translate(a)
        t = [(k, v) for k, v in sorted(ct[b].items(), key=lambda i: i[1], reverse=True)]
        r.append(t[0][0])
    return r


def optimize_remap_remap3(vir):
    from remap3 import remap

    return [remap[c] for c in vir]


def optimize_remap_remap6(vir):
    from remap6 import remap

    mvir = vir + ["XXX"]
    vax = []
    for i in range(len(vir)):
        vax.append(remap[mvir[i] + mvir[i + 1]])
    return vax


class Runner:
    def compute_match(self, seq0, seq1, seq2, compare_bases=False):
        if len(seq0) != len(seq1) or len(seq1) != len(seq2):
            raise ValueError("length mismatch")
        good = 0
        if compare_bases:  # False = compare_codons
            seq0 = "".join(seq0)
            seq1 = "".join(seq1)
            seq2 = "".join(seq2)
        for a, b, c in zip(seq0, seq1, seq2):
            if b != c:
                # print(a, b, c, "BAD")
                pass
            else:
                good += 1
                # print(a, b, c, "OK")
                pass
        return good / len(seq1)

    def average_runs(self, virus, vaccine, func, iters=20):
        vir = open(f"codons/{virus}.txt").read().splitlines()
        vax = open(f"codons/{vaccine}.txt").read().splitlines()

        score_codons, score_bases = 0.0, 0.0
        for _ in range(iters):
            v = func(vir)
            score_codons += self.compute_match(vir, vax, v)
            score_bases += self.compute_match(vir, vax, v, compare_bases=True)
        score_codons /= iters
        score_bases /= iters
        print(
            f"{func.__name__}({virus}->{vaccine}): codons={score_codons:.2%} bases={score_bases:.2%}"
        )


run = Runner()

vir = "virus"
for vax in ["tozinameran", "zorecimeran"]:
    run.average_runs(vir, vax, optimize_dnachisel)
    run.average_runs(vir, vax, optimize_remap_pct)
    run.average_runs(vir, vax, optimize_remap_remap3)
    run.average_runs(vir, vax, optimize_remap_remap6)
