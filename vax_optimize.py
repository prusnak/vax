from dnachisel import *
from dnachisel import biotools
import python_codon_tables as pct

species = "h_sapiens"
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


def optimize_remap(vir):
    ct = pct.get_codons_table(species)
    r = []
    for a in vir:
        b = biotools.translate(a)
        t = [(k, v) for k, v in sorted(ct[b].items(), key=lambda i: i[1], reverse=True)]
        r.append(t[0][0])
    return r


class Runner:
    def __init__(self):
        self.vir = open("codons_virus.txt").read().splitlines()
        self.vax = open("codons_tozinameran.txt").read().splitlines()
        # self.vax = open('codons_zorecimeran.txt').read().splitlines()

    def compute_match(self, seq0, seq1, seq2):
        if len(seq0) != len(seq1) or len(seq1) != len(seq2):
            raise ValueError("length mismatch")
        good = 0
        # compare whole codons (False) or individual bases (True)
        if False:
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

    def average_runs(self, func, iters=10):
        score = (
            sum(
                [
                    self.compute_match(self.vir, self.vax, func(self.vir))
                    for _ in range(iters)
                ]
            )
            / iters
        )
        print(f"{func.__name__} : {score:.2%}")
        return score


run = Runner()

run.average_runs(optimize_dnachisel, 5)

run.average_runs(optimize_remap, 1)
