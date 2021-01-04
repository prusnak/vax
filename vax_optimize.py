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


def optimize_remap_harpel(vir):
    remap = {
        "ATG": "ATG",
        "TTT": "TTC",
        "GTT": "GTG",
        "CTT": "CTG",
        "TTA": "CTG",
        "TTG": "CTG",
        "CCA": "CCC",
        "CTA": "CTG",
        "GTC": "GTG",
        "TCT": "AGC",
        "AGT": "AGC",
        "CAG": "CAG",
        "TGT": "TGC",
        "AAT": "AAC",
        "ACA": "ACC",
        "ACC": "ACC",
        "AGA": "AGA",
        "ACT": "ACC",
        "CAA": "CAG",
        "CCC": "CCC",
        "CCT": "CCT",
        "GCA": "GCC",
        "TAC": "TAC",
        "TTC": "TTC",
        "CGT": "AGA",
        "GGT": "GGC",
        "TAT": "TAC",
        "GAC": "GAC",
        "AAA": "AAG",
        "TCC": "AGC",
        "TCA": "AGC",
        "CAT": "CAC",
        "TGG": "TGG",
        "GCT": "GCC",
        "ATA": "ATC",
        "GGG": "GGA",
        "AAG": "AAG",
        "AGG": "CGG",
        "GAT": "GAC",
        "AAC": "AAC",
        "GAG": "GAG",
        "GGC": "GGC",
        "ATT": "ATC",
        "TCG": "AGC",
        "GAA": "GAG",
        "CAC": "CAC",
        "GCG": "GCC",
        "TGC": "TGC",
        "GGA": "GGC",
        "GTG": "GTG",
        "ACG": "ACC",
        "CTC": "CTG",
        "GTA": "GTG",
        "ATC": "ATC",
        "GCC": "GCC",
        "AGC": "AGC",
        "CTG": "CTG",
        "CGG": "CGG",
        "CGC": "CGG",
        "TAA": "TGA",
    }
    return [remap[c] for c in vir]


class Runner:
    def __init__(self):
        self.vir = open("codons/virus.txt").read().splitlines()
        self.vax = open("codons/tozinameran.txt").read().splitlines()
        # self.vax = open('codons/zorecimeran.txt').read().splitlines()

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

    def average_runs(self, func, iters=20):
        score_codons, score_bases = 0.0, 0.0
        for _ in range(iters):
            v = func(self.vir)
            score_codons += self.compute_match(self.vir, self.vax, v)
            score_bases += self.compute_match(self.vir, self.vax, v, compare_bases=True)
        score_codons /= iters
        score_bases /= iters
        print(f"{func.__name__} : {score_codons:.2%} {score_bases:.2%}")


run = Runner()

# run.average_runs(optimize_dnachisel)

run.average_runs(optimize_remap_pct)

run.average_runs(optimize_remap_harpel)
