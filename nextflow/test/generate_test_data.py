#!/usr/bin/env python3
"""Generate synthetic test data for MarkerMAG Nextflow pipeline testing.

Uses one real 16S rRNA sequence from the bundled SILVA file so that the
barrnap-based POLISH_16S step accepts it.  MAG contigs are random sequence.

Produces:
  data/16s/test_16S.fasta   – 1 real 16S sequence (from SILVA_16S_order.fasta)
  data/synthetic_mags/test_mag.fa – 2 random contigs (3000 + 2000 bp)
  data/test_R1.fasta        – 160 reads named test_N.1
  data/test_R2.fasta        – 160 reads named test_N.2

Naming follows MarkerMAG's pairing convention (PREFIX_INDEX.1 / .2).
60 pairs bridge 16S→MAG, 60 bridge MAG→16S, 20 pure-16S, 20 pure-MAG.
"""
import random
from pathlib import Path

random.seed(42)
BASES    = 'ACGT'
BASE_DIR = Path(__file__).parent
SILVA    = Path(__file__).parent.parent.parent / 'MarkerMAG' / 'SILVA_16S_order.fasta'


def rseq(n):
    return ''.join(random.choice(BASES) for _ in range(n))


def write_fasta(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as fh:
        for name, seq in records:
            fh.write(f'>{name}\n{seq}\n')


def read_first_fasta(path):
    """Return (header, sequence) for the first record in a FASTA file."""
    header, seqlines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header:
                    break
                header = line[1:]
            else:
                seqlines.append(line)
    return header, ''.join(seqlines)


# ── Reference sequences ────────────────────────────────────────────────────────
# 16S: first real sequence from the bundled SILVA file
silva_header, seq_16s = read_first_fasta(SILVA)
seq_16s = seq_16s.upper().replace('U', 'T')   # RNA→DNA if needed

# MAG: two random contigs (entirely distinct from the 16S)
mag_ctg1 = rseq(3000)
mag_ctg2 = rseq(2000)

write_fasta(BASE_DIR / 'data/16s/test_16S.fasta',
            [('test_16S_1', seq_16s)])
write_fasta(BASE_DIR / 'data/synthetic_mags/test_mag.fa',
            [('contig_1', mag_ctg1), ('contig_2', mag_ctg2)])

# ── Read pairs ─────────────────────────────────────────────────────────────────
READ_LEN = 150
r1, r2   = [], []
idx      = 1

# 60 bridging pairs: R1 from 16S body, R2 from MAG contig_1
for _ in range(60):
    p16 = random.randint(0, len(seq_16s) - READ_LEN)
    pm  = random.randint(0, len(mag_ctg1) - READ_LEN)
    r1.append((f'test_{idx}.1', seq_16s[p16:p16 + READ_LEN]))
    r2.append((f'test_{idx}.2', mag_ctg1[pm:pm + READ_LEN]))
    idx += 1

# 60 reverse bridging: R1 from MAG, R2 from 16S
for _ in range(60):
    p16 = random.randint(0, len(seq_16s) - READ_LEN)
    pm  = random.randint(0, len(mag_ctg1) - READ_LEN)
    r1.append((f'test_{idx}.1', mag_ctg1[pm:pm + READ_LEN]))
    r2.append((f'test_{idx}.2', seq_16s[p16:p16 + READ_LEN]))
    idx += 1

# 20 pure-16S pairs (noise)
for _ in range(20):
    p1 = random.randint(0, len(seq_16s) - READ_LEN)
    p2 = random.randint(0, len(seq_16s) - READ_LEN)
    r1.append((f'test_{idx}.1', seq_16s[p1:p1 + READ_LEN]))
    r2.append((f'test_{idx}.2', seq_16s[p2:p2 + READ_LEN]))
    idx += 1

# 20 pure-MAG pairs (noise)
for _ in range(20):
    p1 = random.randint(0, len(mag_ctg2) - READ_LEN)
    p2 = random.randint(0, len(mag_ctg2) - READ_LEN)
    r1.append((f'test_{idx}.1', mag_ctg2[p1:p1 + READ_LEN]))
    r2.append((f'test_{idx}.2', mag_ctg2[p2:p2 + READ_LEN]))
    idx += 1

write_fasta(BASE_DIR / 'data/test_R1.fasta', r1)
write_fasta(BASE_DIR / 'data/test_R2.fasta', r2)

print(f'Generated {idx - 1} read pairs')
print(f'  16S:  data/16s/test_16S.fasta  ({len(seq_16s)} bp, from SILVA: {silva_header[:60]})')
print(f'  MAGs: data/synthetic_mags/test_mag.fa (contig_1: 3000 bp, contig_2: 2000 bp)')
print(f'  R1:   data/test_R1.fasta')
print(f'  R2:   data/test_R2.fasta')
