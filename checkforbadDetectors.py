# JRP — find common detectors within angles and across angles
from collections import defaultdict, Counter

FNAME = "angle_combinations.txt"

# angles of interest; we'll normalize/round angles in the file to 1 decimal
BAD_ANGLES = [63.4, 87.7, 90.0, 92.3, 116.6, 180.0]
ANGLE_DP = 1  # decimal places for matching floating angles (tune if needed)

def norm_angle(a, dp=ANGLE_DP):
    return round(float(a), dp)

# --- load file and group by angle ---
angle_to_pairs = defaultdict(list)            # {angle: [(d1,d2), ...]}
pair_seen = set()                             # to avoid counting exact duplicate lines

with open(FNAME) as f:
    for line in f:
        parts = line.split()
        if len(parts) != 3:
            continue
        ang = norm_angle(parts[0])
        d1, d2 = int(parts[1]), int(parts[2])

        # treat pairs as unordered (119,144) == (144,119)
        pair = (min(d1, d2), max(d1, d2), ang)
        if pair in pair_seen:
            continue
        pair_seen.add(pair)

        angle_to_pairs[ang].append((d1, d2))

# --- per-angle repeats + cross-angle analysis ---
bad_norm = [norm_angle(a) for a in BAD_ANGLES]

det_count_by_angle = {}                       # {angle: Counter({det: count, ...})}
det_to_angles = defaultdict(set)              # {det: {angle1, angle2, ...}}
pair_to_angles = defaultdict(set)             # {(d1,d2): {angle1, ...}}

for ang in bad_norm:
    pairs = angle_to_pairs.get(ang, [])
    c = Counter()
    for d1, d2 in pairs:
        c.update([d1, d2])
        det_to_angles[d1].add(ang)
        det_to_angles[d2].add(ang)
        pair_to_angles[(min(d1,d2), max(d1,d2))].add(ang)
    det_count_by_angle[ang] = c

# --- report ---
for ang in bad_norm:
    pairs = angle_to_pairs.get(ang, [])
    if not pairs:
        continue
    print(f"\nAngle {ang:.1f}° → {len(pairs)} pairs")
    for d1, d2 in pairs:
        print(f"  {d1:3d} – {d2:3d}")
    repeats = {d: n for d, n in det_count_by_angle[ang].items() if n > 1}
    if repeats:
        print("  Detectors repeated at this angle:", dict(sorted(repeats.items())))
    else:
        print("  No detector repeats at this angle.")

# detectors showing up at 2+ distinct *bad* angles
multi_angle_dets = {d: sorted(list(angs)) for d, angs in det_to_angles.items() if len(angs) >= 2}
if multi_angle_dets:
    print("\nDetectors appearing across multiple bad angles:")
    for d, angs in sorted(multi_angle_dets.items()):
        print(f"  Det {d:3d} in angles: {', '.join(f'{a:.1f}°' for a in angs)}")

# top detectors by total occurrences across all bad angles
totals = Counter()
for ang in bad_norm:
    totals.update(det_count_by_angle.get(ang, {}))
top = totals.most_common(10)
if top:
    print("\nTop detectors across all bad angles (by appearances):")
    for det, n in top:
        print(f"  Det {det:3d}: {n} appearances")

# optional: pairs that show up at multiple bad angles (can hint at geometry culprits)
pair_multi = {pair: sorted(list(angs)) for pair, angs in pair_to_angles.items() if len(angs) >= 2}
if pair_multi:
    print("\nDetector pairs repeated across multiple bad angles:")
    for (d1, d2), angs in sorted(pair_multi.items()):
        print(f"  {d1:3d}–{d2:3d} in angles: {', '.join(f'{a:.1f}°' for a in angs)}")
