# JRP
from collections import defaultdict
from typing import Dict, List, Sequence, Set, Tuple, Optional, Iterable
import math
import numpy as np
import time

starttime = time.time()
maxMultiplicity = 3  # largest considered multiplicity (maxMultiplicity is included) 
#   3 -- runtime ~ 1.7min          4 -- runtime ~  40min               5 -- runtime ~ 1.5-5h             6 -- runtime ~ 12-48h
# for a 114Cd run, for two clusters
    # cluster of size 1 contributes  24 million counts (~74%)
    # cluster of size 2 contributes 8 million counts (~24%)
    # cluster of size 3 contributes 740k counts (~2%)
    # cluster of size 4 contributes 93k counts (~0.3%)
    # cluster of size 5 contributes 10k counts (~0.03%)
    # cluster of size 6 contributes 2500 counts (~0.007%)

NDET = 162
cluster_verbose = True
cluster_pair_verbose = False
angle_precision = 1
nearestNeighbors = {}
clusterAngles = []
def from_file_idx(x: int) -> int:
    return -1 if x == 200 else x  # keep 0..161 as-is

theta= [90.0 , 99.6795 , 105.786 , 90.0 , 74.2137 , 80.3205 , 90.0 , 108.001 , 117.284 , 120.001 , 106.458 , 90 , 73.5423 , 59.9987 , 62.716 , 71.9993 , 97.41 , 115.379 , 125.971 , 131.842 , 133.908 , 119.471 , 106.458 , 90 , 73.5423 , 60.529 , 46.0923 , 48.1581 , 54.0294 , 64.6209 , 82.59 , 89.9979 , 105.097 , 121.719 , 134.479 , 144.002 , 149.499 , 148.287 , 133.907 , 120.003 , 105.786 , 90 , 74.2137 , 60.0005 , 46.0934 , 31.7134 , 30.501 , 35.9954 , 45.5211 , 58.281 , 74.9026 , 97.4085 , 115.378 , 134.478 , 150.535 , 161.868 , 164.908 , 149.496 , 131.839 , 117.284 , 99.6795 , 80.3205 , 62.716 , 48.1612 , 30.5037 , 15.0917 , 18.1315 , 29.4654 , 45.5223 , 64.6224 , 82.5915 , 90 , 107.998 , 125.968 , 143.998 , 161.864 , 180 , 161.864 , 144.002 , 125.968 , 108.002 , 90 , 72.0019 , 54.0322 , 36.0019 , 18.1361 , 0 , 18.1361 , 35.9981 , 54.0322 , 71.9981 , 99.6795 , 117.284 , 131.839 , 149.496 , 164.908 , 161.868 , 150.535 , 134.478 , 115.378 , 97.4085 , 82.5915 , 64.6224 , 45.5223 , 29.4654 , 18.1315 , 15.0918 , 30.5037 , 48.1612 , 62.716 , 80.3205 , 90 , 105.786 , 120 , 133.907 , 148.287 , 149.499 , 144.005 , 134.479 , 121.719 , 105.097 , 90.0021 , 74.9026 , 58.281 , 45.5211 , 35.9976 , 30.501 , 31.7134 , 46.0935 , 59.9966 , 74.2137 , 90 , 106.458 , 119.471 , 133.908 , 131.842 , 125.971 , 115.379 , 97.41 , 82.59 , 64.6209 , 54.0294 , 48.1581 , 46.0923 , 60.529 , 73.5423 , 90 , 106.458 , 120.001 , 117.284 , 108.001 , 90 , 71.9993 , 62.716 , 59.9987 , 73.5423 , 90 , 105.786 , 99.6795 , 80.3205 , 74.2137 , 90]
phi= [0.0 , 346.422 , 5.27057 , 16.6216 , 5.27057 , 346.422 , 331.184 , 333.434 , 350.352 , 10.8129 , 23.9913 , 31.7189 , 23.9913 , 10.8129 , 350.352 , 333.434 , 318.313 , 319.236 , 336.211 , 353.747 , 18.2257 , 31.7205 , 39.4487 , 46.8184 , 39.4487 , 31.7205 , 18.2257 , 353.747 , 336.211 , 319.236 , 318.313 , 301.713 , 301.713 , 301.712 , 315.341 , 333.428 , 359.312 , 31.7231 , 45.2171 , 52.6262 , 58.1694 , 63.44 , 58.1694 , 52.6279 , 45.2171 , 31.7231 , 359.312 , 333.434 , 315.341 , 301.712 , 301.713 , 285.113 , 284.188 , 288.081 , 301.709 , 326.181 , 31.7284 , 64.1331 , 69.6945 , 73.0877 , 77.0176 , 77.0176 , 73.0877 , 69.6945 , 64.1331 , 31.7284 , 326.181 , 301.709 , 288.081 , 284.188 , 285.113 , 272.256 , 270 , 267.213 , 270 , 277.264 , 0 , 97.2639 , 90 , 87.2127 , 90 , 92.2556 , 90 , 87.2126 , 90 , 97.2639 , 0 , 277.264 , 270 , 267.213 , 270 , 257.018 , 253.088 , 249.694 , 244.133 , 211.728 , 146.181 , 121.709 , 108.081 , 104.188 , 105.113 , 105.113 , 104.188 , 108.081 , 121.709 , 146.181 , 211.728 , 244.133 , 249.694 , 253.088 , 257.018 , 243.44 , 238.169 , 232.628 , 225.217 , 211.723 , 179.312 , 153.434 , 135.341 , 121.712 , 121.713 , 121.713 , 121.713 , 121.712 , 135.341 , 153.428 , 179.312 , 211.723 , 225.217 , 232.626 , 238.169 , 226.818 , 219.449 , 211.72 , 198.226 , 173.747 , 156.211 , 139.236 , 138.313 , 138.313 , 139.236 , 156.211 , 173.747 , 198.226 , 211.72 , 219.449 , 211.719 , 203.991 , 190.813 , 170.352 , 153.434 , 151.184 , 153.434 , 170.352 , 190.813 , 203.991 , 196.622 , 185.271 , 166.422 , 166.422 , 185.271 , 180.0]


def all_clusters_list(clusters_by_size, min_size=1, max_size=None):
    sizes = sorted(k for k in clusters_by_size.keys()
                   if k >= min_size and (max_size is None or k <= max_size))
    flat = []
    for k in sizes:
        flat.extend(sorted(clusters_by_size[k]))  # clusters are tuples
    return flat


def build_angle_lookupTable(theta, phi):
    lut = np.zeros((NDET, NDET))
    for i in range(NDET):
        for j in range(NDET):
            if i == j:
                lut[i, j] = 0.0
                continue
            t1 = math.radians(theta[i])
            t2 = math.radians(theta[j])
            dphi = math.radians(abs(phi[i] - phi[j]))
            cosang = math.cos(t1)*math.cos(t2) + math.sin(t1)*math.sin(t2)*math.cos(dphi)
            cosang = max(-1.0, min(1.0, cosang))  # clamp for safety
            lut[i, j] = math.degrees(math.acos(cosang))
    return lut

angle_lut = build_angle_lookupTable(theta, phi)


def GetAngle(det1, det2, theta, phi):
    theta_1_rad = math.radians(theta[det1])
    theta_2_rad = math.radians(theta[det2])
    delta_phi_rad = math.radians(abs(phi[det1] - phi[det2]))

    cos_angle_btwn = (
        math.cos(theta_1_rad) * math.cos(theta_2_rad)
        + math.sin(theta_1_rad) * math.sin(theta_2_rad) * math.cos(delta_phi_rad)
    )

    # Clamp for numerical safety
    cos_angle_btwn = max(-1.0, min(1.0, cos_angle_btwn))
    angle_btwn_deg = math.degrees(math.acos(cos_angle_btwn))
    return angle_btwn_deg



with open("cluster_angle_counts_verbose.txt", "w") as output:
    
    with open("/home/jr2514/DANCE/DANCE_Analysis/Config/DetectorMatrix.txt", "r") as f:
        for line in f:
            if not line.strip():
                continue
            seed_det, nn1, nn2, nn3, nn4, nn5, nn6 = map(int, line.split())
            seed = from_file_idx(seed_det)
            nns  = [from_file_idx(x) for x in (nn1, nn2, nn3, nn4, nn5, nn6)]
            nns  = sorted({n for n in nns if 0 <= n < NDET})  # drop sentinel, dedupe
            if 0 <= seed < NDET:
                nearestNeighbors[seed] = nns

    # live mask (adjust as needed)
    live = [True]*NDET
    for dead in (76, 86):
        live[dead] = False

    # build adjacency (filter out dead neighbors)
    adj = []
    for i in range(NDET):
        neighbors = [n for n in nearestNeighbors.get(i, []) if live[n]]
        # de-dup + sort for stability
        adj.append(sorted(set(neighbors)))

    # enforce symmetry
    # sanity: enforce symmetry one time after you build `adj`
    for i in range(NDET):
        for j in adj[i]:
            if i not in adj[j]:
                adj[j] = sorted(set(adj[j] + [i]))

    '''# --- sanity prints to verify against your file snippet ---
    print("\n--- Building adjacency list (adj) ---",file=output)
    for i in range(0, 13):  # print first 13 to match your snippet
        print(f"Detector {i}: neighbors -> {adj[i]}",file=output)'''

        # ---------- Helpers now use adj (not nearestNeighbors directly) ----------

    #print(adj)

    def enumerate_sizek_for_seed(seed: int,
                             k: int,
                             adj: List[List[int]],
                             live: Sequence[bool],
                             cluster_verbose:bool) -> Set[Tuple[int, ...]]:
        """
        Enumerate all unique connected clusters of EXACTLY size k that contain `seed`.
        Uses functional recursion (no in-place mutation) so prints reflect the true branch.
        """
        out: Set[Tuple[int, ...]] = set()
        if not (0 <= seed < NDET) or not live[seed] or k < 1:
            return out
    
        if cluster_verbose:
            print(f"\nEnumerating size-{k} clusters for seed {seed}", file=output)
    
        def frontier(current_cluster: List[int]) -> List[int]:
            inC = set(current_cluster)
            f: Set[int] = set()
            for d in current_cluster:
                for nb in adj[d]:
                    if live[nb] and nb not in inC:
                        f.add(nb)
            f_list = sorted(f)
            if cluster_verbose:
                print(f"  Frontier for {current_cluster}: {f_list}", file=output)
            return f_list
    
        def grow(cluster_so_far: List[int], depth: int):
            size = len(cluster_so_far)
            if size == k:
                cl = tuple(cluster_so_far)  # already sorted
                if cl not in out:
                    out.add(cl)
                    if cluster_verbose:
                        print(f"{'  ' * depth}Cluster (size {size}): {list(cl)}", file=output)
                return
    
            f = frontier(cluster_so_far)
            for nb in f:
                next_cluster = sorted(cluster_so_far + [nb])  # <-- NO IN-PLACE MUTATION
                #if cluster_verbose:
                    #print(f"{'  ' * depth}→ Add {nb} to {cluster_so_far} → {next_cluster}", file=output)
                grow(next_cluster, depth + 1)
    
        grow([seed], depth=0)
        return out


    def enumerate_sizes_up_to_k_for_seed(seed: int,
                                        k: int,
                                        adj: List[List[int]],
                                        live: Sequence[bool],
                                        cluster_verbose:bool) -> Dict[int, Set[Tuple[int, ...]]]:
        """
        Convenience: enumerate all connected clusters of sizes 1..k (each contains `seed`).
        Uses the same engine, calling enumerate_sizek_for_seed for each size.
        """
        results: Dict[int, Set[Tuple[int, ...]]] = {}
        for sz in range(1, k + 1):
            results[sz] = enumerate_sizek_for_seed(seed, sz, adj, live, cluster_verbose)
        return results
        

        

    def enumerate_whole_ball(maxMultiplicity: int,
                            adj: List[List[int]],
                            live: Sequence[bool]) -> Dict[int, Set[Tuple[int, ...]]]:
        merged: Dict[int, Set[Tuple[int, ...]]] = defaultdict(set)
        for seed in range(NDET):
            if not live[seed]:
                continue
            for size in range(1, maxMultiplicity + 1):
                clusters = enumerate_sizek_for_seed(seed, size, adj, live, cluster_verbose)
                merged[size].update(clusters)
        return merged


    allclusters = []
    for seed in range(162):
        for k in range(1, maxMultiplicity+1):
            clusters = enumerate_sizek_for_seed(seed, k, adj, live, cluster_verbose)
            allclusters.append(clusters)
            if cluster_verbose:
                print(f"\nTotal clusters for seed {seed}, size {k}: {len(clusters)}", file=output)
    # Whole-ball enumeration for sizes 1..maxMultiplicity (quiet)
    merged = enumerate_whole_ball(maxMultiplicity, adj, live)
    for sz in range(1, maxMultiplicity+1):
        if cluster_verbose:
            print(f"WHOLE BALL: unique clusters size {sz} = {len(merged.get(sz, set()))}", file=output)






    Cluster = Tuple[int, ...]  # sorted tuple of detector IDs

    def clusters_are_disjoint(c1: Cluster, c2: Cluster) -> bool:
        return not (set(c1) & set(c2))

    def clusters_are_non_touching(c1: Cluster, c2: Cluster, adj: List[List[int]]) -> bool:
        # No edge between any member of c1 and any member of c2
        s2 = set(c2)
        for a in c1:
            if s2 & set(adj[a]):
                return False
        return True

    def all_disjoint_pairs_two_loops(
            clusters: List[Cluster],
            *,
            adj: Optional[List[List[int]]] = None,
            require_non_touching: bool = False,
            output=None,
            echo: bool = False,
        ) -> List[Tuple[Cluster, Cluster]]: 
        """
        Loop over all unordered pairs (i<j). Keep pairs that are disjoint.
        If require_non_touching=True, also require no NN edges between clusters.
        """
        pairs: List[Tuple[Cluster, Cluster]] = []
        for i, c1 in enumerate(clusters):
            for j in range(i + 1, len(clusters)):
                c2 = clusters[j]
                if not clusters_are_disjoint(c1, c2):
                    continue
                if require_non_touching:
                    assert adj is not None, "adj required for non-touching check"
                    if not clusters_are_non_touching(c1, c2, adj):
                        continue
                pairs.append((c1, c2))
                if echo and output is not None:
                    rule = "non-touching" if require_non_touching else "disjoint"
                    #print(f"{rule} pair: {list(c1)} | {list(c2)}", file=output)
        return pairs

    
    clusters_by_size = enumerate_whole_ball(maxMultiplicity, adj, live)
    all_clusters = all_clusters_list(clusters_by_size)  # or all_clusters_set(...)
    #print(f"Total clusters (sizes 1..{maxMultiplicity}): {len(all_clusters)}", file=output)
    #print(all_clusters)

    pairs_non_touching = all_disjoint_pairs_two_loops(all_clusters, adj=adj, require_non_touching=True, output=output, echo=False)
    print(f"Total disjoint & non-touching pairs: {len(pairs_non_touching)}")

    if cluster_pair_verbose:
        for idx in range(10):
            print(f"pairs_non_touching[{idx}]: {pairs_non_touching[idx]}   cluster1", pairs_non_touching[idx][0], "   cluster2", pairs_non_touching[idx][1], file=output)
   

    def accumulate_clusterpair_angle_counts(
                                                                        cluster_pairs: Iterable[Tuple[Cluster, Cluster]],
                                                                        angle_lut,                    # shape (NDET, NDET), degrees
                                                                        angle_precision,
                                                                        unique_per_pair: bool = False,
                                                                        ) -> Dict[int, int]:
        """
        Build angle → count map over cluster pairs.
        angle key = int(round(angle_deg * scale)), where scale = 10**angle_precision
        If unique_per_pair=True: each angle bin counted at most once for each cluster pair.
        """
        counts: Dict[int, int] = defaultdict(int)
        scale = 10**angle_precision
        for (c1, c2) in cluster_pairs:
            seen_angles: Set[int] = set() if unique_per_pair else None
            # cross all detectors in the two clusters
            for a in c1:
                for b in c2:
                    # (clusters are disjoint / non-touching already; a!=b guaranteed)
                    ang = round(angle_lut[a][b], angle_precision)
                    if unique_per_pair:
                        if ang in seen_angles:
                            continue
                        seen_angles.add(ang)
                    counts[ang] += 1
        return counts

endtime = time.time()
print(f"Analysis initialization elapsed time {endtime - starttime:.3f} seconds")
starttime = time.time()
angle_counts = accumulate_clusterpair_angle_counts(pairs_non_touching, angle_lut, angle_precision)

with open(f"cluster_angle_counts_output_Mcr1to{maxMultiplicity}.txt", "w") as clOut:
    for ang in sorted(angle_counts):
        print(f"{ang:.1f} deg : {angle_counts[ang]}", file=clOut)
endtime = time.time()
print(f"Getting angle counts (for cluster size 1 - {maxMultiplicity}) elapsed time {endtime - starttime:.3f} seconds")