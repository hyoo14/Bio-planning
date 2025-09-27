

# Parse reactions file and get product2reactant and reactant2product

import csv
from collections import defaultdict

def key_from_side(side: str) -> tuple:
    """'.'로 분리 → 공백 제거 → 정렬 → tuple 키(멀티셋)"""
    parts = [p.strip() for p in side.split(".") if p.strip()]
    parts.sort()
    return tuple(parts)  # 예: ("A","A","B")

def parse_reaction_line(line: str):
    """SMILES 반응식 'A.B>>C' 형태를 (reactants, products)로 분리"""
    if ">>" not in line:
        return None
    left, right = line.split(">>", 1)
    return left.strip(), right.strip()

def build_maps_from_csv(csv_file):
    r2p = defaultdict(list)  # { reactants_key: [(id, product)] }
    p2r = defaultdict(list)  # { products_key: [(id, reactants)] }

    with open(csv_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)  # , delimiter="\t"  # 탭 구분 CSV
        for row in reader:
            rid = row["id"]
            smiles = row["smiles_original"]
            parsed = parse_reaction_line(smiles)
            if not parsed:
                continue
            L, R = parsed
            kL = key_from_side(L)
            kR = key_from_side(R)

            r2p[kL].append((rid, R))
            p2r[kR].append((rid, L))

    return dict(r2p), dict(p2r)

def lookup_products(reactants_str: str, r2p: dict):
    """주어진 reactants로 만들 수 있는 (id, product) 리스트 반환"""
    return r2p.get(key_from_side(reactants_str), [])

def lookup_reactants(products_str: str, p2r: dict):
    """주어진 products로 만들 수 있는 (id, reactants) 리스트 반환"""
    return p2r.get(key_from_side(products_str), [])

csv_file =  OUTPUT_FILE
reactants_to_products, products_to_reactants = build_maps_from_csv(csv_file)



# filtering reaction files

from rdkit import Chem

# SMARTS 패턴 (이미 정의된 것 그대로 사용)
patt_amino = Chem.MolFromSmarts("C(C(=O)O)N")
patt_phosphate = Chem.MolFromSmarts("P(=O)(O)O")
patt_isoprenoid = Chem.MolFromSmarts("C=C-C=C")


def parse_line(line: str):
    line = line.strip()
    if not line or line.startswith("#") or ">>" not in line:
        return None
    left, right = line.split(">>", 1)
    return left.strip(), right.strip()

def has_motif(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return (mol.HasSubstructMatch(patt_amino) )# or
            # mol.HasSubstructMatch(patt_phosphate) or
            # mol.HasSubstructMatch(patt_isoprenoid))

def filter_reaction_file(in_path, out_path):
    kept = []
    with open(in_path, encoding="utf-8") as fin:
        for line in fin:
            parsed = parse_line(line)   # 네가 위에 정의한 parse_line 재사용
            if not parsed:
                continue
            L, R = parsed
            # reactants, products 각각 motif 체크
            keep = False
            for r in L.split("."):
                if has_motif(r.strip()):
                    keep = True
                    break
            if not keep:
                for p in R.split("."):
                    if has_motif(p.strip()):
                        keep = True
                        break
            if keep:
                kept.append(line.strip())

    # 결과 저장
    with open(out_path, "w", encoding="utf-8") as fout:
        for line in kept:
            fout.write(line + "\n")
    print(f"{len(kept)}개 반응 저장 완료 → {out_path}")



# 사용 예시
filter_reaction_file(
    "/content/drive/MyDrive/_September_hard/task_04/test_data_01/train_canonicalized_filtered.txt",
    "/content/drive/MyDrive/_September_hard/task_04/test_data_01/train_canonicalized_filtered_amino_only.txt"
)





# check whether target can reach to the ground truth

from collections import deque
from pathlib import Path

def is_reachable_from_target(target: str, goal: str, products_to_reactants: dict, max_depth: int = 15) -> bool:
    """target에서 시작해서 goal reactant까지 거슬러 올라갈 수 있는지 확인"""
    queue = deque([(target, 0)])
    visited = set()
    while queue:
        current, depth = queue.popleft()
        if depth > max_depth:
            continue
        if current in visited:
            continue
        visited.add(current)

        if current.strip() == goal.strip():  # 최종 reactant에 도달하면 성공
            return True

        parents = products_to_reactants.get(key_from_side(current), [])
        for parent in parents:
            queue.append((parent, depth + 1))

    return False


def check_test_gt(gt_path, products_to_reactants):
    lines = Path(gt_path).read_text(encoding="utf-8").splitlines()
    total = len(lines)
    found = 0

    for i, line in enumerate(lines, 1):
        cols = line.strip().split("\t")
        if len(cols) < 6:
            continue  # target과 goal 둘 다 있어야 함
        target = cols[3]          # 5번째 컬럼이 target
        goal = cols[-1]           # 마지막 컬럼이 최종 reactant

        print(f"length of cols: {len(cols)} target: {target}, goal: {goal}")
        
        ok = is_reachable_from_target(target, goal, products_to_reactants)
        print(f" is it ok? {ok}")
        if ok:
            found += 1
        print(f"[{i}/{total}] target={target[:40]}... goal={goal[:30]} => {'OK' if ok else 'FAIL'}")


      

        

    print(f"\n총 {total}개 중 {found}개 reachable")


# 실행 예시
check_test_gt("/content/drive/MyDrive/_September_hard/task_04/test_data_01/test_gt.txt", products_to_reactants)



# filtering test file(with ground truth)

from pathlib import Path
import pandas as pd

def filter_test_gt_by_check(gt_path, products_to_reactants, out_path="/content/drive/MyDrive/_September_hard/task_04/test_data_01/test_gt_filtered_amino_only.txt"):
    lines = Path(gt_path).read_text(encoding="utf-8").splitlines()
    kept = []
    total = len(lines)
    found = 0

    for i, line in enumerate(lines, 1):
        cols = line.strip().split("\t")
        if len(cols) < 6:
            continue  # target과 goal 둘 다 있어야 함
        target = cols[3]          # 5번째 컬럼 ########   이게 문제였음!! 4번째 칼럼임
        goal = cols[-1]           # 마지막 컬럼
        ok = is_reachable_from_target(target, goal, products_to_reactants)
        if ok:
            kept.append(cols)
            found += 1
        print(f"[{i}/{total}] target={target[:40]}... goal={goal[:30]} => {'OK' if ok else 'FAIL'}")

    # DataFrame으로 저장
    df_filtered = pd.DataFrame(kept)
    df_filtered.to_csv(out_path, sep="\t", index=False, header=False)
    print(f"\n총 {total}개 중 {found}개 reachable → {out_path}")
    return df_filtered
