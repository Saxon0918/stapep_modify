import os
from pymol import cmd

# 设置输入输出目录
pdb_folder = "/home/d3008/Documents/zhr/1OU8/1OU8_100_1/RFdiffusion1/"
output_folder = "/home/d3008/Documents/zhr/1OU8/1OU8_100_1/RFdiffusion/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

cmd.reinitialize()

for filename in os.listdir(pdb_folder):
    if filename.endswith(".pdb"):
        filepath = os.path.join(pdb_folder, filename)
        name = os.path.splitext(filename)[0]

        # 载入结构
        cmd.load(filepath, name)

        # =========================
        # Step 1. 用临时对象分别提取三部分
        # part1: 原A链1-109 -> 将来并入新B链前段
        # part2: 原B链 -> 将来并入新B链后段
        # part3: 原A链110-end -> 将来作为新A链
        # =========================
        part1 = f"{name}_part1"
        part2 = f"{name}_part2"
        part3 = f"{name}_part3"
        new_obj = f"{name}_reordered"

        cmd.create(part1, f"{name} and chain A and resi 1-109")
        cmd.create(part2, f"{name} and chain B")
        cmd.create(part3, f"{name} and chain A and not resi 1-109")

        # =========================
        # Step 2. 统一链名
        # part1 和 part2 都设为 B
        # part3 设为 A
        # =========================
        cmd.alter(part1, "chain='B'")
        cmd.alter(part2, "chain='B'")
        cmd.alter(part3, "chain='A'")
        cmd.sort()

        # =========================
        # Step 3. 按顺序拼接
        # 顺序必须是：
        # 原A 1-109 -> 原B -> 原A剩余
        # =========================
        cmd.create(new_obj, f"({part1}) or ({part2}) or ({part3})")
        cmd.sort()

        # =========================
        # Step 4. 重新编号新的B链
        # B链顺序：part1 在前，part2 在后
        # 需要整体连续编号
        # =========================
        stored_b = {"atoms": [], "counter": 1}

        # 先收集 B链中每个残基，按出现顺序记录
        cmd.iterate(
            f"{new_obj} and chain B and name CA",
            "stored_b['atoms'].append(resi)",
            space={"stored_b": stored_b}
        )

        # 去重并保持顺序
        old_resi_order_b = []
        for r in stored_b["atoms"]:
            if r not in old_resi_order_b:
                old_resi_order_b.append(r)

        # 由于 part1 和 part2 里可能都有相同 resi 号，
        # 直接按 resi 改会冲突，所以先分两段分别处理
        # 先给 part1 和 part2 打临时 segi 标记
        cmd.alter(f"{part1}", "segi='P1'")
        cmd.alter(f"{part2}", "segi='P2'")
        cmd.alter(f"{part3}", "segi='P3'")

        cmd.create(new_obj, f"({part1}) or ({part2}) or ({part3})")
        cmd.sort()

        # 重新给 B链连续编号：先 P1，再 P2
        counter = 1

        stored_p1 = {"resis": []}
        cmd.iterate(
            f"{new_obj} and segi P1 and name CA",
            "stored_p1['resis'].append(resi)",
            space={"stored_p1": stored_p1}
        )
        p1_resis = []
        for r in stored_p1["resis"]:
            if r not in p1_resis:
                p1_resis.append(r)

        for old_resi in p1_resis:
            cmd.alter(
                f"{new_obj} and segi P1 and resi {old_resi}",
                f"resi='{counter}'"
            )
            counter += 1

        stored_p2 = {"resis": []}
        cmd.iterate(
            f"{new_obj} and segi P2 and name CA",
            "stored_p2['resis'].append(resi)",
            space={"stored_p2": stored_p2}
        )
        p2_resis = []
        for r in stored_p2["resis"]:
            if r not in p2_resis:
                p2_resis.append(r)

        for old_resi in p2_resis:
            cmd.alter(
                f"{new_obj} and segi P2 and resi {old_resi}",
                f"resi='{counter}'"
            )
            counter += 1

        # =========================
        # Step 5. 重新编号新的A链（原A链110+）
        # 从1开始重新编号
        # =========================
        stored_p3 = {"resis": [], "counter": 1}
        cmd.iterate(
            f"{new_obj} and segi P3 and name CA",
            "stored_p3['resis'].append(resi)",
            space={"stored_p3": stored_p3}
        )

        p3_resis = []
        for r in stored_p3["resis"]:
            if r not in p3_resis:
                p3_resis.append(r)

        for i, old_resi in enumerate(p3_resis, start=1):
            cmd.alter(
                f"{new_obj} and segi P3 and resi {old_resi}",
                f"resi='{i}'"
            )

        # 去掉临时 segi
        cmd.alter(new_obj, "segi=''")
        cmd.sort()

        # 保存
        output_path = os.path.join(output_folder, f"{name}.pdb")
        cmd.save(output_path, new_obj)
        print(f"Processed and saved: {output_path}")

        # 清理对象
        cmd.delete(name)
        cmd.delete(part1)
        cmd.delete(part2)
        cmd.delete(part3)
        cmd.delete(new_obj)

cmd.quit()
