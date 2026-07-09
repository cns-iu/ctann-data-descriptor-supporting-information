from __future__ import annotations

import csv
import textwrap
from collections import Counter, defaultdict, deque
from datetime import datetime
from pathlib import Path

import matplotlib


matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "DejaVu Sans"

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "ctann-v7-3.csv"
RUN_ID = datetime.now().strftime("%Y%m%d-%H%M%S")
OUTPUT_DIR = ROOT / "outputs" / "ctaa-v7-curve-spread2x-circle" / RUN_ID
SVG_OUT = OUTPUT_DIR / "ctaa-v7-curve-spread2x-circle.svg"
PNG_OUT = OUTPUT_DIR / "ctaa-v7-curve-spread2x-circle.png"

KEEP_SOURCES = {"azimuth", "frmatch", "ribca"}

COMBO_COLORS = {
    ("azimuth",): "#e53935",
    ("frmatch",): "#1e88e5",
    ("ribca",): "#fdd835",
    ("azimuth", "ribca"): "#fb8c00",
    ("azimuth", "frmatch"): "#8e24aa",
    ("frmatch", "ribca"): "#43a047",
    ("azimuth", "frmatch", "ribca"): "#6d4c41",
}


def normalize_combo(sources: set[str]) -> tuple[str, ...]:
    return tuple(sorted(sources))


def wrap_label(label: str, width: int = 24) -> list[str]:
    lines = textwrap.wrap(label, width=width, break_long_words=False, break_on_hyphens=False)
    return lines or [label]


def main() -> None:
    rows = list(csv.DictReader(INPUT.open(newline="", encoding="utf-8-sig")))

    node_sources: dict[str, set[str]] = defaultdict(set)
    node_label: dict[str, str] = {}
    parents: dict[str, set[str]] = defaultdict(set)
    children: dict[str, set[str]] = defaultdict(set)
    edge_counts: Counter[tuple[str, str]] = Counter()
    row_count = 0

    for row in rows:
        source = row.get("CT/1 - Sources", "").strip().lower()
        if source not in KEEP_SOURCES:
            continue
        row_count += 1

        path: list[tuple[str, str]] = []
        for idx in range(1, 13):
            node_id = row.get(f"AS/{idx}/ID", "").strip()
            node_lab = row.get(f"AS/{idx}/LABEL", "").strip()
            if node_id:
                path.append((node_id, node_lab))

        prev_id: str | None = None
        for node_id, node_lab in path:
            node_sources[node_id].add(source)
            node_label[node_id] = node_lab
            if prev_id and prev_id != node_id:
                parents[node_id].add(prev_id)
                children[prev_id].add(node_id)
                edge_counts[(prev_id, node_id)] += 1
            prev_id = node_id

    indeg = {node: len(parents[node]) for node in node_sources}
    queue = deque([node for node, deg in indeg.items() if deg == 0])
    depth: dict[str, int] = {node: 0 for node in queue}

    while queue:
        node = queue.popleft()
        for child in children[node]:
            indeg[child] -= 1
            if indeg[child] == 0:
                depth[child] = max(depth.get(child, 0), depth[node] + 1)
                queue.append(child)

    changed = True
    while changed:
        changed = False
        for node in node_sources:
            if parents[node]:
                best = max((depth[p] + 1 for p in parents[node] if p in depth), default=0)
                if best > depth.get(node, -1):
                    depth[node] = best
                    changed = True

    ordered: list[str] = []
    seen: set[str] = set()

    def dfs(node: str) -> None:
        if node in seen:
            return
        seen.add(node)
        ordered.append(node)
        for child in sorted(children[node], key=lambda c: (-edge_counts[(node, c)], node_label.get(c, ""), c)):
            dfs(child)

    roots = [node for node, deg in indeg.items() if deg == 0]
    if not roots:
        roots = [min(node_sources)]
    for root in sorted(roots, key=lambda n: (depth.get(n, 0), node_label.get(n, ""), n)):
        dfs(root)
    for node in node_sources:
        if node not in seen:
            dfs(node)

    layers: dict[int, list[str]] = defaultdict(list)
    for node in node_sources:
        layers[depth.get(node, 0)].append(node)
    for d in layers:
        layers[d].sort(key=lambda n: (ordered.index(n), node_label.get(n, ""), n))

    max_depth = max(depth.values()) if depth else 0
    layer_size = max(len(nodes) for nodes in layers.values()) if layers else 1

    left_margin = 50
    top_margin = 130
    right_legend_space = 640
    column_dx = 600
    row_dy = 96
    node_d = 28
    label_gap = 10
    label_line_height = 11
    width = left_margin + max_depth * column_dx + 320 + right_legend_space
    height = top_margin + max(1, layer_size - 1) * row_dy + 180

    positions: dict[str, tuple[float, float]] = {}
    for depth_level, nodes in layers.items():
        total_height = (len(nodes) - 1) * row_dy
        start_y = top_margin + (max(layer_size - 1, len(nodes) - 1) - total_height / row_dy) * 0.5 * row_dy
        if len(nodes) == layer_size:
            start_y = top_margin
        for idx, node in enumerate(nodes):
            positions[node] = (left_margin + depth_level * column_dx, start_y + idx * row_dy)

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, FancyArrowPatch

    fig = plt.figure(figsize=(width / 120, height / 120), dpi=120)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")

    for parent, child in edge_counts:
        if parent not in positions or child not in positions:
            continue
        px, py = positions[parent]
        cx, cy = positions[child]
        start = (px + node_d / 2, py + node_d / 2)
        end = (cx, cy + node_d / 2)
        span = max(40.0, abs(end[0] - start[0]))
        direction = 1 if end[1] >= start[1] else -1
        rad = min(0.35, 0.08 + span / 2000.0) * direction
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-",
                connectionstyle=f"arc3,rad={rad}",
                mutation_scale=10,
                linewidth=0.9 if edge_counts[(parent, child)] < 5 else 1.25,
                color="#4a4a4a",
                alpha=0.42 if edge_counts[(parent, child)] < 5 else 0.62,
                zorder=1,
            )
        )

    for node, (x, y) in positions.items():
        fill = COMBO_COLORS.get(normalize_combo(node_sources[node]), "#bdbdbd")
        ax.add_patch(
            Circle(
                (x + node_d / 2, y + node_d / 2),
                radius=node_d / 2,
                linewidth=1.0,
                edgecolor="#202124",
                facecolor=fill,
                zorder=2,
            )
        )

        label = node_label.get(node, node)
        lines = wrap_label(label, width=24)
        if len(lines) > 4:
            lines = textwrap.wrap(label, width=22, break_long_words=False, break_on_hyphens=False)[:4]
        text_y = y + node_d + label_gap
        for idx, line in enumerate(lines):
            ax.text(
                x + node_d / 2,
                text_y + idx * label_line_height,
                line,
                fontsize=8.4,
                color="#111111",
                ha="center",
                va="top",
                clip_on=True,
                zorder=3,
            )

    legend_x = left_margin + max_depth * column_dx + 320 + 42
    legend_y = 130
    # ax.text(legend_x, legend_y - 24, "Legend", fontsize=14, fontweight="bold", ha="left", va="top", color="#111111")

    legend_items = [
        (("azimuth",), "Azimuth only"),
        (("ribca",), "RIBCA only"),
        (("frmatch",), "FR-Match only"),
        (("azimuth", "ribca"), "Azimuth + RIBCA"),
        (("azimuth", "frmatch"), "Azimuth + FR-Match"),
        (("frmatch", "ribca"), "RIBCA + FR-Match"),
        (("azimuth", "frmatch", "ribca"), "Azimuth + RIBCA + FR-Match"),
    ]
    combo_counts = Counter(normalize_combo(v) for v in node_sources.values())

    item_y = legend_y
    for combo, label in legend_items:
        color = COMBO_COLORS.get(combo, "#bdbdbd")
        count = combo_counts.get(combo, 0)
        ax.add_patch(
            Circle(
                (legend_x + 8, item_y - 5),
                radius=7,
                linewidth=1.0,
                edgecolor="#202124",
                facecolor=color,
            )
        )
        ax.text(legend_x + 22, item_y - 1, f"{label} ({count} nodes)", fontsize=9.5, ha="left", va="top", color="#222222")
        item_y += 26

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(SVG_OUT, format="svg")
    fig.savefig(PNG_OUT, format="png", dpi=160)
    plt.close(fig)
    print(f"Saved {SVG_OUT}")
    print(f"Saved {PNG_OUT}")


if __name__ == "__main__":
    main()
