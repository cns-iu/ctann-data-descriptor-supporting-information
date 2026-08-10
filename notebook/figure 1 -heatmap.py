from __future__ import annotations

import csv
from collections import Counter, defaultdict, deque
from datetime import datetime
from pathlib import Path

import matplotlib


matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 10.5

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent.parent / "ctann" / "ctann-v9.csv"

RUN_ID = datetime.now().strftime("%Y%m%d-%H%M%S")
OUTPUT_DIR = ROOT.parent / "vis"
SVG_OUT = OUTPUT_DIR / "ctann-v9-top3-azimuth-pan-human-crosswalk-fonts-node1_5x-labelright-leveldistance250.svg"
PNG_OUT = OUTPUT_DIR / "ctann-v9-top3-azimuth-pan-human-crosswalk-fonts-node1_5x-labelright-leveldistance250.png"

INCLUDED_SOURCES = {"azimuth", "pan-human-azimuth"}
SOURCE_DISPLAY_NAMES = {
    "azimuth": "Azimuth",
    "pan-human-azimuth": "Pan-Human Azimuth",
}
AZIMUTH_NODE_FILL_COLOR = "#E31A1C"
PAN_HUMAN_AZIMUTH_NODE_FILL_COLOR = "#1E90FF"
OVERLAP_NODE_FILL_COLOR = "#7B2CBF"
LABEL_FONT_SIZE = 10.5


def main() -> None:
    rows = list(csv.DictReader(INPUT.open(newline="", encoding="utf-8-sig")))

    # node -> set of normalized sources (for coloring logic)
    node_sources: dict[str, set[str]] = defaultdict(set)

    # node -> sources where it is the final, most-specific node in a path
    leaf_sources: dict[str, set[str]] = defaultdict(set)

    # node -> label
    node_label: dict[str, str] = {}

    # graph structure
    parents: dict[str, set[str]] = defaultdict(set)
    children: dict[str, set[str]] = defaultdict(set)
    edge_counts: Counter[tuple[str, str]] = Counter()

    # for legend / summary
    source_counts: Counter[str] = Counter()

    row_count = 0

    for row in rows:
        source_raw = row.get("CT/1 - Sources", "").strip()
        source_norm = source_raw.lower()

        if source_norm not in INCLUDED_SOURCES:
            continue

        source_display = SOURCE_DISPLAY_NAMES[source_norm]

        row_count += 1
        source_counts[source_display] += 1

        path: list[tuple[str, str]] = []

        # The figure deliberately ends at the third AS level.
        for idx in range(1, 4):
            node_id = row.get(f"AS/{idx}/ID", "").strip()
            node_lab = row.get(f"AS/{idx}/LABEL", "").strip()
            if node_id:
                path.append((node_id, node_lab))

        if path:
            leaf_sources[path[-1][0]].add(source_norm)

        prev_id: str | None = None
        for node_id, node_lab in path:
            node_sources[node_id].add(source_norm)
            node_label[node_id] = node_lab

            if prev_id and prev_id != node_id:
                parents[node_id].add(prev_id)
                children[prev_id].add(node_id)
                edge_counts[(prev_id, node_id)] += 1

            prev_id = node_id

    indeg = {node: len(parents[node]) for node in node_sources}
    roots = [node for node, deg in indeg.items() if deg == 0]

    if not roots and node_sources:
        roots = [min(node_sources)]

    # Use shortest root distance so a node stays in the leftmost valid layer.
    depth: dict[str, int] = {}
    queue = deque()

    for root in sorted(roots, key=lambda n: (node_label.get(n, ""), n)):
        depth[root] = 0
        queue.append(root)

    while queue:
        node = queue.popleft()

        for child in sorted(
            children[node],
            key=lambda c: (
                -edge_counts[(node, c)],
                node_label.get(c, ""),
                c,
            ),
        ):
            candidate = depth[node] + 1
            if child not in depth or candidate < depth[child]:
                depth[child] = candidate
                queue.append(child)

    # Build one primary parent tree for layout so each subtree can fan out cleanly.
    primary_parent: dict[str, str] = {}

    for node in node_sources:
        if not parents[node]:
            continue

        ranked = sorted(
            parents[node],
            key=lambda p: (
                -edge_counts[(p, node)],
                depth.get(p, 999),
                node_label.get(p, ""),
                p,
            ),
        )
        primary_parent[node] = ranked[0]

    primary_children: dict[str, list[str]] = defaultdict(list)

    for child, parent in primary_parent.items():
        primary_children[parent].append(child)

    for parent in primary_children:
        primary_children[parent].sort(
            key=lambda c: (
                depth.get(c, 999),
                -len(children[c]),
                node_label.get(c, ""),
                c,
            )
        )

    # Count leaves so each subtree gets enough vertical room.
    leaf_weight: dict[str, int] = {}

    def compute_leaf_weight(node: str) -> int:
        kids = primary_children.get(node, [])
        if not kids:
            leaf_weight[node] = 1
            return 1

        total = sum(compute_leaf_weight(child) for child in kids)
        leaf_weight[node] = max(1, total)
        return leaf_weight[node]

    for root in sorted(roots, key=lambda n: (node_label.get(n, ""), n)):
        compute_leaf_weight(root)

    y_pos: dict[str, float] = {}
    cursor = 0.0
    leaf_step = 3.0
    root_gap = 2.0

    def assign_y(node: str) -> None:
        nonlocal cursor

        kids = primary_children.get(node, [])
        if not kids:
            y_pos[node] = cursor
            cursor += leaf_step
            return

        for child in kids:
            assign_y(child)

        y_pos[node] = sum(y_pos[child] for child in kids) / len(kids)

    ordered_roots = sorted(roots, key=lambda n: (node_label.get(n, ""), n))

    for idx, root in enumerate(ordered_roots):
        assign_y(root)
        if idx < len(ordered_roots) - 1:
            cursor += root_gap

    # Any detached nodes fall back below the main layout.
    for node in sorted(
        node_sources,
        key=lambda n: (
            depth.get(n, 999),
            node_label.get(n, ""),
            n,
        ),
    ):
        if node not in y_pos:
            y_pos[node] = cursor
            cursor += leaf_step

    max_depth = max(depth.values()) if depth else 0
    min_y = min(y_pos.values()) if y_pos else 0.0
    max_y = max(y_pos.values()) if y_pos else 1.0

    left_margin = 110
    top_margin = 80
    right_legend_space = 460
    right_label_space = 310
    column_dx = 250
    row_dy = 8.5
    node_d = 21
    label_gap = 4

    width = left_margin + max_depth * column_dx + right_label_space + right_legend_space
    layout_height = max(1.0, max_y - min_y) * row_dy
    height = top_margin + layout_height + 160

    vertical_offset = top_margin + max(
        0.0,
        (height - top_margin - 140 - layout_height) / 2,
    )

    positions: dict[str, tuple[float, float]] = {}

    for node in node_sources:
        x = left_margin + depth.get(node, 0) * column_dx
        y = vertical_offset + (y_pos[node] - min_y) * row_dy
        positions[node] = (x, y)

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    fig = plt.figure(figsize=(width / 120, height / 120), dpi=120)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")

    # Draw edges
    for parent, child in edge_counts:
        if parent not in positions or child not in positions:
            continue

        px, py = positions[parent]
        cx, cy = positions[child]

        x1 = px + node_d / 2
        y1 = py + node_d / 2
        x2 = cx + node_d / 2
        y2 = cy + node_d / 2

        is_primary = primary_parent.get(child) == parent

        ax.plot(
            [x1, x2],
            [y1, y2],
            color="#4a4a4a",
            linewidth=1.0 if is_primary else 0.8,
            alpha=0.56 if is_primary else 0.22,
            zorder=1,
        )

    # Draw nodes
    for node, (x, y) in positions.items():
        center = (x + node_d / 2, y + node_d / 2)
        sources = node_sources[node]
        if sources == INCLUDED_SOURCES:
            fill_color = OVERLAP_NODE_FILL_COLOR
        elif "azimuth" in sources:
            fill_color = AZIMUTH_NODE_FILL_COLOR
        else:
            fill_color = PAN_HUMAN_AZIMUTH_NODE_FILL_COLOR

        ax.add_patch(
            Circle(
                center,
                radius=node_d / 2,
                linewidth=0,
                edgecolor="none",
                facecolor=fill_color,
                zorder=2,
            )
        )

        label = node_label.get(node, node)

        ax.text(
            x + node_d + label_gap,
            y + node_d / 2,
            label,
            fontsize=LABEL_FONT_SIZE,
            color="#111111",
            ha="left",
            va="center",
            clip_on=False,
            zorder=3,
        )

    # Legend
    legend_x = left_margin + max_depth * column_dx + right_label_space + 35
    legend_y = 70

    legend_items = [
        (
            AZIMUTH_NODE_FILL_COLOR,
            f"Azimuth",
        ),
        (
            PAN_HUMAN_AZIMUTH_NODE_FILL_COLOR,
            f"Pan-Human Azimuth",
        ),
        (
            OVERLAP_NODE_FILL_COLOR,
            f"Both",
        ),
    ]

    item_y = legend_y

    for color, label in legend_items:
        ax.add_patch(
            Circle(
                (legend_x + node_d / 2, item_y),
                radius=node_d / 2,
                linewidth=0,
                edgecolor="none",
                facecolor=color,
            )
        )

        ax.text(
            legend_x + node_d + 8,
            item_y,
            label,
            fontsize=LABEL_FONT_SIZE,
            ha="left",
            va="center",
            color="#222222",
        )

        item_y += 36

    # for source, count in sorted(source_counts.items(), key=lambda x: x[0].lower()):
    #     ax.text(
    #         legend_x,
    #         item_y,
    #         f"{source}: {count} rows",
    #         fontsize=18,
    #         ha="left",
    #         va="top",
    #         color="#444444",
    #     )
    #     item_y += 21

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig.savefig(SVG_OUT, format="svg")
    fig.savefig(PNG_OUT, format="png", dpi=180)

    plt.close(fig)

    print(f"Saved {SVG_OUT}")
    print(f"Saved {PNG_OUT}")


if __name__ == "__main__":
    main()
