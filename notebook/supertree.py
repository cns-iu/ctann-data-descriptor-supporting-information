from __future__ import annotations

import csv
from collections import Counter, defaultdict, deque
from datetime import datetime
from pathlib import Path

import matplotlib


matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "DejaVu Sans"

INPUT = "ctann-v8-3.csv"

OUTPUT_DIR = "vis"
SVG_OUT = OUTPUT_DIR / "ctann-v8-all-sources-martin-blue.svg"
PNG_OUT = OUTPUT_DIR / "ctann-v8-all-sources-martin-blue.png"

DEFAULT_NODE_FILL_COLOR = "#7A7B78"
MARTIN_NODE_FILL_COLOR = "#1E90FF"


def main() -> None:
    rows = list(csv.DictReader(INPUT.open(newline="", encoding="utf-8-sig")))

    node_sources: dict[str, set[str]] = defaultdict(set)

    node_label: dict[str, str] = {}

    # graph structure
    parents: dict[str, set[str]] = defaultdict(set)
    children: dict[str, set[str]] = defaultdict(set)
    edge_counts: Counter[tuple[str, str]] = Counter()

    source_counts: Counter[str] = Counter()

    row_count = 0

    for row in rows:
        # Include ALL CT/1 - Sources, no filtering
        source_raw = row.get("CT/1 - Sources", "").strip()
        source_display = source_raw if source_raw else "(blank)"
        source_norm = source_display.lower()

        row_count += 1
        source_counts[source_display] += 1

        path: list[tuple[str, str]] = []

        for idx in range(1, 13):
            node_id = row.get(f"AS/{idx}/ID", "").strip()
            node_lab = row.get(f"AS/{idx}/LABEL", "").strip()
            if node_id:
                path.append((node_id, node_lab))

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

    left_margin = 760
    top_margin = 140
    right_legend_space = 620
    column_dx = 292.5
    row_dy = 8.5
    node_d = 14
    label_gap = 2

    width = left_margin + max_depth * column_dx + 360 + right_legend_space
    layout_height = max(1.0, max_y - min_y) * row_dy
    height = top_margin + layout_height + 220

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
        if "martin" in node_sources[node]:
            fill_color = MARTIN_NODE_FILL_COLOR
        else:
            fill_color = DEFAULT_NODE_FILL_COLOR

        ax.add_patch(
            Circle(
                (x + node_d / 2, y + node_d / 2),
                radius=node_d / 2,
                linewidth=0,
                edgecolor="none",
                facecolor=fill_color,
                zorder=2,
            )
        )

        label = node_label.get(node, node)

        ax.text(
            x + node_d / 2,
            y + node_d + label_gap,
            label,
            fontsize=6.8,
            color="#111111",
            ha="center",
            va="top",
            clip_on=False,
            zorder=3,
        )

    # Legend
    legend_x = 40
    legend_y = 50

    martin_node_count = sum(1 for sources in node_sources.values() if "martin" in sources)
    non_martin_node_count = len(node_sources) - martin_node_count

    legend_items = [
        (
            DEFAULT_NODE_FILL_COLOR,
            f"All non-Martin nodes ({non_martin_node_count} nodes)",
        ),
        (
            MARTIN_NODE_FILL_COLOR,
            f'Martin nodes from CT/1 - Sources = "Martin" ({martin_node_count} nodes)',
        ),
    ]

    item_y = legend_y

    for color, label in legend_items:
        ax.add_patch(
            Circle(
                (legend_x + 8, item_y - 5),
                radius=7,
                linewidth=0,
                edgecolor="none",
                facecolor=color,
            )
        )

        ax.text(
            legend_x + 22,
            item_y - 1,
            label,
            fontsize=19,
            ha="left",
            va="top",
            color="#222222",
        )

        item_y += 28

    item_y += 16

    ax.text(
        legend_x,
        item_y,
        "Included CT/1 - Sources:",
        fontsize=19,
        fontweight="bold",
        ha="left",
        va="top",
        color="#222222",
    )

    item_y += 24

    for source, count in sorted(source_counts.items(), key=lambda x: x[0].lower()):
        ax.text(
            legend_x,
            item_y,
            f"{source}: {count} rows",
            fontsize=18,
            ha="left",
            va="top",
            color="#444444",
        )
        item_y += 21

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig.savefig(SVG_OUT, format="svg")
    fig.savefig(PNG_OUT, format="png", dpi=180)

    plt.close(fig)

    print(f"Saved {SVG_OUT}")
    print(f"Saved {PNG_OUT}")


if __name__ == "__main__":
    main()
