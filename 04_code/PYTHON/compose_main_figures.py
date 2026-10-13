"""Compose the four main figures with publication-safe Pillow layouts.

The statistical source data are preserved in ``03_source_data``. This module
only controls drawing, labels, colours and panel composition.
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from PIL import JpegImagePlugin  # register JPEG encoder used by Pillow PDF output


def font(size: int, bold: bool = False):
    custom = os.environ.get("FIGURE_FONT_BOLD" if bold else "FIGURE_FONT_REGULAR")
    windows_font_dir = Path(os.environ.get("WINDIR", "")) / "Fonts"
    paths = [
        windows_font_dir / ("arialbd.ttf" if bold else "arial.ttf"),
        windows_font_dir / ("calibrib.ttf" if bold else "calibri.ttf"),
        Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
    ]
    if custom:
        paths.insert(0, Path(custom))
    for path in paths:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    raise FileNotFoundError("A scalable font is required; set FIGURE_FONT_REGULAR and FIGURE_FONT_BOLD.")


def fit(im: Image.Image, max_w: int, max_h: int) -> Image.Image:
    scale = min(max_w / im.width, max_h / im.height)
    return im.resize((max(1, int(im.width * scale)), max(1, int(im.height * scale))), Image.Resampling.LANCZOS)


def save(canvas: Image.Image, out: Path, stem: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    canvas.save(out / f"{stem}.png", dpi=(300, 300))
    canvas.save(out / f"{stem}.tiff", compression="tiff_lzw", dpi=(300, 300))
    canvas.save(out / f"{stem}.pdf", resolution=300.0)


def draw_label(draw, xy, value: str, fnt, fill="black") -> None:
    """Render O3 with a real raised baseline for O and lowered small 3."""
    x, y = xy
    parts = value.replace("O\u2083", "O3").split("O3")
    for i, part in enumerate(parts):
        draw.text((x, y), part, font=fnt, fill=fill)
        x += draw.textlength(part, font=fnt)
        if i < len(parts) - 1:
            draw.text((x, y), "O", font=fnt, fill=fill)
            x += draw.textlength("O", font=fnt)
            small = font(round(fnt.size * 0.7))
            draw.text((x, y + round(fnt.size * 0.36)), "3", font=small, fill=fill)
            x += draw.textlength("3", font=small)


def label_width(draw, value: str, fnt) -> float:
    parts = value.replace("O\u2083", "O3").split("O3")
    return sum(draw.textlength(p, font=fnt) for p in parts) + (len(parts) - 1) * (
        draw.textlength("O", font=fnt) + draw.textlength("3", font=font(round(fnt.size * 0.7))))


def wrap_to_width(draw, value: str, fnt, max_width: int) -> list[str]:
    result = []
    for paragraph in value.split("\n"):
        line = ""
        for word in paragraph.split():
            candidate = f"{line} {word}".strip()
            if line and label_width(draw, candidate, fnt) > max_width:
                result.append(line)
                line = word
            else:
                line = candidate
        result.append(line)
    return result


def text_center(draw, box, value, fnt, fill="black", padding=14):
    x0, y0, x1, y1 = box
    lines = wrap_to_width(draw, value, fnt, x1 - x0 - 2 * padding)
    line_h = round(fnt.size * 1.2)
    top = (y0 + y1 - len(lines) * line_h) / 2
    for i, line in enumerate(lines):
        draw_label(draw, ((x0 + x1 - label_width(draw, line, fnt)) / 2,
                          top + i * line_h), line, fnt, fill)


def wrap_words(value: str, width: int = 30) -> list[str]:
    words = value.split()
    lines, current = [], ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if current and len(candidate) > width:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def draw_arrow(draw, start, end, fill="#555555", width=5, dashed=False):
    x1, y1 = start
    x2, y2 = end
    if dashed:
        distance = max(abs(x2 - x1), abs(y2 - y1))
        n = max(1, int(distance / 28))
        for i in range(0, n, 2):
            a, b = i / n, min(1, (i + 1) / n)
            draw.line((x1 + (x2-x1)*a, y1 + (y2-y1)*a,
                       x1 + (x2-x1)*b, y1 + (y2-y1)*b), fill=fill, width=width)
    else:
        draw.line((x1, y1, x2, y2), fill=fill, width=width)
    if abs(x2-x1) >= abs(y2-y1):
        s = 20 if x2 >= x1 else -20
        pts = [(x2, y2), (x2-s, y2-11), (x2-s, y2+11)]
    else:
        s = 20 if y2 >= y1 else -20
        pts = [(x2, y2), (x2-11, y2-s), (x2+11, y2-s)]
    draw.polygon(pts, fill=fill)


def draw_poly_arrow(draw, points, fill="#555555", width=5, dashed=False):
    for a, b in zip(points, points[1:]):
        draw_arrow(draw, a, b, fill=fill, width=width, dashed=dashed)


def draw_route(draw, points, fill="#555555", width=5, dashed=False):
    for start, end in zip(points, points[1:]):
        if dashed:
            dx, dy = end[0] - start[0], end[1] - start[1]
            count = max(1, int((dx * dx + dy * dy) ** 0.5 / 22))
            for i in range(0, count, 2):
                a, b = i / count, min(1, (i + 1) / count)
                draw.line((start[0] + dx * a, start[1] + dy * a,
                           start[0] + dx * b, start[1] + dy * b), fill=fill, width=width)
        else:
            draw.line((*start, *end), fill=fill, width=width)
    x1, y1 = points[-2]
    x2, y2 = points[-1]
    dx, dy = x2 - x1, y2 - y1
    length = max(1, (dx * dx + dy * dy) ** 0.5)
    ux, uy = dx / length, dy / length
    head = [(x2, y2),
            (x2 - 20 * ux - 10 * uy, y2 - 20 * uy + 10 * ux),
            (x2 - 20 * ux + 10 * uy, y2 - 20 * uy - 10 * ux)]
    draw.polygon(head, fill=fill)


def bezier(start, control1, control2, end, count=45):
    points = []
    for i in range(count + 1):
        t = i / count
        s = 1 - t
        points.append((s**3 * start[0] + 3*s*s*t*control1[0] +
                       3*s*t*t*control2[0] + t**3*end[0],
                       s**3 * start[1] + 3*s*s*t*control1[1] +
                       3*s*t*t*control2[1] + t**3*end[1]))
    return points


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def draw_flow_panel(draw, root: Path, x: int, y: int, w: int, h: int) -> None:
    rows = read_rows(root / "03_source_data" / "FigureS1_participant_flow.csv")
    draw.text((x, y), "a", font=font(54, True), fill="black")
    draw.text((x + 75, y), "Participant flow and sample selection", font=font(50, True), fill="black")
    left, right = x + 30, x + w - 40
    main_right = x + int(w * 0.56)
    excluded_left = x + int(w * 0.64)
    box_h = 205
    top = y + 110
    gap = max(30, int((h - 150 - len(rows) * box_h) / max(1, len(rows) - 1)))
    screening_fill, final_fill, excluded_fill = "#d9e8f5", "#a9cbe5", "#b1c5d6"
    border = "#39739d"
    for idx, row in enumerate(rows):
        by = top + idx * (box_h + gap)
        step = int(row["step"])
        is_final = step == 8
        fill = final_fill if is_final else screening_fill
        draw.rounded_rectangle((left, by, main_right, by + box_h), radius=16, fill=fill, outline=border, width=5)
        size = 46
        name_lines = wrap_to_width(draw, row["criterion"], font(size), main_right - left - 38)
        while (len(name_lines) + 1) * round(size * 1.15) > box_h - 18 and size > 34:
            size -= 2
            name_lines = wrap_to_width(draw, row["criterion"], font(size), main_right - left - 38)
        if step == 7:
            final_line = f"n = {int(row['n_obs']):,} obs; 3,284 participants"
        elif step == 8:
            final_line = f"n = {int(row['n_obs']):,} obs; 3,281 participants"
        else:
            final_line = f"n = {int(row['n_obs']):,}"
        label_lines = name_lines + [final_line]
        line_h = round(size * 1.15)
        total_height = len(label_lines) * line_h
        yy = by + (box_h - total_height) / 2
        for j, line in enumerate(label_lines):
            fnt = font(size, bold=(j == len(label_lines) - 1))
            draw_label(draw, ((left + main_right - label_width(draw, line, fnt)) / 2,
                              yy + j * line_h), line, fnt)
        if idx < len(rows) - 1:
            draw_arrow(draw, ((left + main_right)//2, by + box_h + 8),
                       ((left + main_right)//2, by + box_h + gap - 8), width=4)
        if step >= 2:
            excluded = f"Excluded: {int(row['excluded_from_previous']):,}\n{row['exclusion_reason']}"
            ey = by + 5
            draw.rounded_rectangle((excluded_left, ey, right, ey + box_h - 10), radius=12,
                                   fill=excluded_fill, outline="#6e879b", width=3)
            excluded_size = 38
            while len(wrap_to_width(draw, excluded, font(excluded_size), right - excluded_left - 52)) * round(excluded_size * 1.2) > box_h - 34 and excluded_size > 30:
                excluded_size -= 2
            text_center(draw, (excluded_left + 18, ey + 8, right - 18, ey + box_h - 18),
                        excluded, font(excluded_size), fill="#233746")
            draw_arrow(draw, (main_right + 15, by + box_h // 2), (excluded_left - 12, by + box_h // 2),
                       fill="#6e879b", width=3)


def draw_timing_panel(draw, x: int, y: int, w: int, h: int) -> None:
    draw.text((x, y), "b", font=font(54, True), fill="black")
    draw.text((x + 75, y), "Study timing", font=font(50, True), fill="black")
    x0, x1, base = x + 130, x + w - 100, y + 300
    draw.line((x0, base, x1, base), fill="#333333", width=6)
    years = [2007, 2010, 2011, 2015]
    xs = [x0 + int((x1 - x0) * (year - 2007) / 8) for year in years]
    for year, xx in zip(years, xs):
        draw.line((xx, base - 28, xx, base + 28), fill="#333333", width=5)
        bb = draw.textbbox((0, 0), str(year), font=font(44))
        draw.text((xx - (bb[2]-bb[0]) / 2, base + 45), str(year), font=font(44), fill="black")
    draw.line((xs[0], base - 105, xs[1], base - 105), fill="#2b8cbe", width=25)
    draw_label(draw, (xs[0], base - 200), "4-year mean residential EVI / O3", font(48, True), "#1d5f85")
    draw_arrow(draw, (xs[1], base - 105), (xs[2] - 30, base - 28), fill="#555555", width=4)
    for xx in (xs[2], xs[3]):
        draw.ellipse((xx - 20, base - 20, xx + 20, base + 20), fill="#333333")
    draw.text((xs[2] - 105, base + 105), "BAA wave 1", fill="black", font=font(44))
    draw.text((xs[3] - 105, base + 105), "BAA wave 2", fill="black", font=font(44))


def draw_node(draw, box, label, fill="#e8f1f8", outline="#39739d", dashed=False, size=44):
    x0, y0, x1, y1 = box
    draw.rounded_rectangle((x0, y0, x1, y1), radius=18, fill=fill, outline=outline, width=5)
    if dashed:
        draw.rounded_rectangle((x0 + 4, y0 + 4, x1 - 4, y1 - 4), radius=14, outline="#687d8b", width=3)
    text_center(draw, box, label, font(size), fill="black")


def draw_dag_panel(draw, x: int, y: int, w: int, h: int) -> None:
    draw.text((x, y), "c", font=font(54, True), fill="black")
    draw.text((x + 75, y), "Directed acyclic graph (analysis set)", font=font(50, True), fill="black")
    conf = (x + 80, y + 510, x + 510, y + 950)
    evi = (x + 650, y + 320, x + 1030, y + 500)
    o3 = (x + 650, y + 1245, x + 1030, y + 1445)
    baa = (x + 1290, y + 620, x + 1590, y + 820)
    med = (x + 1140, y + 1210, x + 1580, y + 1480)
    center_x = lambda box: (box[0] + box[2]) // 2
    center_y = lambda box: (box[1] + box[3]) // 2
    draw_route(draw, [(center_x(conf), conf[1]), (evi[0], center_y(evi))], width=5)
    draw_route(draw, [(center_x(conf), conf[3]), (o3[0], center_y(o3))], width=5)
    draw_route(draw, [(conf[2], center_y(conf)), (baa[0], center_y(baa))], width=5)
    draw_route(draw, [(evi[2], center_y(evi)), (center_x(baa), baa[1])], width=5)
    draw_route(draw, [(o3[2], center_y(o3)), (center_x(baa), baa[3])], width=5)
    draw_route(draw, [(center_x(evi), evi[3]), (center_x(med), med[1])],
               fill="#687d8b", width=4, dashed=True)
    draw_route(draw, [(center_x(med), med[1]), (center_x(baa), baa[3])],
               fill="#687d8b", width=4, dashed=True)
    draw_node(draw, conf, "Measured adjustment set\nAge, sex, education;\nconsumption, retirement,\npension, climate;\nBMI, behaviours,\ncomorbidity, survey year", fill="#eef4f8", size=38)
    draw_node(draw, evi, "Residential EVI")
    draw_node(draw, o3, "Residential O3")
    draw_node(draw, baa, "BAA", fill="#d9e8f5")
    draw_node(draw, med, "Candidate mediators\nMetabolic activity,\nsocial participation", fill="#f7f7f7", outline="#687d8b", dashed=True, size=40)
    note = ("Solid arrows: working adjustment and exposure-outcome assumptions.\n"
            "Dashed arrows: exploratory pathways outside the primary models.")
    text_center(draw, (x + 55, y + h - 140, x + w - 20, y + h - 20),
                note, font(40), fill="#333333")


def figure1(root: Path, out: Path) -> None:
    canvas = Image.new("RGB", (3600, 3000), "white")
    draw = ImageDraw.Draw(canvas)
    draw_flow_panel(draw, root, 70, 55, 1570, 2830)
    draw_timing_panel(draw, 1790, 70, 1650, 980)
    draw_dag_panel(draw, 1790, 1130, 1650, 1760)
    save(canvas, out, "Figure1_study_design")


def draw_means_panel(canvas, draw, rows, box, panel_label, title, show_y_label=False):
    x, y, w, h = box
    draw.text((x, y), panel_label, font=font(52, True), fill="black")
    draw.text((x + 75, y), title, font=font(50, True), fill="black")
    plot_left, plot_right = x + 270, x + w - 110
    plot_top, plot_bottom = y + 180, y + h - 340
    ymin, ymax = -3.0, 2.3
    def py(v):
        return plot_bottom - (v - ymin) / (ymax - ymin) * (plot_bottom - plot_top)
    for tick in [-3, -2, -1, 0, 1, 2]:
        yy = py(tick)
        draw.line((plot_left, yy, plot_right, yy), fill="#e2e2e2" if tick else "#b7b7b7", width=3 if tick else 4)
        bb = draw.textbbox((0, 0), str(tick), font=font(44))
        draw.text((plot_left - 50 - bb[2], yy - 26), str(tick), font=font(44), fill="black")
    draw.line((plot_left, plot_top, plot_left, plot_bottom), fill="#333333", width=4)
    draw.line((plot_left, plot_bottom, plot_right, plot_bottom), fill="#333333", width=4)
    if show_y_label:
        axis_label = "Adjusted marginal mean BAA (years)"
        label = Image.new("RGBA", (1200, 80), (255, 255, 255, 0))
        ImageDraw.Draw(label).text((0, 0), axis_label, font=font(45), fill="#333333")
        rotated = label.rotate(90, expand=True)
        canvas.paste(rotated, (x + 18, round((plot_top + plot_bottom - rotated.height) / 2)), rotated)
    colors = ["#1682bd", "#e2a300", "#16a57a", "#c96aa5"]
    shapes = ["circle", "diamond", "triangle", "square"]
    labels = ["Low EVI / low O3", "Low EVI / high O3", "High EVI / low O3", "High EVI / high O3"]
    inner_left, inner_right = plot_left + 85, plot_right - 85
    xs = [inner_left + (inner_right - inner_left) * i / 3 for i in range(4)]
    for i, row in enumerate(rows):
        xx = xs[i]
        mean, lci, uci = float(row["adjusted_mean"]), float(row["lci"]), float(row["uci"])
        ymean, yl, yu = py(mean), py(lci), py(uci)
        draw.line((xx, yu, xx, yl), fill="#666666", width=6)
        draw.line((xx - 18, yu, xx + 18, yu), fill="#666666", width=5)
        draw.line((xx - 18, yl, xx + 18, yl), fill="#666666", width=5)
        r = 18
        if shapes[i] == "circle":
            draw.ellipse((xx-r, ymean-r, xx+r, ymean+r), fill=colors[i])
        elif shapes[i] == "diamond":
            draw.polygon([(xx, ymean-r-3), (xx+r+3, ymean), (xx, ymean+r+3), (xx-r-3, ymean)], fill=colors[i])
        elif shapes[i] == "triangle":
            draw.polygon([(xx, ymean-r-4), (xx+r+3, ymean+r+5), (xx-r-3, ymean+r+5)], fill=colors[i])
        else:
            draw.rectangle((xx-r, ymean-r, xx+r, ymean+r), fill=colors[i])
        first, second = labels[i].split(" / ")
        fnt = font(44)
        draw_label(draw, (xx - label_width(draw, first, fnt)/2, plot_bottom + 28), first, fnt)
        draw_label(draw, (xx - label_width(draw, second, fnt)/2, plot_bottom + 84), second, fnt)
    return colors, shapes, labels


def draw_contrast_forest(draw, rows, box):
    """Render adjusted pairwise contrasts and the DiD estimate from public CSVs."""
    x, y, w, h = box
    draw.text((x, y), "c", font=font(48, True), fill="black")
    draw.text((x + 68, y), "Adjusted contrasts and interaction", font=font(44, True), fill="black")
    label_x = x + 30
    plot_left, plot_right = x + 630, x + 1460
    text_x = x + 1510
    xmin, xmax = -3.5, 3.5
    axis = lambda value: plot_left + (float(value) - xmin) / (xmax - xmin) * (plot_right - plot_left)
    keys = [
        ("Low EVI/high O3 vs reference", "Low EVI / high O3 - ref."),
        ("High EVI/low O3 vs reference", "High EVI / low O3 - ref."),
        ("High EVI/high O3 vs reference", "High EVI / high O3 - ref."),
        ("O3 contrast within low EVI", "O3 within low EVI"),
        ("O3 contrast within high EVI", "O3 within high EVI"),
        ("Difference-in-differences interaction", "DiD interaction"),
    ]
    rows_by = {(r["scale"], r["contrast"]): r for r in rows}
    panel_starts = {"1KM": y + 135, "500m": y + 625}
    colours = {"pair": "#365f7d", "within": "#16a57a", "did": "#d27b2c"}
    for scale, scale_label in (("1KM", "1 km"), ("500m", "500 m")):
        top = panel_starts[scale]
        draw.text((label_x, top - 60), scale_label, font=font(35, True), fill="#233746")
        for j, (key, label) in enumerate(keys):
            row = rows_by[(scale, key)]
            yy = top + j * 76
            kind = "did" if key.startswith("Difference") else ("within" if "within" in key else "pair")
            colour = colours[kind]
            draw_label(draw, (label_x, yy - 24), label, font(40), fill="#333333")
            draw.line((axis(row["lci"]), yy, axis(row["uci"]), yy), fill=colour, width=5)
            draw.line((axis(row["lci"]), yy - 12, axis(row["lci"]), yy + 12), fill=colour, width=4)
            draw.line((axis(row["uci"]), yy - 12, axis(row["uci"]), yy + 12), fill=colour, width=4)
            xx = axis(row["estimate"])
            if kind == "did":
                draw.polygon([(xx, yy - 15), (xx - 14, yy + 12), (xx + 14, yy + 12)], fill=colour)
            else:
                draw.ellipse((xx - 11, yy - 11, xx + 11, yy + 11), fill=colour)
            p = float(row["p"])
            p_text = "P < 0.001" if p < 0.001 else f"P = {p:.3f}"
            draw.text((text_x, yy - 24), f"{float(row['estimate']):+.2f} [{float(row['lci']):+.2f}, {float(row['uci']):+.2f}]  {p_text}",
                      font=font(38), fill="#333333")

    axis_y = y + h - 130
    draw.line((plot_left, axis_y, plot_right, axis_y), fill="#333333", width=3)
    draw.line((axis(0), y + 105, axis(0), axis_y), fill="#b8b8b8", width=3)
    for tick in (-3, -2, -1, 0, 1, 2, 3):
        xx = axis(tick)
        draw.line((xx, axis_y, xx, axis_y + 12), fill="#333333", width=3)
        draw.text((xx - 12, axis_y + 18), str(tick), font=font(35), fill="#333333")
    draw.text((plot_left + 180, axis_y + 58), "Adjusted contrast in BAA (years)", font=font(31), fill="#333333")
    draw.text((label_x, y + h - 5), "Reference: low EVI / low O3. Circles: contrasts; triangles: DiD interaction.",
              font=font(32), fill="#555555")


def draw_precision_panel(draw, rows, box):
    """Show DiD precision, detectable interaction size and the prespecified margin."""
    x, y, w, h = box
    draw.text((x, y), "d", font=font(48, True), fill="black")
    draw.text((x + 68, y), "DiD precision", font=font(44, True), fill="black")
    xmin, xmax = -3.8, 3.8
    left, right = x + 150, x + w - 55
    axis = lambda value: left + (float(value) - xmin) / (xmax - xmin) * (right - left)
    top, bottom = y + 230, y + 620
    # The ±1-year band is a prespecified equivalence margin, not a result.
    draw.rectangle((axis(-1), top - 55, axis(1), bottom + 15), fill="#f3f0df")
    draw.text((axis(-1) + 8, top - 95), "±1 y margin", font=font(35), fill="#796f3e")
    draw.line((axis(0), top - 55, axis(0), bottom + 15), fill="#777777", width=3)
    for bound in (-1, 1):
        draw.line((axis(bound), top - 40, axis(bound), bottom), fill="#b8ae74", width=2)
    for i, scale in enumerate(("1KM", "500m")):
        row = next(r for r in rows if r["scale"] == scale)
        yy = top + i * 205
        colour = "#d27b2c"
        draw.text((x + 25, yy - 20), "1 km" if scale == "1KM" else "500 m", font=font(31, True), fill="#233746")
        draw.line((axis(row["lci"]), yy, axis(row["uci"]), yy), fill=colour, width=6)
        draw.line((axis(row["lci"]), yy - 14, axis(row["lci"]), yy + 14), fill=colour, width=4)
        draw.line((axis(row["uci"]), yy - 14, axis(row["uci"]), yy + 14), fill=colour, width=4)
        xx = axis(row["estimate"])
        draw.polygon([(xx, yy - 16), (xx - 15, yy + 13), (xx + 15, yy + 13)], fill=colour)
        draw.text((x + 25, yy + 35), f"{float(row['estimate']):+.2f} [{float(row['lci']):+.2f}, {float(row['uci']):+.2f}]",
                  font=font(35), fill="#333333")
    axis_y = bottom + 20
    draw.line((left, axis_y, right, axis_y), fill="#333333", width=3)
    for tick in (-3, -2, -1, 0, 1, 2, 3):
        xx = axis(tick)
        draw.line((xx, axis_y, xx, axis_y + 10), fill="#333333", width=3)
        draw.text((xx - 10, axis_y + 15), str(tick), font=font(35), fill="#333333")
    draw.text((left + 70, axis_y + 52), "DiD estimate (years)", font=font(28), fill="#333333")
    mde = {r["scale"]: float(r["mde_80pct"]) for r in rows}
    tost = {r["scale"]: float(r["tost_p"]) for r in rows}
    note = (f"80% MDE: {mde['1KM']:.2f} y (1 km); {mde['500m']:.2f} y (500 m)\n"
            f"TOST P (±1 y): {tost['1KM']:.3f}; {tost['500m']:.3f}\n"
            "Shaded band is the prespecified margin; equivalence was not established.")
    lines = wrap_to_width(draw, note, font(32), w - 50)
    for j, line in enumerate(lines):
        draw.text((x + 25, y + h - 260 + j * 43), line, font=font(32), fill="#555555")


def figure3(root: Path, out: Path) -> None:
    rows = read_rows(root / "03_source_data" / "Figure3a_adjusted_means.csv")
    contrasts = read_rows(root / "03_source_data" / "Figure3b_contrasts.csv")
    precision = read_rows(root / "03_source_data" / "Figure3b_DiD_precision.csv")
    rows_1km = sorted([r for r in rows if r["scale"] == "1KM"], key=lambda r: int(r["group_order"]))
    rows_500 = sorted([r for r in rows if r["scale"] == "500m"], key=lambda r: int(r["group_order"]))
    for frame in [rows_1km, rows_500]:
        if [int(row["group_order"]) for row in frame] != [1, 2, 3, 4]:
            raise ValueError("Expected four unique ordered exposure groups")
        if any(int(row["n"]) != 6471 for row in frame):
            raise ValueError("Joint-model frame differs from 6,471 observations")

    expected_contrasts = {
        "Low EVI/high O3 vs reference",
        "High EVI/low O3 vs reference",
        "High EVI/high O3 vs reference",
        "O3 contrast within low EVI",
        "O3 contrast within high EVI",
        "Difference-in-differences interaction",
    }
    if {r["contrast"] for r in contrasts} != expected_contrasts:
        raise ValueError("Figure 3 contrast source does not contain the expected six contrasts")
    if {r["scale"] for r in precision} != {"1KM", "500m"}:
        raise ValueError("Figure 3 precision source must contain both spatial scales")

    canvas = Image.new("RGB", (3600, 3200), "white")
    draw = ImageDraw.Draw(canvas)

    # Hero row: the descriptive four-cell exposure pattern estimates.
    colors, shapes, labels = draw_means_panel(
        canvas, draw, rows_1km, (80, 65, 1650, 1500), "a", "1 km residential buffer", show_y_label=True)
    draw_means_panel(
        canvas, draw, rows_500, (1870, 65, 1650, 1500), "b", "500 m residential buffer")

    # Compact colour/shape key between the descriptive and inferential rows.
    legend_y = 1605
    draw.text((80, legend_y), "Joint exposure pattern", font=font(38, True), fill="#333333")
    for i, (color, shape, label) in enumerate(zip(colors, shapes, labels)):
        xx = 690 + (i % 2) * 790
        yy = legend_y + (i // 2) * 60
        if shape == "circle": draw.ellipse((xx, yy, xx+22, yy+22), fill=color)
        elif shape == "diamond": draw.polygon([(xx+11, yy-3), (xx+25, yy+11), (xx+11, yy+25), (xx-3, yy+11)], fill=color)
        elif shape == "triangle": draw.polygon([(xx+11, yy-4), (xx+26, yy+24), (xx-4, yy+24)], fill=color)
        else: draw.rectangle((xx, yy, xx+22, yy+22), fill=color)
        draw_label(draw, (xx + 34, yy - 8), label, font(35))

    draw_contrast_forest(draw, contrasts, (80, 1740, 2260, 1260))
    draw_precision_panel(draw, precision, (2460, 1740, 1080, 1260))
    draw.text((80, 3070),
              "Points: adjusted marginal means; whiskers: 95% CI. Contrast panels use the common joint-model frame (6,471 observations).",
              font=font(34), fill="#333333")
    save(canvas, out, "Figure3_joint_patterns")


def draw_forest_panel(draw, rows, top: int) -> None:
    left, right = 880, 2390
    x_min, x_max = -1.25, 2.15
    axis = lambda value: left + (value - x_min) / (x_max - x_min) * (right - left)
    groups = [("EVI", "0.053 at 1 km; 0.048 at 500 m"), ("O3", "6.8 at 1 km; 7.3 at 500 m, µg/m³")]
    row_names = {("1KM", "unweighted"): "1 km, unweighted",
                 ("1KM", "ipw"): "1 km, IPW",
                 ("500m", "unweighted"): "500 m, unweighted",
                 ("500m", "ipw"): "500 m, IPW"}
    for group_index, (exposure, iqr) in enumerate(groups):
        group_top = top + 90 + group_index * 760
        draw_label(draw, (150, group_top), f"{exposure}  (IQR = {iqr})",
                   font(52, True), "#233746")
        for row_index, key in enumerate(row_names):
            row = next(r for r in rows if r["exposure"] == exposure and
                       (r["scale"], r["weight_spec"]) == key)
            yy = group_top + 145 + row_index * 140
            primary = key[1] == "unweighted"
            color = "#1d5f85" if primary else "#687d8b"
            draw.text((150, yy - 30), row_names[key], font=font(44), fill="#333333")
            estimate = float(row["beta_iqr"])
            lower, upper = float(row["lci_iqr"]), float(row["uci_iqr"])
            draw.line((axis(lower), yy, axis(upper), yy), fill=color, width=6)
            for bound in (lower, upper):
                draw.line((axis(bound), yy - 16, axis(bound), yy + 16), fill=color, width=5)
            xx = axis(estimate)
            if primary:
                draw.ellipse((xx - 17, yy - 17, xx + 17, yy + 17), fill=color)
            else:
                draw.polygon([(xx, yy - 21), (xx - 20, yy + 17),
                              (xx + 20, yy + 17)], fill=color)
            result = f"{estimate:+.2f} [{lower:+.2f}, {upper:+.2f}], P = {float(row['p']):.3f}"
            draw.text((2470, yy - 28), result, font=font(42), fill="#333333")
    axis_y = top + 1630
    draw.line((left, top + 165, left, axis_y), fill="#333333", width=4)
    draw.line((left, axis_y, right, axis_y), fill="#333333", width=4)
    draw.line((axis(0), top + 165, axis(0), axis_y), fill="#969696", width=4)
    for tick in (-1, 0, 1, 2):
        xx = axis(tick)
        draw.line((xx, axis_y, xx, axis_y + 17), fill="#333333", width=4)
        draw.text((xx - 15, axis_y + 26), str(tick), font=font(42), fill="#333333")
    draw.text((1040, axis_y + 90), "Adjusted BAA difference per IQR increase (years)",
              font=font(44), fill="#333333")


def figure2(root: Path, out: Path) -> None:
    from render_public_distributions import raincloud_panel
    raincloud = raincloud_panel(root)
    raincloud = fit(raincloud, 3380, 4000)
    forest_top = 130 + raincloud.height + 105
    canvas = Image.new("RGB", (3600, forest_top + 1850), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((85, 45), "a", font=font(54, True), fill="black")
    draw.text((165, 45), "Observed BAA by joint exposure pattern", font=font(50, True), fill="black")
    canvas.paste(raincloud, ((3600 - raincloud.width) // 2, 130))
    draw.text((85, forest_top - 65), "b", font=font(54, True), fill="black")
    draw.text((165, forest_top - 65), "Continuous-exposure associations", font=font(50, True), fill="black")
    draw_forest_panel(draw, read_rows(root / "03_source_data" / "Figure2_continuous_forest.csv"), forest_top)
    save(canvas, out, "Figure2_primary_associations")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path.cwd())
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    root = args.root.resolve()
    out = args.out.resolve() if args.out else root / "02_figures" / "main"
    figure1(root, out)
    figure2(root, out)
    figure3(root, out)
    print(f"Wrote revised Figure 1–3 composites to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
