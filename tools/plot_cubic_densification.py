"""Rebuild the saved Hermite cubic densification chart from solved JSON evidence.

The 5x5 series changes axial resolution only. 7x7x17 is plotted as an
unconnected transverse-refinement observation, NOT a sixth axial h-level.
Neither series measures a pointwise stress-field error or a certified bound.

From the repository root::

    python -m tools.plot_cubic_densification
    python -m tools.plot_cubic_densification --out outputs/cubic_densification.svg
"""
from __future__ import annotations

import argparse
from pathlib import Path

from tools.section_refinement import read_case

ROOT = Path(__file__).resolve().parents[1]
AXIAL_LEVELS = (9, 17, 33, 65, 129)
TRANSVERSE_LEVEL = (7, 7, 17)


def _saved_result(mesh: tuple[int, int, int]) -> Path:
    name = f"cubic_{mesh[0]}x{mesh[1]}x{mesh[2]}"
    standard = ROOT / "samples" / f"{name}.json"
    if standard.is_file():
        return standard
    historical = (ROOT / "samples" / "source_sections_20261026" /
                  f"{name}_section_only.json")
    if historical.is_file():
        return historical
    raise FileNotFoundError(f"No recorded result/section pair for {mesh}; no value will be invented")


def load_measured_cases():
    """Load validated source-section results, keeping plotting data out of code."""
    axial = [read_case(_saved_result((5, 5, n))) for n in AXIAL_LEVELS]
    transverse = read_case(_saved_result(TRANSVERSE_LEVEL))
    for item in [*axial, transverse]:
        if item['positions'] != [0.0]:
            raise ValueError(f"Expected axial-midpoint sections: {item['result_file']}")
    return axial, transverse


def render(out: Path) -> Path:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    axial, other = load_measured_cases()
    xs = list(range(len(axial)))
    n = [int(case['mesh'][2]) for case in axial]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.8, 5.7), layout='constrained')
    shear = [c['maximum_midspan_shear_vs_beam_pct'] for c in axial]
    displacement = [c['tip_vs_beam_pct'] for c in axial]
    moment = [c['maximum_midspan_moment_vs_beam_pct'] for c in axial]
    ax.plot(xs, shear, marker='o', lw=2.3, color='#087087', label='5×5 axial series')
    bx.plot(xs, displacement, marker='o', lw=2.1, color='#b94e20', label='5×5 tip displacement')
    bx.plot(xs, moment, marker='s', lw=2.1, color='#6b52b5', label='5×5 section moment')
    # Independently measured transverse refinement: share the nz=17 tick,
    # but NEVER connect it to the axial series.
    index = n.index(TRANSVERSE_LEVEL[2])
    kw = dict(marker='D', s=91, facecolors='white', linewidths=2.4, zorder=5)
    ax.scatter([index], [other['maximum_midspan_shear_vs_beam_pct']],
               edgecolors='#087087', **kw)
    bx.scatter([index], [other['tip_vs_beam_pct']], edgecolors='#b94e20', **kw)
    bx.scatter([index], [other['maximum_midspan_moment_vs_beam_pct']],
               edgecolors='#6b52b5', **kw)
    ax.annotate('7×7×17\n2.080%', (index, other['maximum_midspan_shear_vs_beam_pct']),
                xytext=(17, 24), textcoords='offset points', fontsize=9,
                arrowprops=dict(arrowstyle='-', color='#777777'))
    bx.text(.55, .93, '◇ 7×7×17 transverse refinement\n'
            f"tip {other['tip_vs_beam_pct']:.3f}% · moment "
            f"{other['maximum_midspan_moment_vs_beam_pct']:.3f}%",
            transform=bx.transAxes, fontsize=9, va='top', ha='left',
            bbox=dict(facecolor='white', alpha=.92, edgecolor='none', pad=3))
    ax.set(yscale='log', title='Midpoint section-shear resultant',
           ylabel='Absolute discrepancy vs 1D beam reference (%)')
    bx.set(title='Displacement and section moment',
           ylabel='Absolute discrepancy vs 1D beam reference (%)')
    for a in (ax, bx):
        a.set_xticks(xs, [str(v) for v in n])
        a.set(xlabel='Axial vertex count nz (5×5 transverse grid)')
        a.grid(True, axis='y', alpha=.22)
        a.spines[['top','right']].set_visible(False)
        a.set_xlim(-.18, len(xs)-.82)
    ax.set_ylim(0.0014, 48)
    bx.set_ylim(0, 3.1)
    for x, value in zip(xs, shear):
        ax.annotate(f'{value:.4g}%', (x, value), xytext=(0, 10),
                    textcoords='offset points', ha='center', fontsize=8)
    ax.legend(loc='upper right', frameon=False, fontsize=9)
    bx.legend(loc='upper left', frameon=False, fontsize=9)
    fig.suptitle('Cubic Hermite: measured mesh-density trends (saved solved results)',
                 fontsize=15, weight='semibold')
    fig.supxlabel('Open diamonds = separate 7×7×17 transverse refinement; not connected to the 5×5 axial series.\n'
                  'Section-midpoint integrals, not pointwise stress or certified 3D error.', fontsize=9)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'docs/figures/cubic_densification.svg')
    args = parser.parse_args(argv)
    print('Updated figure from measured JSON:', render(args.out))


if __name__ == '__main__':
    main()
