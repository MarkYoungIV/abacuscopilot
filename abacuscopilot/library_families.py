"""Pseudopotential / orbital library "families".

A *family* is a named set of pseudopotential + numerical-orbital directories
that ABACUS should use for a calculation.  The bundled ``SG15`` (+ lanthanide)
libraries stay the default and are auto-detected from ``PP-Orb/``; additional
series such as ``ABACUS-APNS-PPORBs-v1`` are **not** bundled — the user points
to their own downloaded directories via ``config.yaml``
(``libraries.families.<id>.pseudo_dir`` etc.) and selects the series from the
interactive INPUT flow.

Known series are registered in :data:`KNOWN_FAMILIES`.  Adding a future series
only requires a new entry here (a label, its config path keys, and — when the
orbital files are split into sub-variants such as APNS ``efficiency`` /
``precision`` — the variant list).

The "active" choice is persisted to ``config['libraries']['family']``
(e.g. ``sg15`` / ``apns/efficiency`` / ``apns/precision``), and the matching
directories are materialized into the flat ``pseudo_library`` /
``orbital_library`` lists that every auto-copy / resolution step already reads.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from abacuscopilot.config import save_config
from abacuscopilot.core.standards import is_lcao_basis

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FamilyVariant:
    """One orbital sub-variant of a family (e.g. APNS 'efficiency')."""

    id: str
    label: str
    # Config path to the orbital dir, e.g.
    # ("libraries","families","apns","orbital_dirs","efficiency")
    orbital_key: tuple[str, ...]


@dataclass(frozen=True)
class FamilySpec:
    """A registered, user-configurable library family."""

    id: str
    label: str
    # Config path to the pseudopotential dir, e.g.
    # ("libraries","families","apns","pseudo_dir")
    pseudo_key: tuple[str, ...]
    variants: dict[str, FamilyVariant] = field(default_factory=dict)


# Extensible registry: future series (ONCV, PBE-dojo, ...) get an entry here.
KNOWN_FAMILIES: dict[str, FamilySpec] = {
    "apns": FamilySpec(
        id="apns",
        label="ABACUS-APNS-PPORBs-v1",
        pseudo_key=("libraries", "families", "apns", "pseudo_dir"),
        variants={
            "efficiency": FamilyVariant(
                id="efficiency",
                label="efficiency (relax / MD / quick forces)",
                orbital_key=("libraries", "families", "apns", "orbital_dirs", "efficiency"),
            ),
            "precision": FamilyVariant(
                id="precision",
                label="precision (high-accuracy: bands / TDDFT / GW)",
                orbital_key=("libraries", "families", "apns", "orbital_dirs", "precision"),
            ),
        },
    ),
    # Dojo NC-FR: PseudoDojo norm-conserving FULLY-RELATIVISTIC PPs (ONCVPSP,
    # relativistic="full", has_so) + their matching numerical orbitals — the set
    # for SOC (INPUT: lspinorb 1) or any fully-relativistic run.  The orbital
    # package ships one folder per element and tier, {El}_{SZ,DZP,TZDP}, all
    # three tiers under the same Orbitals_v2.0 root; the selected tier is what
    # matters (materialized to the shared root, resolution restricted by tier).
    "dojoncfr": FamilySpec(
        id="dojoncfr",
        label="Dojo-NC-FR (fully-relativistic — for SOC)",
        pseudo_key=("libraries", "families", "dojoncfr", "pseudo_dir"),
        variants={
            # dzp first: it is the recommended default (Enter/blank on a fresh
            # pick chooses the first listed variant).
            "dzp": FamilyVariant(
                id="dzp",
                label="DZP (double-zeta + polarization) — recommended",
                orbital_key=("libraries", "families", "dojoncfr", "orbital_dirs", "dzp"),
            ),
            "sz": FamilyVariant(
                id="sz",
                label="SZ (single-zeta) — fastest",
                orbital_key=("libraries", "families", "dojoncfr", "orbital_dirs", "sz"),
            ),
            "tzdp": FamilyVariant(
                id="tzdp",
                label="TZDP (triple-zeta + double polarization) — most accurate",
                orbital_key=("libraries", "families", "dojoncfr", "orbital_dirs", "tzdp"),
            ),
        },
    ),
}


# PP-Orb/ drop-in layout for each registered family: the top folder a user's
# own copy of the series lives under when dropped into the package PP-Orb/
# directory, plus the relative sub-paths to the pseudopotential and (per
# variant) orbital roots inside it.  Used to auto-register a series that is
# already physically present, instead of prompting for its paths by hand.
#   key: family id -> (top folder, pseudo sub-path, {variant: orbital sub-path})
# A shared orbital root (Dojo keeps SZ/DZP/TZDP all under one Orbitals_v2.0)
# simply repeats the same sub-path for every variant.
_PP_ORB_LAYOUT: dict[str, tuple[str, str, dict[str, str]]] = {
    "apns": (
        "ABACUS-APNS-PPORBs-v1",
        "apns-pseudopotentials-v1",
        {
            "efficiency": "apns-orbitals-efficiency-v1",
            "precision": "apns-orbitals-precision-v1",
        },
    ),
    "dojoncfr": (
        "Dojo-NC-FR",
        "Pseudopotential",
        {"sz": "Orbitals_v2.0", "dzp": "Orbitals_v2.0", "tzdp": "Orbitals_v2.0"},
    ),
}


def _pp_orb_root() -> Path:
    """The package-level ``PP-Orb/`` directory (bundled libs + dropped series).

    Kept as a function (not a module constant) so tests can monkeypatch it:
    the developer's real PP-Orb/ physically contains Dojo-NC-FR / APNS, which
    would otherwise leak into tests that must exercise the prompt path.
    """
    return Path(__file__).resolve().parent.parent / "PP-Orb"


def _auto_detect_pp_orb(spec: FamilySpec) -> dict[str, Any] | None:
    """Locate a registered family's folder under the package PP-Orb/.

    Returns ``{"pseudo_dir": str, "orbital_dirs": {variant: str}}`` for the
    dirs that physically exist (paths resolved), or ``None`` when the family's
    canonical top folder is not present under PP-Orb/.  Only missing config
    paths are ever filled from the result, so a manual registration elsewhere
    still wins.
    """
    layout = _PP_ORB_LAYOUT.get(spec.id)
    if layout is None:
        return None
    top, pseudo_sub, orb_subs = layout
    base = _pp_orb_root() / top
    if not base.is_dir():
        return None
    out: dict[str, Any] = {}
    pseudo = base / pseudo_sub
    if pseudo.is_dir():
        out["pseudo_dir"] = str(pseudo.resolve())
    orbital_dirs: dict[str, str] = {}
    for vid, sub in orb_subs.items():
        d = base / sub
        if d.is_dir():
            orbital_dirs[vid] = str(d.resolve())
    if orbital_dirs:
        out["orbital_dirs"] = orbital_dirs
    return out if (out.get("pseudo_dir") or out.get("orbital_dirs")) else None


# ---------------------------------------------------------------------------
# Small config accessors
# ---------------------------------------------------------------------------


def _nested(config: dict, key: tuple[str, ...]) -> Any:
    cur: Any = config
    for k in key:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(k)
    return cur


def current_family(config: dict) -> str:
    """The family id currently in effect (default ``sg15``)."""
    libs = config.get("libraries", {})
    return libs.get("family", "sg15") or "sg15"


def _family_root(config: dict) -> dict:
    libs = config.setdefault("libraries", {})
    fams = libs.setdefault("families", {})
    if not isinstance(fams, dict):
        fams = {}
        libs["families"] = fams
    return fams


def is_user_family(family: str) -> bool:
    """True when *family* refers to a registered, non-bundled series.

    Used to decide whether a missing element is a hard error (registered user
    families: never silently fall back / mix) vs. the legacy SG15 soft warning.
    """
    if not family or family == "sg15":
        return False
    return family.split("/", 1)[0] in KNOWN_FAMILIES


def family_label(family: str) -> str:
    """Human label for a family id, e.g. ``apns/precision`` -> APNS label + variant."""
    if not family or family == "sg15":
        return "SG15 (bundled default)"
    if family == "custom":
        return "custom (hand-configured lists)"
    head, _, variant = family.partition("/")
    spec = KNOWN_FAMILIES.get(head)
    if spec is None:
        return family
    if variant and variant in spec.variants:
        return f"{spec.label} — {spec.variants[variant].label}"
    return spec.label


def hint_no_library_dirs(console, config: dict) -> None:
    """Point the user at PP-Orb/ when no pseudopotential/orbital library is usable.

    Fires only when the active library lists are empty AND no external series is
    registered (so no family path could cover the structure either).  Guides the
    user to drop their own ``*.upf`` / ``*.orb`` library folders under
    ``PP-Orb/`` (each top-level folder there is auto-detected; see the shipped
    ``PP-Orb/README.md``), or to register an external series such as
    ``ABACUS-APNS-PPORBs-v1`` under ``libraries.families`` and pick it below.
    """
    from abacuscopilot.config import _registered_family_dirs

    libs = config.get("libraries", {})
    if libs.get("pseudo_library") or libs.get("orbital_library"):
        return
    if _registered_family_dirs(config):
        return  # a family is configured — the picker can activate it

    pp_orb = _pp_orb_root()
    console.print()
    console.print("  [yellow]No pseudopotential/orbital libraries are available yet.[/yellow]")
    console.print("    Drop library folders (containing *.upf / *.orb) under:")
    console.print(f"      [bold]{pp_orb}[/bold]")
    console.print("    Each top-level folder there is auto-detected (searched recursively) —")
    console.print("    see PP-Orb/README.md in the package root.  Or register an external")
    console.print("    series (e.g. ABACUS-APNS-PPORBs-v1) under libraries.families in the")
    console.print("    config and choose it from the list below.")
    console.print()


# ---------------------------------------------------------------------------
# Materialize: family id  ->  the two flat library lists (+ persist family key)
# ---------------------------------------------------------------------------


def materialize_family(config: dict, family_id: str, variant: str | None = None) -> bool:
    """Point ``pseudo_library``/``orbital_library`` at *family_id*'s directories.

    ``sg15`` re-derives the bundled auto-detected dirs.  A registered user
    family reads its configured paths from ``libraries.families``; *variant*
    selects an orbital sub-set when the family defines one (and LCAO needs it).

    Updates ``config`` in place (caller persists with save_config) and returns
    False (leaving the config untouched) if a required path is missing/invalid.
    """
    libs = config.setdefault("libraries", {})

    if family_id == "sg15":
        # Re-derive the bundled auto-detected libraries (PP-Orb under the pkg
        # root), excluding any directory the user registered as an external
        # family (e.g. an ABACUS-APNS copy dropped under PP-Orb/).
        from abacuscopilot.config import _detect_library_dirs, _registered_family_dirs

        exclude = _registered_family_dirs(config)
        libs["pseudo_library"] = _detect_library_dirs(".upf", exclude)
        libs["orbital_library"] = _detect_library_dirs(".orb", exclude)
        libs["family"] = "sg15"
        return True

    spec = KNOWN_FAMILIES.get(family_id)
    if spec is None:
        return False

    fam_root = _family_root(config).setdefault(family_id, {})
    pseudo_dir = str(fam_root.get("pseudo_dir") or "")
    pseudo_path = Path(pseudo_dir).expanduser() if pseudo_dir else None
    if pseudo_path is None or not pseudo_path.is_dir():
        return False

    # Validate every required path *before* mutating, so a failed materialize
    # leaves the previously active libraries untouched.
    orb_path: Path | None = None
    if variant is not None:
        v = spec.variants.get(variant)
        if v is None:
            return False
        orb_dir = str(_nested(config, v.orbital_key) or "")
        if not orb_dir:
            return False
        orb_path = Path(orb_dir).expanduser()
        if not orb_path.is_dir():
            return False

    libs["pseudo_library"] = [str(pseudo_path.resolve())]
    if variant is not None and orb_path is not None:
        libs["orbital_library"] = [str(orb_path.resolve())]
        libs["family"] = f"{family_id}/{variant}"
    else:
        libs["orbital_library"] = []
        libs["family"] = family_id
    return True


# ---------------------------------------------------------------------------
# Interactive picker (INPUT / full-calculation-setup flows)
# ---------------------------------------------------------------------------


def pick_library_family(
    console,
    config: dict,
    *,
    basis_type: str,
    interactive: bool = True,
) -> bool:
    """Offer to switch the pseudopotential/orbital library family.

    Runs only when *interactive*; otherwise the persisted family is used as-is
    (returns False).  Enter at the first prompt keeps the current family.

    The dialog mirrors the requested expansion:
      1. current family (default SG15) -> Enter keeps it;
      2. choosing a user series such as ABACUS-APNS-PPORBs-v1 reveals its
         orbital sub-variants (efficiency / precision) — asked when the basis
         needs orbitals (LCAO/lcao_in_pw) and no sub-variant is recorded yet;
      3. unconfigured paths are prompted for and validated on the spot (the
         family is listed as "paths not set yet" until configured).

    Returns True if the family changed (config was updated + saved).
    """
    if not interactive:
        return False

    need_orb = is_lcao_basis(basis_type)
    current = current_family(config)
    _print_header(console, family_label(current))
    # Fresh install with nothing usable under PP-Orb/ (or registered): tell the
    # user where to drop their own libraries before listing the choices.
    hint_no_library_dirs(console, config)

    options: list[tuple[str, str]] = [("sg15", "SG15 (bundled default)")]
    options += [
        (fid, _option_text(config, fid, spec)) for fid, spec in KNOWN_FAMILIES.items()
    ]

    for i, (fid, text) in enumerate(options, 1):
        mark = "  [dim](current)[/dim]" if fid == current else ""
        console.print(f"    [green]{i})[/green]  {text}{mark}")

    console.print("  [dim]Enter = keep current family (no change).[/dim]")
    ans = _ask_index(console, "Choose a library family", "keep current").strip()
    if not ans:
        return False
    try:
        choice_id = options[int(ans) - 1][0]
    except (ValueError, IndexError):
        console.print("  [red]![/red] Invalid choice — keeping current family.")
        return False

    if choice_id == "sg15":
        materialize_family(config, "sg15")
        save_config(config)
        console.print("  [green]✓ Library family set to SG15 (bundled).[/green]")
        console.print()
        return True

    spec = KNOWN_FAMILIES[choice_id]
    variant: str | None = None
    cur_variant = current.split("/", 1)[1] if current.startswith(f"{choice_id}/") else None
    if spec.variants and need_orb:
        # Always offer the orbital sets so the choice stays discoverable /
        # switchable.  The currently active one is marked and is the Enter
        # default, so re-picking the family can never silently flip
        # efficiency<->precision (an early bug had a fixed "1" default that
        # downgraded precision users on an accidental Enter).
        variant = _ask_variant(console, spec, default=cur_variant)
        if variant is None:
            return False               # user aborted
    elif cur_variant is not None:
        variant = cur_variant          # keep the current pin (e.g. PW flow)

    if not _ensure_family_paths(console, config, spec, variant):
        return False

    # Dojo tiers ship several cutoff radii of the SAME basis (identical zeta);
    # let the user pick which copy to use.  The persisted policy
    # (libraries.rcut_policy) drives every later resolution of this family.
    if choice_id == "dojoncfr":
        _ask_rcut_policy(console, config)

    if not materialize_family(config, choice_id, variant):
        console.print("  [red]![/red] Could not activate the selected family (paths invalid?).")
        return False

    save_config(config)
    libs = config.get("libraries", {})
    console.print(
        f"  [green]✓ Library family set to {family_label(libs.get('family', choice_id))}.[/green]"
    )
    console.print(f"    PP:     [dim]{libs.get('pseudo_library', [])}[/dim]")
    if libs.get("orbital_library"):
        console.print(f"    Orbitals: [dim]{libs.get('orbital_library')}[/dim]")
    if choice_id == "apns":
        console.print(
            "  [yellow]note:[/yellow] APNS PPs recommend ecutwfc up to ~150 Ry for some "
            "elements — raise ecutwfc accordingly or run task 107."
        )
    elif choice_id == "dojoncfr":
        from abacuscopilot.preprocessing.system_tasks import (
            RCUT_POLICY_MAX,
            RCUT_POLICY_MIN,
            _rcut_policy,
        )

        pol = _rcut_policy(config)
        pol_label = {
            RCUT_POLICY_MIN: "smallest rcut (fastest)",
            RCUT_POLICY_MAX: "largest rcut (most complete)",
        }.get(pol, "closest to 7 au (recommended)")
        console.print(f"    rcut default: [dim]{pol_label}[/dim] (config libraries.rcut_policy)")
        console.print(
            "  [yellow]note:[/yellow] Dojo-NC-FR are fully-relativistic norm-conserving "
            "PPs (relativistic=\"full\", has_so=1) — required for spin-orbit coupling "
            "(INPUT: lspinorb 1). They also work for non-SOC runs: ABACUS reduces them "
            "to scalar-relativistic automatically."
        )
    console.print()
    return True


def _option_text(config: dict, family_id: str, spec: FamilySpec) -> str:
    if _family_configured(config, family_id):
        return spec.label
    if _auto_detect_pp_orb(spec):
        return f"{spec.label}  (found in PP-Orb — will auto-register)"
    return f"{spec.label}  (paths not set yet — will prompt)"


def _family_configured(config: dict, family_id: str) -> bool:
    """True once the family has at least one usable directory configured."""
    fam = _family_root(config).get(family_id)
    if not isinstance(fam, dict):
        return False

    def _valid(v: Any) -> bool:
        return bool(isinstance(v, str) and v and Path(v).expanduser().is_dir())

    if _valid(fam.get("pseudo_dir")):
        return True
    orb = fam.get("orbital_dirs")
    if isinstance(orb, dict):
        return any(_valid(v) for v in orb.values())
    return False


def _print_header(console, current_label: str) -> None:
    console.print()
    console.print("[bold]Pseudopotential / orbital library family[/bold]")
    console.print(f"  Current: {current_label}")


def _ask_index(console, text: str, default_desc: str) -> str:
    """Ask for an option number; blank returns '' (keep current)."""
    return console.input(f"  {text} [{default_desc}]: ").strip()


def _ask_variant(console, spec: FamilySpec, default: str | None = None) -> str | None:
    """Ask which orbital sub-variant of *spec* to use.

    *default* names the variant to mark "(current)" and return on Enter/blank.
    """
    items = list(spec.variants.values())
    idx = {v.id: i for i, v in enumerate(items)}
    if default is not None and default not in idx:
        default = None
    console.print(f"  {spec.label} has these orbital sets:")
    for i, v in enumerate(items, 1):
        mark = "  [dim](current)[/dim]" if v.id == default else ""
        console.print(f"    [green]{i})[/green]  {v.label}{mark}")
    if default is None:
        default = items[0].id
    ans = _ask_index(console, "Which orbital set",
                     f"{idx[default] + 1} = {default}").strip()
    if not ans:
        return default
    try:
        return items[int(ans) - 1].id
    except (ValueError, IndexError):
        console.print("  [red]![/red] Invalid choice — aborted.")
        return None


def _ask_rcut_policy(console, config: dict) -> None:
    """Ask which cutoff-radius copy of a Dojo basis to use (persists on change).

    A Dojo-NC-FR tier ships the *same* zeta basis at several radii (6-12 au);
    this only disambiguates between those.  Enter/blank keeps the current
    policy (default "7" = closest to the bundled SG15 standard orbitals).
    Stores the choice under ``config['libraries']['rcut_policy']``.
    """
    from abacuscopilot.preprocessing.system_tasks import (
        RCUT_POLICY_DEFAULT,
        RCUT_POLICY_MAX,
        RCUT_POLICY_MIN,
        _rcut_policy,
    )

    current = _rcut_policy(config)
    options = [
        (RCUT_POLICY_DEFAULT,
         "closest to 7 au (recommended — matches bundled SG15 orbitals)"),
        (RCUT_POLICY_MIN, "smallest rcut (fastest / most compact)"),
        (RCUT_POLICY_MAX, "largest rcut (most complete / closest to PW limit)"),
    ]
    idx = {v: i for i, (v, _) in enumerate(options)}
    console.print()
    console.print("  Dojo-NC-FR ships each basis at several cutoff radii (6-12 au):")
    for i, (value, label) in enumerate(options, 1):
        mark = "  [dim](current)[/dim]" if value == current else ""
        console.print(f"    [green]{i})[/green]  {label}{mark}")
    ans = _ask_index(console, "Which rcut copy to use per element",
                     f"{idx[current] + 1} = {current}").strip()
    if not ans:
        return  # keep current
    try:
        value = options[int(ans) - 1][0]
    except (ValueError, IndexError):
        console.print("  [red]![/red] Invalid choice — keeping current rcut policy.")
        return
    if value != current:
        config.setdefault("libraries", {})["rcut_policy"] = value


def _ensure_family_paths(console, config: dict, spec: FamilySpec, variant: str | None) -> bool:
    """Register any family path that is missing — auto-detecting first.

    If the family's canonical folder already sits under PP-Orb/ (e.g. the user
    copied Dojo-NC-FR into PP-Orb/ on a fresh machine), its paths are filled in
    automatically and no prompt is shown.  Otherwise the missing path(s) are
    asked for and validated on the spot.  A stale/broken registered path also
    falls back to auto-detection / prompting.
    """
    root = _family_root(config).setdefault(spec.id, {})

    def _valid_dir(value: Any) -> bool:
        return bool(isinstance(value, str) and value
                    and Path(value).expanduser().is_dir())

    # 1. Auto-register every path the PP-Orb drop-in can provide.
    found = _auto_detect_pp_orb(spec)
    auto_registered = False
    if found:
        if not _valid_dir(root.get("pseudo_dir")) and found.get("pseudo_dir"):
            root["pseudo_dir"] = found["pseudo_dir"]
            auto_registered = True
        for v in spec.variants.values():
            key = v.orbital_key
            if not _valid_dir(_nested(config, key)) and found.get("orbital_dirs", {}).get(v.id):
                container = _nested_container(root, key)
                container[key[-1]] = found["orbital_dirs"][v.id]
                auto_registered = True
    if auto_registered:
        console.print()
        console.print("  [green]✓ " + spec.label
                      + " found under PP-Orb/ — paths auto-registered.[/green]")
        console.print()

    # 2. Whatever is still missing (or remains invalid) gets prompted.
    pseudo = str(root.get("pseudo_dir") or "")
    if (not pseudo or not Path(pseudo).expanduser().is_dir()) and not _ask_path(
        console, f"Path to {spec.id} pseudopotentials (contains *.upf)",
        root, "pseudo_dir",
    ):
        return False
    if not root.get("pseudo_dir"):
        return False

    if variant is not None:
        v = spec.variants[variant]
        cur = str(_nested(config, v.orbital_key) or "")
        if (not cur or not Path(cur).expanduser().is_dir()) and not _ask_path(
            console, f"Path to {spec.id} orbital basis — {v.label}",
            _nested_container(root, v.orbital_key), v.orbital_key[-1],
        ):
            return False
        if not _nested(config, v.orbital_key):
            return False
    return True


def _nested_container(root: dict, key: tuple[str, ...]) -> dict:
    """Locate/build the dict under the family root that owns the leaf of *key*.

    *key* is a full config path (e.g. (…,"apns","orbital_dirs","efficiency"));
    *root* is the family dict (…,"families","apns").  Walks the segments
    between the family id and the leaf, creating dicts as needed.
    """
    cur: Any = root
    for k in key[3:-1]:  # e.g. ("orbital_dirs",) — relative to the family root
        nxt = cur.get(k)
        if not isinstance(nxt, dict):
            nxt = {}
            cur[k] = nxt
        cur = nxt
    return cur


def _ask_path(console, question: str, container: dict, key: str) -> bool:
    """Prompt for an existing directory, storing it under container[key]."""
    for _ in range(3):
        console.print(f"  [yellow]{question}[/yellow]")
        ans = console.input("    Directory (or blank to abort): ").strip()
        if not ans:
            return False
        p = Path(ans).expanduser()
        if p.is_dir():
            container[key] = str(p.resolve())
            return True
        console.print(f"    [red]![/red] Not a directory: {p}")
    return False
