#!/usr/bin/env python3
"""
QCxMS2 Fragment Cache CLI

Command-line interface for managing the fragment cache.

Usage:
    python -m database_integration.cache.cli seed <directory> [--smiles SMILES]
    python -m database_integration.cache.cli lookup --mz 77.039 [--formula C6H5]
    python -m database_integration.cache.cli stats
    python -m database_integration.cache.cli export <output.json>
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from .fragment_cache import FragmentCache, get_default_cache


def cmd_seed(args: argparse.Namespace) -> int:
    """Seed cache from QCxMS2 output directories."""
    cache = FragmentCache(args.db)
    
    total_added = 0
    for directory in args.directories:
        path = Path(directory)
        if not path.exists():
            print(f"Warning: Directory not found: {directory}")
            continue
        
        # Check if it's a QCxMS2 output (has allfragments file)
        if (path / "allfragments").exists():
            added = cache.seed_from_qcxms2_output(
                str(path),
                parent_smiles=args.smiles,
                calculation_level=args.level
            )
            print(f"Added {added} fragments from {path.name}")
            total_added += added
        else:
            # Maybe it's a parent directory with multiple calculations
            for subdir in path.iterdir():
                if subdir.is_dir() and (subdir / "allfragments").exists():
                    added = cache.seed_from_qcxms2_output(
                        str(subdir),
                        parent_smiles=args.smiles,
                        calculation_level=args.level
                    )
                    print(f"Added {added} fragments from {subdir.name}")
                    total_added += added
    
    print(f"\nTotal fragments added: {total_added}")
    return 0


def cmd_lookup(args: argparse.Namespace) -> int:
    """Look up a fragment in the cache."""
    cache = FragmentCache(args.db)
    
    result = cache.find_fragment(
        mz=args.mz,
        formula=args.formula,
        mz_tolerance=args.tolerance
    )
    
    if result:
        print(f"Found cached fragment:")
        print(f"  m/z:         {result.mz:.4f}")
        print(f"  Formula:     {result.formula}")
        print(f"  Barrier:     {result.barrier_kcal_mol:.2f} kcal/mol ({result.barrier_ev:.4f} eV)")
        if result.de_kcal_mol:
            print(f"  ΔE:          {result.de_kcal_mol:.2f} kcal/mol")
        print(f"  Type:        {result.fragment_type}")
        print(f"  Atoms:       {result.atom_count}")
        print(f"  Charge:      {result.charge}+")
        if result.parent_formula:
            print(f"  Parent:      {result.parent_formula}")
        print(f"  Source:      {result.source_directory}")
        print(f"  Accesses:    {result.access_count}")
        return 0
    else:
        print(f"No cached fragment found for m/z={args.mz:.4f}", end="")
        if args.formula:
            print(f", formula={args.formula}", end="")
        print()
        return 1


def cmd_search(args: argparse.Namespace) -> int:
    """Search for fragments by various criteria."""
    cache = FragmentCache(args.db)
    
    if args.formula:
        results = cache.find_fragments_by_formula(args.formula, limit=args.limit)
        print(f"Fragments with formula {args.formula}:")
    elif args.mz_range:
        mz_min, mz_max = map(float, args.mz_range.split("-"))
        results = cache.find_fragments_by_mz_range(mz_min, mz_max, limit=args.limit)
        print(f"Fragments in m/z range {mz_min:.2f}-{mz_max:.2f}:")
    else:
        print("Error: Must specify --formula or --mz-range")
        return 1
    
    if not results:
        print("  (no results)")
        return 0
    
    print(f"\n{'m/z':>10} {'Formula':<10} {'Barrier (kcal/mol)':<20} {'Type':<15}")
    print("-" * 60)
    for frag in results:
        barrier_str = f"{frag.barrier_kcal_mol:.2f}" if frag.barrier_kcal_mol else "N/A"
        print(f"{frag.mz:>10.4f} {frag.formula:<10} {barrier_str:<20} {frag.fragment_type:<15}")
    
    print(f"\nTotal: {len(results)} fragments")
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    """Show cache statistics."""
    cache = FragmentCache(args.db)
    stats = cache.get_statistics()
    
    print("Fragment Cache Statistics")
    print("=" * 40)
    print(f"Database:           {cache.db_path}")
    print(f"Size:               {stats['db_size_mb']:.2f} MB")
    print(f"Total fragments:    {stats['total_fragments']}")
    print(f"Unique formulas:    {stats['unique_formulas']}")
    print(f"Unique m/z values:  {stats['unique_mz_values']}")
    print(f"With barriers:      {stats['fragments_with_barriers']}")
    print(f"Source calculations: {stats['source_calculations']}")
    
    if stats['most_accessed']:
        print("\nMost Accessed Fragments:")
        for item in stats['most_accessed']:
            print(f"  {item['formula']:<10} m/z={item['mz']:.4f}  ({item['access_count']} accesses)")
    
    return 0


def cmd_export(args: argparse.Namespace) -> int:
    """Export cache to JSON."""
    cache = FragmentCache(args.db)
    count = cache.export_to_json(args.output)
    print(f"Exported {count} fragments to {args.output}")
    return 0


def cmd_clear(args: argparse.Namespace) -> int:
    """Clear the cache."""
    if not args.confirm:
        print("This will delete all cached fragments!")
        response = input("Type 'yes' to confirm: ")
        if response.lower() != 'yes':
            print("Aborted.")
            return 1
    
    cache = FragmentCache(args.db)
    cache.clear(confirm=True)
    print("Cache cleared.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="QCxMS2 Fragment Cache Manager",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--db", 
        default=None,
        help="Path to cache database (default: /mnt/beegfs/software/qcxms2/fragment_cache.db)"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # seed command
    seed_parser = subparsers.add_parser("seed", help="Seed cache from QCxMS2 output")
    seed_parser.add_argument("directories", nargs="+", help="QCxMS2 output directories")
    seed_parser.add_argument("--smiles", help="Parent molecule SMILES")
    seed_parser.add_argument("--level", default="gfn2", help="Calculation level (default: gfn2)")
    
    # lookup command
    lookup_parser = subparsers.add_parser("lookup", help="Look up a fragment")
    lookup_parser.add_argument("--mz", type=float, required=True, help="Target m/z")
    lookup_parser.add_argument("--formula", help="Molecular formula")
    lookup_parser.add_argument("--tolerance", type=float, default=0.001, help="m/z tolerance")
    
    # search command
    search_parser = subparsers.add_parser("search", help="Search fragments")
    search_parser.add_argument("--formula", help="Search by formula")
    search_parser.add_argument("--mz-range", help="Search by m/z range (e.g., '50-100')")
    search_parser.add_argument("--limit", type=int, default=50, help="Max results")
    
    # stats command
    stats_parser = subparsers.add_parser("stats", help="Show cache statistics")
    
    # export command
    export_parser = subparsers.add_parser("export", help="Export cache to JSON")
    export_parser.add_argument("output", help="Output JSON file")
    
    # clear command
    clear_parser = subparsers.add_parser("clear", help="Clear the cache")
    clear_parser.add_argument("--confirm", action="store_true", help="Skip confirmation")
    
    args = parser.parse_args()
    
    if args.command == "seed":
        return cmd_seed(args)
    elif args.command == "lookup":
        return cmd_lookup(args)
    elif args.command == "search":
        return cmd_search(args)
    elif args.command == "stats":
        return cmd_stats(args)
    elif args.command == "export":
        return cmd_export(args)
    elif args.command == "clear":
        return cmd_clear(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
