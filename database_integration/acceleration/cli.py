#!/usr/bin/env python3
"""
QCxMS2 Acceleration CLI

Command-line interface for fragment prioritization and cache lookup.

Usage:
    python -m database_integration.acceleration.cli prioritize --smiles "c1ccccc1"
    python -m database_integration.acceleration.cli prioritize --smiles "c1ccccc1" --output constraints.inp
    python -m database_integration.acceleration.cli prioritize --spectrum query.mgf --library library.mgf
"""

import argparse
import sys
from pathlib import Path


def cmd_prioritize(args: argparse.Namespace) -> int:
    """Prioritize fragments for a molecule."""
    from .fragment_prioritizer import FragmentPrioritizer
    
    prioritizer = FragmentPrioritizer()
    
    if args.smiles:
        print(f"Prioritizing fragments for: {args.smiles}")
        print()
        
        if args.spectrum_guided:
            # Use DreaMS with spectrum-based search
            fragments = prioritizer.from_smiles(
                args.smiles,
                top_k_fragments=args.top_k
            )
        else:
            # Use chemical rules (fast, no dependencies)
            fragments = prioritizer.from_common_rules(
                args.smiles,
                top_k=args.top_k
            )
    
    elif args.spectrum:
        print(f"Prioritizing fragments from spectrum: {args.spectrum}")
        print()
        
        fragments = prioritizer.from_spectrum(
            args.spectrum,
            library_path=args.library,
            top_k_fragments=args.top_k
        )
    
    else:
        print("Error: Must specify --smiles or --spectrum")
        return 1
    
    # Display results
    print(f"{'Rank':>4} {'m/z':>10} {'Formula':>10} {'Priority':>10} {'Annotation':<25}")
    print("-" * 65)
    
    for i, frag in enumerate(fragments, 1):
        formula = frag.formula or ''
        annotation = frag.annotation or ''
        print(f"{i:>4} {frag.mz:>10.4f} {formula:>10} {frag.priority_score:>10.3f} {annotation:<25}")
    
    print()
    
    # Write output files if requested
    if args.output:
        if args.output.endswith('.inp'):
            prioritizer.write_crest_constraints(
                fragments, 
                args.output,
                max_fragments=args.max_constraints
            )
        elif args.output.endswith('.md'):
            prioritizer.write_priority_report(fragments, args.output)
        elif args.output.endswith('.json'):
            import json
            with open(args.output, 'w') as f:
                json.dump([{
                    'mz': frag.mz,
                    'priority_score': frag.priority_score,
                    'formula': frag.formula,
                    'annotation': frag.annotation
                } for frag in fragments], f, indent=2)
            print(f"Wrote {len(fragments)} fragments to {args.output}")
    
    return 0


def cmd_cache_lookup(args: argparse.Namespace) -> int:
    """Look up cached fragments."""
    from ..cache import FragmentCache
    
    cache = FragmentCache(args.db)
    
    if args.mz:
        result = cache.find_fragment(mz=args.mz, formula=args.formula)
        if result:
            print(f"Found cached fragment:")
            print(f"  m/z:     {result.mz:.4f}")
            print(f"  Formula: {result.formula}")
            print(f"  Barrier: {result.barrier_kcal_mol:.2f} kcal/mol")
            return 0
        else:
            print(f"No cached fragment found for m/z={args.mz:.4f}")
            return 1
    
    elif args.formula:
        results = cache.find_fragments_by_formula(args.formula)
        print(f"Fragments with formula {args.formula}:")
        for frag in results[:20]:
            print(f"  m/z={frag.mz:.4f}, barrier={frag.barrier_kcal_mol:.2f} kcal/mol")
        return 0
    
    else:
        print("Error: Must specify --mz or --formula")
        return 1


def cmd_integrate(args: argparse.Namespace) -> int:
    """
    Combined workflow: prioritize + cache lookup.
    
    Prioritizes fragments using DreaMS/rules, then checks cache for
    pre-computed barriers.
    """
    from .fragment_prioritizer import FragmentPrioritizer
    from ..cache import FragmentCache
    
    prioritizer = FragmentPrioritizer()
    cache = FragmentCache(args.db)
    
    print(f"Integrated fragment analysis for: {args.smiles}")
    print()
    
    # Step 1: Prioritize fragments
    fragments = prioritizer.from_common_rules(args.smiles, top_k=args.top_k)
    
    # Step 2: Look up cached barriers
    print(f"{'Rank':>4} {'m/z':>10} {'Formula':>10} {'Priority':>10} {'Cached Barrier':>15} {'Status':<10}")
    print("-" * 75)
    
    cached_count = 0
    for i, frag in enumerate(fragments, 1):
        formula = frag.formula or ''
        
        # Look up in cache
        cached = cache.find_fragment(mz=frag.mz, formula=frag.formula if frag.formula else None, mz_tolerance=0.1)
        
        if cached:
            barrier_str = f"{cached.barrier_kcal_mol:.1f} kcal/mol"
            status = "CACHED"
            cached_count += 1
        else:
            barrier_str = "-"
            status = "COMPUTE"
        
        print(f"{i:>4} {frag.mz:>10.4f} {formula:>10} {frag.priority_score:>10.3f} {barrier_str:>15} {status:<10}")
    
    print()
    print(f"Cache hit rate: {cached_count}/{len(fragments)} ({100*cached_count/len(fragments):.1f}%)")
    
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="QCxMS2 Acceleration Tools",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # prioritize command
    prio_parser = subparsers.add_parser("prioritize", help="Prioritize fragments for a molecule")
    prio_parser.add_argument("--smiles", help="SMILES string of molecule")
    prio_parser.add_argument("--spectrum", help="Path to query spectrum (MGF)")
    prio_parser.add_argument("--library", help="Path to reference library (MGF)")
    prio_parser.add_argument("--spectrum-guided", action="store_true", 
                            help="Use DreaMS for structure-based search (requires database)")
    prio_parser.add_argument("--top-k", type=int, default=20, help="Number of fragments to return")
    prio_parser.add_argument("--output", "-o", help="Output file (.inp for CREST, .md for report, .json for data)")
    prio_parser.add_argument("--max-constraints", type=int, default=15, help="Max constraints in CREST file")
    
    # cache lookup command
    cache_parser = subparsers.add_parser("cache", help="Look up cached fragments")
    cache_parser.add_argument("--db", default=None, help="Cache database path")
    cache_parser.add_argument("--mz", type=float, help="m/z to look up")
    cache_parser.add_argument("--formula", help="Formula to search")
    
    # integrate command
    int_parser = subparsers.add_parser("integrate", help="Combined prioritize + cache lookup")
    int_parser.add_argument("--smiles", required=True, help="SMILES string of molecule")
    int_parser.add_argument("--db", default=None, help="Cache database path")
    int_parser.add_argument("--top-k", type=int, default=20, help="Number of fragments to analyze")
    
    args = parser.parse_args()
    
    if args.command == "prioritize":
        return cmd_prioritize(args)
    elif args.command == "cache":
        return cmd_cache_lookup(args)
    elif args.command == "integrate":
        return cmd_integrate(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
