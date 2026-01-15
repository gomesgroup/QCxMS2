"""
Fragment Cache for QCxMS2

SQLite-based cache for storing computed fragments and their NEB barriers.
Provides fast m/z-based lookup and structure matching for accelerating
repeated QCxMS2 calculations.

Features:
- m/z indexing with configurable precision (default: 4 decimal places)
- Formula-based lookup
- Structure fingerprint matching
- Automatic seeding from QCxMS2 output directories
- Statistics and cache management
"""

import sqlite3
import json
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple, Union

from .qcxms2_output_parser import (
    QCxMS2OutputParser,
    QCxMS2Output,
    FragmentData
)


@dataclass
class FragmentEntry:
    """A cached fragment entry."""
    
    # Primary identification
    id: Optional[int] = None
    mz: float = 0.0
    mz_rounded: float = 0.0  # Rounded to cache precision
    formula: str = ""
    
    # Energetics
    barrier_kcal_mol: Optional[float] = None
    barrier_ev: Optional[float] = None
    de_kcal_mol: Optional[float] = None
    sumreac_kcal_mol: Optional[float] = None
    
    # Structure
    xyz_content: Optional[str] = None
    atom_count: int = 0
    structure_hash: Optional[str] = None
    
    # Fragment type and metadata
    fragment_type: str = "unknown"
    charge: int = 1
    multiplicity: int = 1
    
    # Parent molecule info
    parent_formula: Optional[str] = None
    parent_mz: Optional[float] = None
    parent_smiles: Optional[str] = None
    
    # Calculation info
    calculation_level: str = "gfn2"
    source_directory: Optional[str] = None
    source_type: str = "qcxms2"  # "qcxms2", "literature", "ml_predicted"
    
    # Timestamps
    created_at: Optional[str] = None
    accessed_at: Optional[str] = None
    access_count: int = 0
    
    @property
    def barrier_ev_computed(self) -> Optional[float]:
        """Get barrier in eV (computed from kcal/mol if needed)."""
        if self.barrier_ev is not None:
            return self.barrier_ev
        if self.barrier_kcal_mol is not None:
            return self.barrier_kcal_mol * 0.0433641  # kcal/mol to eV
        return None


class FragmentCache:
    """
    SQLite-based cache for QCxMS2 fragments.
    
    Provides fast lookup of pre-computed fragment barriers by m/z and formula,
    enabling significant speedups for repeated calculations.
    
    Example:
        cache = FragmentCache("/path/to/cache.db")
        
        # Seed from existing calculations
        cache.seed_from_qcxms2_output("/path/to/benzene_calculation")
        
        # Look up cached barrier
        result = cache.find_fragment(mz=77.039, formula="C6H5")
        if result:
            print(f"Using cached barrier: {result.barrier_kcal_mol} kcal/mol")
        else:
            # Run NEB calculation
            ...
    """
    
    # Conversion factors
    KCAL_TO_EV = 0.0433641
    EV_TO_KCAL = 23.0605
    
    def __init__(
        self,
        db_path: str = None,
        mz_precision: int = 4,
        auto_create: bool = True
    ):
        """
        Initialize fragment cache.
        
        Args:
            db_path: Path to SQLite database file. Default: auto-generated
            mz_precision: Decimal places for m/z rounding (default: 4)
            auto_create: Create database if it doesn't exist
        """
        if db_path is None:
            db_path = self._default_db_path()
        
        self.db_path = Path(db_path)
        self.mz_precision = mz_precision
        
        if auto_create:
            self._init_database()
    
    @staticmethod
    def _default_db_path() -> str:
        """Get default database path."""
        return "/mnt/beegfs/software/qcxms2/fragment_cache.db"
    
    def _init_database(self) -> None:
        """Initialize database schema."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        with sqlite3.connect(str(self.db_path)) as conn:
            conn.executescript("""
                -- Main fragments table
                CREATE TABLE IF NOT EXISTS fragments (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    
                    -- Mass/formula identification
                    mz REAL NOT NULL,
                    mz_rounded REAL NOT NULL,
                    formula TEXT NOT NULL,
                    
                    -- Energetics
                    barrier_kcal_mol REAL,
                    barrier_ev REAL,
                    de_kcal_mol REAL,
                    sumreac_kcal_mol REAL,
                    
                    -- Structure
                    xyz_content TEXT,
                    atom_count INTEGER,
                    structure_hash TEXT,
                    
                    -- Fragment type
                    fragment_type TEXT,
                    charge INTEGER DEFAULT 1,
                    multiplicity INTEGER DEFAULT 1,
                    
                    -- Parent molecule
                    parent_formula TEXT,
                    parent_mz REAL,
                    parent_smiles TEXT,
                    
                    -- Calculation info
                    calculation_level TEXT DEFAULT 'gfn2',
                    source_directory TEXT,
                    source_type TEXT DEFAULT 'qcxms2',
                    
                    -- Usage tracking
                    created_at TEXT,
                    accessed_at TEXT,
                    access_count INTEGER DEFAULT 0
                );
                
                -- Indexes for fast lookup
                CREATE INDEX IF NOT EXISTS idx_mz_rounded 
                    ON fragments(mz_rounded);
                CREATE INDEX IF NOT EXISTS idx_formula 
                    ON fragments(formula);
                CREATE INDEX IF NOT EXISTS idx_mz_formula 
                    ON fragments(mz_rounded, formula);
                CREATE INDEX IF NOT EXISTS idx_structure_hash 
                    ON fragments(structure_hash);
                CREATE INDEX IF NOT EXISTS idx_parent_formula 
                    ON fragments(parent_formula);
                
                -- Metadata table
                CREATE TABLE IF NOT EXISTS cache_metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT,
                    updated_at TEXT
                );
                
                -- Sources table (track seeded calculations)
                CREATE TABLE IF NOT EXISTS sources (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    source_path TEXT UNIQUE,
                    parent_formula TEXT,
                    fragment_count INTEGER,
                    calculation_type TEXT,
                    seeded_at TEXT
                );
            """)
            
            # Set metadata
            self._set_metadata(conn, "mz_precision", str(self.mz_precision))
            self._set_metadata(conn, "created_at", datetime.now().isoformat())
    
    def _set_metadata(self, conn: sqlite3.Connection, key: str, value: str) -> None:
        """Set a metadata value."""
        conn.execute(
            """
            INSERT OR REPLACE INTO cache_metadata (key, value, updated_at)
            VALUES (?, ?, ?)
            """,
            (key, value, datetime.now().isoformat())
        )
    
    def _round_mz(self, mz: float) -> float:
        """Round m/z to cache precision."""
        return round(mz, self.mz_precision)
    
    def add_fragment(self, fragment: FragmentEntry) -> int:
        """
        Add a fragment to the cache.
        
        Args:
            fragment: FragmentEntry to add
            
        Returns:
            ID of inserted fragment
        """
        # Compute derived values
        fragment.mz_rounded = self._round_mz(fragment.mz)
        if fragment.barrier_kcal_mol is not None and fragment.barrier_ev is None:
            fragment.barrier_ev = fragment.barrier_kcal_mol * self.KCAL_TO_EV
        if fragment.created_at is None:
            fragment.created_at = datetime.now().isoformat()
        
        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.execute(
                """
                INSERT INTO fragments (
                    mz, mz_rounded, formula,
                    barrier_kcal_mol, barrier_ev, de_kcal_mol, sumreac_kcal_mol,
                    xyz_content, atom_count, structure_hash,
                    fragment_type, charge, multiplicity,
                    parent_formula, parent_mz, parent_smiles,
                    calculation_level, source_directory, source_type,
                    created_at, accessed_at, access_count
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    fragment.mz, fragment.mz_rounded, fragment.formula,
                    fragment.barrier_kcal_mol, fragment.barrier_ev,
                    fragment.de_kcal_mol, fragment.sumreac_kcal_mol,
                    fragment.xyz_content, fragment.atom_count, fragment.structure_hash,
                    fragment.fragment_type, fragment.charge, fragment.multiplicity,
                    fragment.parent_formula, fragment.parent_mz, fragment.parent_smiles,
                    fragment.calculation_level, fragment.source_directory, fragment.source_type,
                    fragment.created_at, fragment.accessed_at, fragment.access_count
                )
            )
            return cursor.lastrowid
    
    def find_fragment(
        self,
        mz: float,
        formula: Optional[str] = None,
        structure_hash: Optional[str] = None,
        mz_tolerance: float = 0.001
    ) -> Optional[FragmentEntry]:
        """
        Find a cached fragment by m/z and optional criteria.
        
        Args:
            mz: Target m/z value
            formula: Optional molecular formula filter
            structure_hash: Optional structure hash for exact matching
            mz_tolerance: Tolerance for m/z matching (default: 0.001 Da)
            
        Returns:
            FragmentEntry if found, None otherwise
        """
        mz_rounded = self._round_mz(mz)
        mz_min = mz_rounded - mz_tolerance
        mz_max = mz_rounded + mz_tolerance
        
        query = "SELECT * FROM fragments WHERE mz_rounded BETWEEN ? AND ?"
        params: List[Any] = [mz_min, mz_max]
        
        if formula:
            query += " AND formula = ?"
            params.append(formula)
        
        if structure_hash:
            query += " AND structure_hash = ?"
            params.append(structure_hash)
        
        # Order by barrier (prefer lower barriers) and access count
        query += " ORDER BY barrier_kcal_mol ASC, access_count DESC LIMIT 1"
        
        with sqlite3.connect(str(self.db_path)) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(query, params)
            row = cursor.fetchone()
            
            if row:
                # Update access tracking
                conn.execute(
                    """
                    UPDATE fragments 
                    SET accessed_at = ?, access_count = access_count + 1
                    WHERE id = ?
                    """,
                    (datetime.now().isoformat(), row["id"])
                )
                return self._row_to_entry(row)
        
        return None
    
    def find_fragments_by_formula(
        self,
        formula: str,
        limit: int = 100
    ) -> List[FragmentEntry]:
        """
        Find all cached fragments with a given formula.
        
        Args:
            formula: Molecular formula to search for
            limit: Maximum results to return
            
        Returns:
            List of matching FragmentEntry objects
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                """
                SELECT * FROM fragments 
                WHERE formula = ?
                ORDER BY barrier_kcal_mol ASC
                LIMIT ?
                """,
                (formula, limit)
            )
            return [self._row_to_entry(row) for row in cursor.fetchall()]
    
    def find_fragments_by_mz_range(
        self,
        mz_min: float,
        mz_max: float,
        limit: int = 100
    ) -> List[FragmentEntry]:
        """
        Find fragments within an m/z range.
        
        Args:
            mz_min: Minimum m/z
            mz_max: Maximum m/z
            limit: Maximum results
            
        Returns:
            List of matching FragmentEntry objects
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                """
                SELECT * FROM fragments 
                WHERE mz_rounded BETWEEN ? AND ?
                ORDER BY mz_rounded ASC
                LIMIT ?
                """,
                (mz_min, mz_max, limit)
            )
            return [self._row_to_entry(row) for row in cursor.fetchall()]
    
    def _row_to_entry(self, row: sqlite3.Row) -> FragmentEntry:
        """Convert database row to FragmentEntry."""
        return FragmentEntry(
            id=row["id"],
            mz=row["mz"],
            mz_rounded=row["mz_rounded"],
            formula=row["formula"],
            barrier_kcal_mol=row["barrier_kcal_mol"],
            barrier_ev=row["barrier_ev"],
            de_kcal_mol=row["de_kcal_mol"],
            sumreac_kcal_mol=row["sumreac_kcal_mol"],
            xyz_content=row["xyz_content"],
            atom_count=row["atom_count"],
            structure_hash=row["structure_hash"],
            fragment_type=row["fragment_type"],
            charge=row["charge"],
            multiplicity=row["multiplicity"],
            parent_formula=row["parent_formula"],
            parent_mz=row["parent_mz"],
            parent_smiles=row["parent_smiles"],
            calculation_level=row["calculation_level"],
            source_directory=row["source_directory"],
            source_type=row["source_type"],
            created_at=row["created_at"],
            accessed_at=row["accessed_at"],
            access_count=row["access_count"]
        )
    
    def seed_from_qcxms2_output(
        self,
        output_dir: str,
        parent_smiles: Optional[str] = None,
        calculation_level: str = "gfn2"
    ) -> int:
        """
        Seed cache from a QCxMS2 calculation output directory.
        
        Args:
            output_dir: Path to QCxMS2 output directory
            parent_smiles: Optional SMILES of parent molecule
            calculation_level: QM level used (default: "gfn2")
            
        Returns:
            Number of fragments added
        """
        output_path = Path(output_dir)
        
        # Check if already seeded
        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.execute(
                "SELECT id FROM sources WHERE source_path = ?",
                (str(output_path.resolve()),)
            )
            if cursor.fetchone():
                return 0  # Already seeded
        
        # Parse output
        parser = QCxMS2OutputParser(output_dir)
        output = parser.parse()
        
        # Add fragments
        added = 0
        for frag_data in output.fragments:
            entry = FragmentEntry(
                mz=frag_data.mz,
                formula=frag_data.formula,
                barrier_kcal_mol=frag_data.barrier_kcal_mol,
                de_kcal_mol=frag_data.de_kcal_mol,
                sumreac_kcal_mol=frag_data.sumreac_kcal_mol,
                xyz_content=frag_data.xyz_content,
                atom_count=frag_data.atom_count,
                structure_hash=frag_data.structure_hash,
                fragment_type=frag_data.fragment_type,
                charge=frag_data.charge,
                multiplicity=frag_data.multiplicity,
                parent_formula=frag_data.parent_formula or output.parent_formula,
                parent_mz=frag_data.parent_mz or output.parent_mz,
                parent_smiles=parent_smiles or output.parent_smiles,
                calculation_level=calculation_level,
                source_directory=frag_data.directory,
                source_type="qcxms2"
            )
            self.add_fragment(entry)
            added += 1
        
        # Record source
        with sqlite3.connect(str(self.db_path)) as conn:
            conn.execute(
                """
                INSERT INTO sources (source_path, parent_formula, fragment_count, 
                                    calculation_type, seeded_at)
                VALUES (?, ?, ?, ?, ?)
                """,
                (
                    str(output_path.resolve()),
                    output.parent_formula,
                    added,
                    output.calculation_type,
                    datetime.now().isoformat()
                )
            )
        
        return added
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get cache statistics.
        
        Returns:
            Dict with cache statistics
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            stats = {}
            
            # Total fragments
            cursor = conn.execute("SELECT COUNT(*) FROM fragments")
            stats["total_fragments"] = cursor.fetchone()[0]
            
            # Unique formulas
            cursor = conn.execute("SELECT COUNT(DISTINCT formula) FROM fragments")
            stats["unique_formulas"] = cursor.fetchone()[0]
            
            # Unique m/z values (at cache precision)
            cursor = conn.execute("SELECT COUNT(DISTINCT mz_rounded) FROM fragments")
            stats["unique_mz_values"] = cursor.fetchone()[0]
            
            # Fragments with barriers
            cursor = conn.execute(
                "SELECT COUNT(*) FROM fragments WHERE barrier_kcal_mol IS NOT NULL"
            )
            stats["fragments_with_barriers"] = cursor.fetchone()[0]
            
            # Sources
            cursor = conn.execute("SELECT COUNT(*) FROM sources")
            stats["source_calculations"] = cursor.fetchone()[0]
            
            # Most accessed
            cursor = conn.execute(
                """
                SELECT formula, mz_rounded, access_count 
                FROM fragments 
                ORDER BY access_count DESC 
                LIMIT 5
                """
            )
            stats["most_accessed"] = [
                {"formula": r[0], "mz": r[1], "access_count": r[2]}
                for r in cursor.fetchall()
            ]
            
            # Database size
            stats["db_size_mb"] = self.db_path.stat().st_size / (1024 * 1024)
            
            return stats
    
    def export_to_json(self, output_path: str) -> int:
        """
        Export cache to JSON file.
        
        Args:
            output_path: Path for JSON output
            
        Returns:
            Number of fragments exported
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute("SELECT * FROM fragments")
            
            fragments = []
            for row in cursor.fetchall():
                entry = self._row_to_entry(row)
                # Convert to dict, excluding None values and large content
                d = asdict(entry)
                # Exclude large xyz_content from export
                if d.get("xyz_content"):
                    d["has_structure"] = True
                    del d["xyz_content"]
                fragments.append(d)
            
            with open(output_path, 'w') as f:
                json.dump({
                    "exported_at": datetime.now().isoformat(),
                    "fragment_count": len(fragments),
                    "fragments": fragments
                }, f, indent=2)
            
            return len(fragments)
    
    def clear(self, confirm: bool = False) -> None:
        """
        Clear all cached fragments.
        
        Args:
            confirm: Must be True to actually clear
        """
        if not confirm:
            raise ValueError("Must pass confirm=True to clear cache")
        
        with sqlite3.connect(str(self.db_path)) as conn:
            conn.execute("DELETE FROM fragments")
            conn.execute("DELETE FROM sources")
            self._set_metadata(conn, "cleared_at", datetime.now().isoformat())


def get_default_cache() -> FragmentCache:
    """Get the default fragment cache instance."""
    return FragmentCache()
