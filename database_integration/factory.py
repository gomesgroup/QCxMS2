"""
Factory functions for creating MS database connectors.
"""

from .base import MSDatabaseConnector
from .config import DatabaseConfig
from .connectors import (
    MassBankConnector,
    MassSpecGymConnector,
    MoNAConnector,
    NISTMSConnector,
    SDBSConnector,
)
from .enums import DatabaseType


def create_ms_connector(
    db_type: DatabaseType, **config_params
) -> MSDatabaseConnector:
    """
    Create and return appropriate MS database connector based on type.

    Args:
        db_type: Type of database (MASSBANK, MONA, NIST_MS, SDBS, etc.)
        **config_params: Additional configuration parameters

    Returns:
        Initialized MS database connector

    Raises:
        ValueError: If database type is not supported

    Example:
        >>> connector = create_ms_connector(DatabaseType.MASSBANK)
        >>> criteria = SpectrumSearchCriteria(formula="C6H12O6", ionization_mode=IonizationMode.EI)
        >>> results = connector.search_spectrum(criteria)
    """
    config = DatabaseConfig(database_type=db_type, name=db_type.value, **config_params)

    if db_type == DatabaseType.MASSBANK:
        return MassBankConnector(config)
    elif db_type == DatabaseType.MONA:
        return MoNAConnector(config)
    elif db_type == DatabaseType.NIST_MS:
        return NISTMSConnector(config)
    elif db_type == DatabaseType.SDBS:
        return SDBSConnector(config)
    elif db_type == DatabaseType.MASSSPECGYM:
        return MassSpecGymConnector(config)
    else:
        raise ValueError(f"Unsupported database type: {db_type}")


def create_massbank_connector(**config_params) -> MassBankConnector:
    """
    Create MassBank connector with optional custom configuration.

    Args:
        **config_params: Configuration parameters

    Returns:
        Initialized MassBank connector
    """
    config = DatabaseConfig.for_massbank(**config_params)
    return MassBankConnector(config)


def create_mona_connector(**config_params) -> MoNAConnector:
    """
    Create MoNA connector with optional custom configuration.

    Args:
        **config_params: Configuration parameters

    Returns:
        Initialized MoNA connector
    """
    config = DatabaseConfig.for_mona(**config_params)
    return MoNAConnector(config)


def create_nist_connector(**config_params) -> NISTMSConnector:
    """
    Create NIST WebBook connector (stub).

    Args:
        **config_params: Configuration parameters

    Returns:
        Initialized NIST connector
    """
    config = DatabaseConfig.for_nist_ms(**config_params)
    return NISTMSConnector(config)


def create_sdbs_connector(**config_params) -> SDBSConnector:
    """
    Create SDBS connector (stub).

    Args:
        **config_params: Configuration parameters

    Returns:
        Initialized SDBS connector
    """
    config = DatabaseConfig.for_sdbs(**config_params)
    return SDBSConnector(config)


def create_massspecgym_connector(**config_params) -> MassSpecGymConnector:
    """
    Create MassSpecGym connector for NeurIPS 2024 benchmark dataset.

    MassSpecGym contains 231K curated MS/MS spectra from 29K molecules.
    The dataset is downloaded from Hugging Face on first access (~2GB).

    Args:
        **config_params: Configuration parameters
            - splits: List of splits to load (default: ["train", "val", "test"])

    Returns:
        Initialized MassSpecGym connector

    Example:
        >>> connector = create_massspecgym_connector()
        >>> stats = connector.get_dataset_statistics()
        >>> print(f"Total spectra: {stats['total_spectra']}")
    """
    config = DatabaseConfig.for_massspecgym(**config_params)
    return MassSpecGymConnector(config)
