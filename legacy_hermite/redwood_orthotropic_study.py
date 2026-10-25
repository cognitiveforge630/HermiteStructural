"""Compatibility entry point for the existing Redwood study in backups/.

This forwards unchanged argument handling; the original source remains intact.
"""
from backups.redwood_orthotropic_study import main

if __name__ == '__main__':
    main()
