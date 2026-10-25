# Publish staging branch from Windows

The GitHub staging branch already exists but **does not yet contain the release source**. GitHub `main` is unchanged.

Extract `HermiteStructural_GitHub_Release_20261024.zip` to a clean folder. From that folder:

```powershell
python -m unittest discover -s tests -v
python .\tools\verify_source_provenance.py
cd .\orthotropic_companion
python -m unittest discover -s tests -v
cd ..
.\Publish-Tutorial001.ps1
```

If PowerShell execution policy blocks it, run from a PowerShell prompt:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\Publish-Tutorial001.ps1
```

The script clones ONLY the dedicated staging branch, shows the staged diff, and requires the literal `PUBLISH` before committing and pushing. It **does not** touch `main`. It uses Git for Windows and existing GitHub credentials.

After pushing, review the linked compare page and open a pull request. **Do not merge before full-size Hermite solver and PyVista tests on Windows**. The historical Hermite result is calibrated and not a rigorous a posteriori certificate.
