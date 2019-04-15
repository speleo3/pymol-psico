## Run Tests

Run all tests, suppress deprecation warnings:

```bash
pytest -W ignore::DeprecationWarning
```

Skip tests which use the internet or external binaries:

```bash
pytest -W ignore::DeprecationWarning -m "not (web or exe)"
```

## Install Test Dependencies

```bash
pip install --editable .
conda install pytest biopython
conda install -c speleo3 csb
conda install -c salilab -c speleo3 dssp
conda install -c kimlab stride
conda install -c bioconda -c speleo3 msms
```
