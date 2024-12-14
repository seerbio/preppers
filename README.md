Preppers is a tool for *pr*otein to *pep*tide mapping, written in *R*u*s*t

## Usage

### As a command line tool

**Coming soon!**

For now, you can clone this repo and run:

```shell
cargo run --example annotate_fasta $PEPTIDES_FILE $FASTA_FILE
```

where `$PEPTIDES_FILE` is a file containing a list of peptides, one per line, and `$FASTA_FILE` is a FASTA file containing protein sequences.

### As a Python library

**Coming soon: CodeArtifact package**

For now, you can clone this repo and run:

```shell
pip install maturin
maturin develop
```

and then in your Python code:

```python
import preppers
```

### As a Rust library

Reference this repository in your `Cargo.toml`:

```toml
[dependencies]
preppers = { git = "ssh://git@github.com/seerbio/preppers.git" }
```
