Preppers is a tool for *pr*otein to *pep*tide mapping, written in *R*u*s*t

## Usage

### As a command line tool

**Coming soon!**

For now, you can clone this repo and run:

```shell
cargo run --example annotate_fasta $PEPTIDES_FILE $FASTA_FILE
```

where `$PEPTIDES_FILE` is a file containing a list of peptides, one per line, and `$FASTA_FILE` is a FASTA file containing protein sequences.

> [!NOTE]
> To set up a Rust environment, follow the instructions [here](https://www.rust-lang.org/tools/install).

### As a Python library

**Coming soon: CodeArtifact package**

For now, you can clone this repo and run:

```shell
maturin develop --release -m python/Cargo.toml
```

and then in your Python code:

```python
import preppers
```

> [!NOTE]
> To install `maturin`, create a virtual environment and run `pip install maturin`.

### As a Rust library

Reference this repository in your `Cargo.toml`:

```toml
[dependencies]
preppers = { git = "ssh://git@github.com/seerbio/preppers.git" }
```

## Running tests

This package includes tests for functionality in Rust and Python.

Rust tests:

```shell
cargo test
```

> [!NOTE]
> To set up a Rust environment, follow the instructions [here](https://www.rust-lang.org/tools/install).

Python tests:
```shell
maturin develop --release -m python/Cargo.toml
cd python
pytest
```