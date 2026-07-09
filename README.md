<img alt="preppers logo" src="./static/preppers-logo.jpeg" height="128" align="left" style="margin: 8px">

**Preppers** is a tool for *pr*otein to *pep*tide mapping, written in *R*u*s*t.
It is designed for extreme speed and efficiency, using the *adaptive radix tree* implementation from [`blart`](https://github.com/declanvk/blart).

## Usage

### As a Python library

For ease of use, Preppers includes Python bindings which can be installed
with `pip`:

To install, run:

```shell
pip install preppers
```

and then in your Python code:

```python
import preppers
```

### As a command line tool

**Coming soon!**

For now, you can clone this repo and run:

```shell
cargo run --example annotate_fasta $PEPTIDES_FILE $FASTA_FILE
```

where `$PEPTIDES_FILE` is a file containing a list of peptides, one per line, and `$FASTA_FILE` is a FASTA file containing protein sequences.

> [!TIP]
> To set up a Rust environment, follow the instructions [here](https://www.rust-lang.org/tools/install).

### As a Rust library

To use Preppers in a Rust project, reference this repository in your `Cargo.toml`:

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

Python tests:
```shell
maturin develop --release -m python/Cargo.toml
cd python
pytest
```

> [!TIP]
> For development, you will need to first install a Rust environment by following the instructions [here](https://www.rust-lang.org/tools/install).
>
> To install `maturin`, create a virtual environment and run `pip install maturin`.

## Releasing a new version

1. Update the version in `Cargo.toml` and `python/Cargo.toml`; these should match! Be sure to update the dependency version as well!

> [!IMPORTANT]
> Be sure to pick an appropriate [semantic version](https://semver.org/) for the new release!

2. Commit the changes and push/merge to the `main` branch (this may require a PR).
3. Create a new release in GitHub with the same version number using the format `vX.Y.Z`.
   Click "auto-generate release notes" to automatically document merged PRs.
4. GitHub Actions will automatically build and publish Python packages.
