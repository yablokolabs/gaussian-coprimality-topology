# Verification

This file records release checks and representative outputs from the finished CLI. The
experiments are finite computations; they do not establish results about the full infinite
coprimality space.

## Quality gates

The following commands completed successfully with no warnings:

```bash
cargo fmt --all -- --check
cargo check --workspace --all-targets
cargo test --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
```

The test run executed 19 deterministic tests: 6 arithmetic tests, 9 factor/topology tests,
and 4 CLI integration tests. The arithmetic and topology suites include exhaustive bounded
property tests for norm multiplicativity, Euclidean remainder descent, factor reconstruction,
and the sigma intersection law.

## CLI demonstrations

### Exact gcd and factor/support

```console
$ cargo run --quiet -p gip-cli -- gcd 91 23 26 -7
gcd(91+23i, 26-7i) = 2+i
$ cargo run --quiet -p gip-cli -- factor 65 0
value: 65
norm: 4225
unit: -1
factorization: -1 * (1+2i) * (2+i) * (2+3i) * (3+2i)
support: {[1+2i], [2+i], [2+3i], [3+2i]}
```

### Indistinguishability

```console
$ cargo run --quiet -p gip-cli -- indistinguishable 2 0 1 1
indistinguishable(2, 1+i): true
support(2): {[1+i]}
support(1+i): {[1+i]}
```

### Specialization and quotient open

```console
$ cargo run --quiet -p gip-cli -- specializes 1 1 3 3
3+3i in closure({1+i}): true
Supp(1+i) subset of Supp(3+3i): true
$ cargo run --quiet -p gip-cli -- quotient-open 6 0 5 0
quotient point Supp(6): {[1+i], [3]}
forbidden F = Supp(5): {[1+2i], [2+i]}
point in O_F: true
```

### Mixed parity and bounded enumeration

```console
$ cargo run --quiet -p gip-cli -- basic-open 1 1 4 3
4+3i in sigma_1+i: true
$ cargo run --quiet -p gip-cli -- basic-open 1 1 4 2
4+2i in sigma_1+i: false
$ cargo run --quiet -p gip-cli -- window 3 1 1
window: re,im in [-3,3], excluding 0
index: 1+i
points: 48
sigma members: 24
sigma nonmembers: 24
distinct quotient points: 10
```
