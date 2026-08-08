# Gaussian coprimality topology

This Rust workspace provides exact Gaussian-integer arithmetic and finite experiments with
the coprimality topology on the nonzero Gaussian integers `Z[i]`. It is useful for checking gcds and
factorizations, comparing prime supports, exploring topological indistinguishability and
specialization, and enumerating reproducible bounded lattice windows.

The workspace contains:

- `gaussian-integers`, a dependency-free reusable library;
- `gip-cli`, a dependency-free command-line interface.

All coordinate arithmetic is checked. Norms are exact `u128` values, and invalid domains,
unsupported factorization sizes, division by zero, and unrepresentable results are reported
as structured errors.

## Mathematical model

For a nonzero Gaussian integer `alpha`, let `Supp(alpha)` be the finite set of associate
classes of Gaussian primes dividing it. The basic open indexed by `alpha` is

```text
sigma_alpha = { beta != 0 : gcd(alpha, beta) is a unit }.
```

Equivalently, `beta` belongs to `sigma_alpha` exactly when the two supports are disjoint.
This gives the identities and classifications implemented here:

- `sigma_alpha ∩ sigma_beta = sigma_(alpha*beta)`;
- two points are topologically indistinguishable exactly when their supports are equal;
- `y` is in the closure of `{x}` exactly when `Supp(x) ⊆ Supp(y)`;
- a Kolmogorov-quotient point is a finite prime-class set `S`;
- its basic open `O_F` contains `S` exactly when `S ∩ F` is empty.

The canonical representative of a nonzero associate class is the unique rotation with
positive real part and nonnegative imaginary part. Thus `[1+i]`, for example, denotes all
four associates of `1+i`.

These definitions and examples follow Souvik Mandal,
[“A Coprimality Topology on the Gaussian Integers: Kolmogorov Quotient and Gaussian Prime Density”](https://arxiv.org/abs/2608.06373v1),
arXiv:2608.06373v1 (6 August 2026).

This software checks finite instances exactly. Bounded enumeration does **not** prove the
paper's results about the infinite space, density, or infinitude of Gaussian primes.

## Build and test

Stable Rust 1.85 or newer is supported.

```bash
cargo build --workspace
cargo test --workspace --all-targets
```

There are no third-party Rust dependencies.

## Library quick start

```rust
use gaussian_integers::{
    factor, prime_support, specializes_to, GaussianError, GaussianInt,
    KolmogorovPoint, QuotientBasicOpen,
};

fn main() -> Result<(), GaussianError> {
    let a = GaussianInt::new(91, 23);
    let b = GaussianInt::new(26, -7);
    assert_eq!(a.gcd(b)?, GaussianInt::new(2, 1));

    let value = GaussianInt::new(65, 0);
    let factors = factor(value)?;
    assert_eq!(factors.reconstruct()?, value);
    assert_eq!(factors.support(), prime_support(value)?);

    let source = GaussianInt::new(1, 1);
    let target = GaussianInt::new(3, 3);
    assert!(specializes_to(source, target)?);

    let point = KolmogorovPoint::from_gaussian(GaussianInt::new(6, 0))?;
    let forbidden = prime_support(GaussianInt::new(5, 0))?;
    assert!(QuotientBasicOpen::new(forbidden).contains(&point));
    Ok(())
}
```

The crate also exposes checked addition, subtraction, multiplication, negation,
conjugation, associates, Euclidean division/remainder, exact divisibility, normalized gcd,
Gaussian-prime validation, the sigma intersection predicate, and closure membership.

## CLI quick start

Coordinates are passed as signed decimal `RE IM` pairs. Discover all commands with:

```bash
cargo run --quiet -p gip-cli -- --help
```

The commands are `gcd`, `factor`, `support`, `coprime`, `basic-open`,
`indistinguishable`, `specializes`, `quotient-open`, and `window`.

## Practical use cases

### 1. Exact gcd, factorization, and support analysis

Compute a normalized Gaussian gcd:

```console
$ cargo run --quiet -p gip-cli -- gcd 91 23 26 -7
gcd(91+23i, 26-7i) = 2+i
```

Inspect the unit, canonical prime factors, and support of `65`:

```console
$ cargo run --quiet -p gip-cli -- factor 65 0
value: 65
norm: 4225
unit: -1
factorization: -1 * (1+2i) * (2+i) * (2+3i) * (3+2i)
support: {[1+2i], [2+i], [2+3i], [3+2i]}
```

### 2. Classify topologically indistinguishable integers

The paper's `2` versus `1+i` example is detected by support equality. Multiplicity is
discarded: `2` is a unit multiple of `(1+i)^2`, but both supports contain one class.

```console
$ cargo run --quiet -p gip-cli -- indistinguishable 2 0 1 1
indistinguishable(2, 1+i): true
support(2): {[1+i]}
support(1+i): {[1+i]}
```

All four units have empty support, so the library likewise classifies them as mutually
indistinguishable.

### 3. Test specialization and quotient basic-open membership

The CLI states the closure direction explicitly: the target is in the closure of the
source when the source support is a subset of the target support.

```console
$ cargo run --quiet -p gip-cli -- specializes 1 1 3 3
3+3i in closure({1+i}): true
Supp(1+i) subset of Supp(3+3i): true
```

For quotient points, `quotient-open` builds `F` as the support of the final coordinate
pair and checks disjointness:

```console
$ cargo run --quiet -p gip-cli -- quotient-open 6 0 5 0
quotient point Supp(6): {[1+i], [3]}
forbidden F = Supp(5): {[1+2i], [2+i]}
point in O_F: true
```

The mixed-parity description of `sigma_(1+i)` is also directly testable:

```console
$ cargo run --quiet -p gip-cli -- basic-open 1 1 4 3
4+3i in sigma_1+i: true
$ cargo run --quiet -p gip-cli -- basic-open 1 1 4 2
4+2i in sigma_1+i: false
```

### 4. Explore a bounded lattice window

`window` visits every nonzero `a+bi` in the square `-R <= a,b <= R`, counts membership
in one basic open, and counts distinct quotient supports. Enumeration order and output are
deterministic.

```console
$ cargo run --quiet -p gip-cli -- window 3 1 1
window: re,im in [-3,3], excluding 0
index: 1+i
points: 48
sigma members: 24
sigma nonmembers: 24
distinct quotient points: 10
```

Add `--classes` to obtain the support-frequency table for small exploratory studies:

```console
$ cargo run --quiet -p gip-cli -- window 2 1 1 --classes
window: re,im in [-2,2], excluding 0
index: 1+i
points: 24
sigma members: 12
sigma nonmembers: 12
distinct quotient points: 4
classes:
  {}: 4
  {[1+i]}: 12
  {[1+2i]}: 4
  {[2+i]}: 4
```

## Range, complexity, and overflow

- `GaussianInt` stores two `i64` coordinates. Its norm always fits exactly in `u128`.
- Checked coordinate operations use `i128` intermediates and return
  `GaussianError::ArithmeticOverflow` when an intermediate or result is not representable.
  This includes a few extreme `i64::MIN` normalization cases.
- Euclidean gcd is not restricted by the factorization ceiling, but it remains checked and
  can report overflow for extreme coordinates.
- Primality, factorization, support, indistinguishability, specialization, and quotient-point
  construction accept `N(z) <= 1,000,000,000`. Above that norm they return
  `FactorizationLimitExceeded`; they do not silently approximate.
- The deterministic factorizer scans canonical Gaussian-prime candidates through norm
  `sqrt(N(z))`. A conservative trial-division bound is `O(N(z)^(3/4))` small-integer
  remainder operations; this favors auditability over large-input performance.
- `window` accepts `0 <= R <= 100` and visits exactly `(2R+1)^2 - 1` nonzero points.
  It factors every point in both output modes to count distinct quotient supports, so runtime
  grows with both the number of points and their norms. `--classes` only prints the
  support-frequency table that the command has already accumulated.

Zero is valid for algebraic gcd and coprimality, but it is rejected by support and topology
APIs because the paper's space excludes zero.

## License

MIT. See [LICENSE](LICENSE).
