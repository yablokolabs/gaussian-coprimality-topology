use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fmt::{self, Write as _};
use std::io::{self, Write as _};
use std::process::ExitCode;

use gaussian_integers::{
    coprime, factor, prime_support, sigma_contains, specializes_to,
    topologically_indistinguishable, GaussianError, GaussianInt, KolmogorovPoint, PrimeClass,
    QuotientBasicOpen,
};

const MAX_WINDOW_RADIUS: i64 = 100;

const HELP: &str = "Gaussian coprimality topology finite experiments

Usage:
  gip-cli gcd A_RE A_IM B_RE B_IM
  gip-cli factor RE IM
  gip-cli support RE IM
  gip-cli coprime A_RE A_IM B_RE B_IM
  gip-cli basic-open INDEX_RE INDEX_IM POINT_RE POINT_IM
  gip-cli indistinguishable A_RE A_IM B_RE B_IM
  gip-cli specializes SOURCE_RE SOURCE_IM TARGET_RE TARGET_IM
  gip-cli quotient-open POINT_RE POINT_IM FORBIDDEN_RE FORBIDDEN_IM
  gip-cli window RADIUS INDEX_RE INDEX_IM [--classes]

Coordinates are signed 64-bit decimal integers. Support-derived commands use the
library's documented finite norm range. Window radius is at most 100.
";

fn main() -> ExitCode {
    let args = env::args().skip(1).collect::<Vec<_>>();
    match run(&args) {
        Ok(output) => {
            if let Err(error) = io::stdout().lock().write_all(output.as_bytes()) {
                if error.kind() != io::ErrorKind::BrokenPipe {
                    let _ = writeln!(io::stderr().lock(), "error: {error}");
                    return ExitCode::FAILURE;
                }
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            let _ = writeln!(io::stderr().lock(), "error: {error}");
            ExitCode::from(2)
        }
    }
}

fn run(args: &[String]) -> Result<String, CliError> {
    let Some(command) = args.first().map(String::as_str) else {
        return Ok(HELP.to_owned());
    };
    if matches!(command, "help" | "--help" | "-h") {
        exact_arg_count(args, 1, "gip-cli --help")?;
        return Ok(HELP.to_owned());
    }

    match command {
        "gcd" => command_gcd(args),
        "factor" => command_factor(args),
        "support" => command_support(args),
        "coprime" => command_coprime(args),
        "basic-open" => command_basic_open(args),
        "indistinguishable" => command_indistinguishable(args),
        "specializes" => command_specializes(args),
        "quotient-open" => command_quotient_open(args),
        "window" => command_window(args),
        unknown => Err(CliError(format!(
            "unknown command '{unknown}' (run 'gip-cli --help')"
        ))),
    }
}

fn command_gcd(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(args, 5, "gip-cli gcd A_RE A_IM B_RE B_IM")?;
    let left = gaussian(args, 1, "A")?;
    let right = gaussian(args, 3, "B")?;
    let gcd = left.gcd(right)?;
    Ok(format!("gcd({left}, {right}) = {gcd}\n"))
}

fn command_factor(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(args, 3, "gip-cli factor RE IM")?;
    let value = gaussian(args, 1, "value")?;
    let factorization = factor(value)?;
    Ok(format!(
        "value: {value}\nnorm: {}\nunit: {}\nfactorization: {}\nsupport: {}\n",
        value.norm(),
        factorization.unit(),
        format_factorization(&factorization),
        format_support(&factorization.support())
    ))
}

fn command_support(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(args, 3, "gip-cli support RE IM")?;
    let value = gaussian(args, 1, "value")?;
    Ok(format!(
        "support({value}): {}\n",
        format_support(&prime_support(value)?)
    ))
}

fn command_coprime(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(args, 5, "gip-cli coprime A_RE A_IM B_RE B_IM")?;
    let left = gaussian(args, 1, "A")?;
    let right = gaussian(args, 3, "B")?;
    Ok(format!(
        "coprime({left}, {right}): {}\n",
        coprime(left, right)?
    ))
}

fn command_basic_open(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(
        args,
        5,
        "gip-cli basic-open INDEX_RE INDEX_IM POINT_RE POINT_IM",
    )?;
    let index = gaussian(args, 1, "index")?;
    let point = gaussian(args, 3, "point")?;
    Ok(format!(
        "{point} in sigma_{index}: {}\n",
        sigma_contains(index, point)?
    ))
}

fn command_indistinguishable(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(args, 5, "gip-cli indistinguishable A_RE A_IM B_RE B_IM")?;
    let left = gaussian(args, 1, "A")?;
    let right = gaussian(args, 3, "B")?;
    Ok(format!(
        "indistinguishable({left}, {right}): {}\nsupport({left}): {}\nsupport({right}): {}\n",
        topologically_indistinguishable(left, right)?,
        format_support(&prime_support(left)?),
        format_support(&prime_support(right)?)
    ))
}

fn command_specializes(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(
        args,
        5,
        "gip-cli specializes SOURCE_RE SOURCE_IM TARGET_RE TARGET_IM",
    )?;
    let source = gaussian(args, 1, "source")?;
    let target = gaussian(args, 3, "target")?;
    let result = specializes_to(source, target)?;
    Ok(format!(
        "{target} in closure({{{source}}}): {result}\nSupp({source}) subset of Supp({target}): {result}\n"
    ))
}

fn command_quotient_open(args: &[String]) -> Result<String, CliError> {
    exact_arg_count(
        args,
        5,
        "gip-cli quotient-open POINT_RE POINT_IM FORBIDDEN_RE FORBIDDEN_IM",
    )?;
    let point_value = gaussian(args, 1, "point")?;
    let forbidden_value = gaussian(args, 3, "forbidden index")?;
    let point = KolmogorovPoint::from_gaussian(point_value)?;
    let forbidden = prime_support(forbidden_value)?;
    let open = QuotientBasicOpen::new(forbidden.clone());
    Ok(format!(
        "quotient point Supp({point_value}): {}\nforbidden F = Supp({forbidden_value}): {}\npoint in O_F: {}\n",
        format_support(point.support()),
        format_support(&forbidden),
        open.contains(&point)
    ))
}

fn command_window(args: &[String]) -> Result<String, CliError> {
    if args.len() != 4 && args.len() != 5 {
        return Err(usage_error(
            "gip-cli window RADIUS INDEX_RE INDEX_IM [--classes]",
        ));
    }
    let show_classes = args.len() == 5;
    if show_classes && args[4] != "--classes" {
        return Err(usage_error(
            "gip-cli window RADIUS INDEX_RE INDEX_IM [--classes]",
        ));
    }
    let radius = parse_i64(&args[1], "radius")?;
    if !(0..=MAX_WINDOW_RADIUS).contains(&radius) {
        return Err(CliError(format!(
            "radius must be between 0 and {MAX_WINDOW_RADIUS}"
        )));
    }
    let index = gaussian(args, 2, "index")?;
    if index.is_zero() {
        return Err(GaussianError::ZeroOutsideTopology.into());
    }

    let mut points = 0_u64;
    let mut members = 0_u64;
    let mut classes = BTreeMap::<BTreeSet<PrimeClass>, u64>::new();
    for re in -radius..=radius {
        for im in -radius..=radius {
            let point = GaussianInt::new(re, im);
            if point.is_zero() {
                continue;
            }
            points += 1;
            if sigma_contains(index, point)? {
                members += 1;
            }
            *classes.entry(prime_support(point)?).or_default() += 1;
        }
    }

    let mut output = format!(
        "window: re,im in [-{radius},{radius}], excluding 0\nindex: {index}\npoints: {points}\nsigma members: {members}\nsigma nonmembers: {}\ndistinct quotient points: {}\n",
        points - members,
        classes.len()
    );
    if show_classes {
        output.push_str("classes:\n");
        for (support, count) in classes {
            writeln!(output, "  {}: {count}", format_support(&support))?;
        }
    }
    Ok(output)
}

fn gaussian(args: &[String], offset: usize, name: &str) -> Result<GaussianInt, CliError> {
    Ok(GaussianInt::new(
        parse_i64(&args[offset], &format!("{name} real coordinate"))?,
        parse_i64(&args[offset + 1], &format!("{name} imaginary coordinate"))?,
    ))
}

fn parse_i64(value: &str, name: &str) -> Result<i64, CliError> {
    value.parse().map_err(|_| {
        CliError(format!(
            "invalid {name} '{value}': expected a signed 64-bit integer"
        ))
    })
}

fn exact_arg_count(args: &[String], expected: usize, usage: &str) -> Result<(), CliError> {
    if args.len() == expected {
        Ok(())
    } else {
        Err(usage_error(usage))
    }
}

fn usage_error(usage: &str) -> CliError {
    CliError(format!("invalid arguments\nusage: {usage}"))
}

fn format_factorization(factorization: &gaussian_integers::Factorization) -> String {
    let mut output = factorization.unit().to_string();
    for power in factorization.factors() {
        write!(output, " * ({})", power.prime().representative())
            .expect("writing to a String cannot fail");
        if power.exponent() > 1 {
            write!(output, "^{}", power.exponent()).expect("writing to a String cannot fail");
        }
    }
    output
}

fn format_support(support: &BTreeSet<PrimeClass>) -> String {
    let entries = support
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    format!("{{{entries}}}")
}

#[derive(Debug)]
struct CliError(String);

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<GaussianError> for CliError {
    fn from(value: GaussianError) -> Self {
        Self(value.to_string())
    }
}

impl From<fmt::Error> for CliError {
    fn from(value: fmt::Error) -> Self {
        Self(value.to_string())
    }
}
