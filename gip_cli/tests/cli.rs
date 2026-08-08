use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_gip-cli"))
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn factor_reports_unit_prime_powers_and_support() {
    let output = run(&["factor", "2", "0"]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "value: 2\nnorm: 4\nunit: -i\nfactorization: -i * (1+i)^2\nsupport: {[1+i]}\n"
    );
}

#[test]
fn paper_relation_commands_are_explicit() {
    let output = run(&["indistinguishable", "2", "0", "1", "1"]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "indistinguishable(2, 1+i): true\nsupport(2): {[1+i]}\nsupport(1+i): {[1+i]}\n"
    );

    let output = run(&["specializes", "1", "1", "3", "3"]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "3+3i in closure({1+i}): true\nSupp(1+i) subset of Supp(3+3i): true\n"
    );
}

#[test]
fn quotient_open_and_window_are_deterministic() {
    let output = run(&["quotient-open", "6", "0", "5", "0"]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "quotient point Supp(6): {[1+i], [3]}\nforbidden F = Supp(5): {[1+2i], [2+i]}\npoint in O_F: true\n"
    );

    let first = run(&["window", "2", "1", "1"]);
    let second = run(&["window", "2", "1", "1"]);
    assert!(first.status.success());
    assert_eq!(first.stdout, second.stdout);
    assert_eq!(
        String::from_utf8(first.stdout).unwrap(),
        "window: re,im in [-2,2], excluding 0\nindex: 1+i\npoints: 24\nsigma members: 12\nsigma nonmembers: 12\ndistinct quotient points: 4\n"
    );
}

#[test]
fn invalid_input_is_an_error_not_a_panic() {
    let output = run(&["factor", "0", "0"]);
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("error: zero is outside"), "{stderr}");
    assert!(!stderr.contains("panicked"), "{stderr}");
}
