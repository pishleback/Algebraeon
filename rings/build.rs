use std::{env, fs, path::Path};

fn main() {
    let conway_polynomials_url =
        "https://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/CPimport.txt";

    let out_dir = env::var("OUT_DIR").unwrap();

    if cfg!(feature = "conway-polynomials-buildtime-fetch") {
        let dest_path = Path::new(&out_dir).join("conway_polynomials.txt");

        match reqwest::blocking::get(conway_polynomials_url) {
            Ok(response) => {
                let content = response.text().expect("Failed to read response. Please report this error at `https://github.com/pishleback/Algebraeon/issues`");
                fs::write(&dest_path, content).expect("Failed to write file. Please report this error at `https://github.com/pishleback/Algebraeon/issues`");
            }
            Err(err) => {
                println!(
                    "cargo::error=Failed to download Conway polynomial file: {}. If this error persists, you can disable downloading the file by dissabling the `conway-polynomials-buildtime-fetch` feature for this crate, for example, by passing `--no-default--features` when building or running with cargo.",
                    err
                );
            }
        }
    }

    if cfg!(feature = "conway-polynomials-runtime-fetch") {
        let dest_path = Path::new(&out_dir).join("conway_polynomials_url.txt");
        let content = String::from(conway_polynomials_url);
        fs::write(&dest_path, content).expect("Failed to write file");
    }

    lalrpop::process_src().unwrap();

    println!("cargo:rerun-if-changed=build.rs");
}
