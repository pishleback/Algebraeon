use std::{env, fs, path::Path};

fn main() {
    let conway_polynomials_url =
        "https://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/CPimport.txt";

    let out_dir = env::var("OUT_DIR").unwrap();

    if cfg!(feature = "conway-polynomials-buildtime-fetch") {
        let dest_path = Path::new(&out_dir).join("conway_polynomials.txt");
        let content = reqwest::blocking::get(conway_polynomials_url)
            .expect("Failed to download file")
            .text()
            .expect("Failed to read response");
        fs::write(&dest_path, content).expect("Failed to write file");
    }

    if cfg!(feature = "conway-polynomials-runtime-fetch") {
        let dest_path = Path::new(&out_dir).join("conway_polynomials_url.txt");
        let content = String::from(conway_polynomials_url);
        fs::write(&dest_path, content).expect("Failed to write file");
    }

    println!("cargo:rerun-if-changed=build.rs");
}
