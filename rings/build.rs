fn main() {
    // let conway_polynomials_url =
    //     "https://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/CPimport.txt";


    // match reqwest::blocking::get(conway_polynomials_url) {
    //     Ok(response) => {
    //         let _content = response.text().expect("Failed to read response. Please report this error at `https://github.com/pishleback/Algebraeon/issues`");
    //     }
    //     Err(err) => {
    //         println!(
    //             "cargo::warning=Unable to check for new Conway polynomials: {}.",
    //             err
    //         );
    //     }
    // }

    lalrpop::process_src().unwrap();
}
