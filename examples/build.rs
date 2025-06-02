use std::{
    env,
    fs::{self, File},
    io::{self, Write},
    path::{Path, PathBuf},
};

fn main() -> io::Result<()> {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out_dir = env::var("OUT_DIR").unwrap();
    let doc_src_path = Path::new(&out_dir).join("generated_docs.rs");

    let mut doc_file = File::create(&doc_src_path)?;

    let md_files = find_md_files("../guide/src")?; // change "docs" to your actual directory

    for md_path in md_files {
        let md_name = md_path.file_stem().unwrap().to_str().unwrap();
        let rel_path = md_path.to_str().unwrap();
        let path = String::from(
            manifest_dir
                .join(rel_path)
                .canonicalize()
                .unwrap()
                .to_str()
                .unwrap(),
        );
        let path = path.replace('\\', "\\\\"); // escape Windows paths
        writeln!(
            doc_file,
            "mod {md_name} {{\n    #![doc = include_str!(\"{path}\")]\n}}\n",
        )?;
    }

    println!("cargo:rerun-if-changed=docs");

    Ok(())
}

fn find_md_files<P: AsRef<Path>>(dir: P) -> io::Result<Vec<PathBuf>> {
    let mut md_files = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_dir() {
            md_files.extend(find_md_files(path)?);
        } else if path.extension().is_some_and(|ext| ext == "md") {
            md_files.push(path);
        }
    }
    Ok(md_files)
}
