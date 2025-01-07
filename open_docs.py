import subprocess, os
env = os.environ.copy()
env["RUSTDOCFLAGS"] = "--html-in-header katex-header.html"
subprocess.run(
    ["cargo", "doc", "--no-deps", "--open", "--document-private-items"], 
    env = env
)

# The katex-header.html in the workspace folder is used for local builds
# The katex-header.html in each crate is used for builds on docs.rs
# It's jank but it works