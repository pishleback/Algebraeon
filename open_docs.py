import subprocess, os
env = os.environ.copy()
env["RUSTDOCFLAGS"] = "--html-in-header katex-header.html"
subprocess.run(
    ["cargo", "doc", "--no-deps", "--open", "--document-private-items"], 
    env = env
)