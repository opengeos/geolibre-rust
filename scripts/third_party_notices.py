#!/usr/bin/env python3
"""Generate THIRD-PARTY-NOTICES.txt for the npm package.

cargo-about renders the crates' license texts. It does not cover Apache NOTICE
files, which Apache-2.0 section 4(d) requires us to pass on, or wasi-libc, which
is not a crate, so both are appended here.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CARGO_ABOUT = ["cargo", "about", "generate", "--locked", "--fail", "--all-features"]


def crate_notices(about):
    """Collect the distinct NOTICE files shipped by the listed crates.

    Args:
        about: The parsed JSON output of cargo-about.

    Returns:
        A dict mapping each distinct NOTICE text to the sorted crates carrying it.
    """
    notices = {}
    for license_ in about["licenses"]:
        for use in license_["used_by"]:
            crate = use["crate"]
            for path in Path(crate["manifest_path"]).parent.iterdir():
                if path.is_file() and re.match(r"(?i)notice", path.name):
                    crates = notices.setdefault(path.read_text(encoding="utf-8").strip(), set())
                    crates.add(f"{crate['name']} {crate['version']}")
    return {text: sorted(crates) for text, crates in notices.items()}


def main():
    """Write the notices to the path given as the only argument."""
    output = Path(sys.argv[1])
    subprocess.run([*CARGO_ABOUT, "-o", str(output), "about.hbs"], cwd=ROOT, check=True)
    about = json.loads(
        subprocess.run(
            [*CARGO_ABOUT, "--format", "json"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
        ).stdout
    )
    rule = "=" * 80
    with output.open("a", encoding="utf-8") as out:
        for text, crates in sorted(crate_notices(about).items()):
            out.write(f"{rule}\nNOTICE from: {', '.join(crates)}\n\n{text}\n\n")
        out.write(f"{rule}\n")
        wasi_libc = ROOT / "licenses" / "wasi-libc-NOTICES.txt"
        out.write(wasi_libc.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
