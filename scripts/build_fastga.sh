#!/usr/bin/env bash
set -euo pipefail
destination="${1:-tools/fastga}"
commit=178afcbba2f79abc1686e3bd2b8a916969467e83
if [ ! -d "$destination/.git" ]; then
    git clone https://github.com/thegenemyers/FASTGA.git "$destination"
fi
git -C "$destination" checkout "$commit"
if ! git -C "$destination" diff --quiet; then
    printf '%s\n' 'Source tree is modified; preserve it and choose a clean destination.' >&2
    exit 1
fi
patch_file="$(cd "$(dirname "$0")" && pwd)/fastga_empty_worker.patch"
git -C "$destination" apply "$patch_file"
make -C "$destination" -j2 FastGA FAtoGDB GIXmake ALNtoPAF GDBtoFA
printf '%s\n' "Explicit build: upstream $commit plus fastga_empty_worker.patch."
printf '%s\n' "Add $(cd "$destination" && pwd) to PATH; configure fastga_build_label accordingly."
