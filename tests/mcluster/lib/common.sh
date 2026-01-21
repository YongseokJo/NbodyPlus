#!/usr/bin/env bash
set -euo pipefail

# Test helper functions for McLuster integration tests

test_repo_root() {
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cd "$script_dir/../../.." && pwd
}

test_status() {
    local name="$1"
    local status="$2"  # PASS or FAIL
    local detail="${3:-}"

    printf "%-40s %s\n" "$name" "$status"
    if [[ -n "$detail" && "$status" == "FAIL" ]]; then
        echo "  $detail"
    fi
}

test_setup_workdir() {
    local test_name="$1"
    local workdir
    workdir=$(mktemp -d -t "abyss_test_${test_name}_XXXXXX")
    echo "$workdir"
}

test_cleanup() {
    local workdir="$1"
    [[ -d "$workdir" ]] && rm -rf "$workdir"
}

test_has_mcluster() {
    local repo_root
    repo_root="$(test_repo_root)"
    [[ -x "$repo_root/src/mcluster" ]] || [[ -x "$repo_root/mcluster/mcluster_sse" ]]
}

test_mcluster_path() {
    local repo_root
    repo_root="$(test_repo_root)"
    if [[ -x "$repo_root/src/mcluster" ]]; then
        echo "$repo_root/src/mcluster"
    elif [[ -x "$repo_root/mcluster/mcluster_sse" ]]; then
        echo "$repo_root/mcluster/mcluster_sse"
    else
        echo ""
    fi
}

test_abyss_path() {
    local repo_root
    repo_root="$(test_repo_root)"
    if [[ -x "$repo_root/abyss.exe" ]]; then
        echo "$repo_root/abyss.exe"
    else
        echo ""
    fi
}

test_die() {
    echo "FATAL: $*" >&2
    exit 1
}
