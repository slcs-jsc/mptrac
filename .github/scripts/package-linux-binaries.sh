#!/usr/bin/env bash

set -euo pipefail

readonly root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly src_dir="${root_dir}/src"
readonly platform="linux-x86_64"

usage() {
  echo "Usage: $0 create [output-directory] | verify ARCHIVE" >&2
  exit 2
}

is_system_library() {
  local name
  name="$(basename "$1")"

  case "${name}" in
    ld-linux-*.so*|libanl.so*|libc.so*|libdl.so*|libm.so*|libnss_*.so*|\
      libpthread.so*|libresolv.so*|librt.so*|libthread_db.so*|libutil.so*)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

dependencies() {
  ldd "$1" | awk '$2 == "=>" && $3 ~ /^\// { print $3 }'
}

require_command() {
  command -v "$1" >/dev/null || {
    echo "Required command not found: $1" >&2
    exit 1
  }
}

create_bundle() {
  local output_dir="${1:-${root_dir}}"
  local commit package stage archive library target count_before count_after
  local package_name copyright_file library_path

  require_command ldd
  require_command patchelf
  require_command sha256sum

  commit="$(git -C "${root_dir}" rev-parse --short=7 HEAD)"
  package="mptrac-${platform}-${commit}"
  stage="$(mktemp -d)"
  trap 'rm -rf "${stage}"' RETURN

  mkdir -p "${stage}/${package}/bin" "${stage}/${package}/lib" \
    "${stage}/${package}/licenses" "${output_dir}"

  while IFS= read -r executable; do
    test -x "${src_dir}/${executable}" || {
      echo "Missing executable: ${src_dir}/${executable}" >&2
      exit 1
    }
    cp "${src_dir}/${executable}" "${stage}/${package}/bin/"
  done < <(make -s -C "${src_dir}" print-exc)

  cp "${root_dir}/README.md" "${root_dir}/COPYING" "${stage}/${package}/"

  # Follow dependencies recursively. Libraries supplied by glibc stay on the
  # host; project dependencies and compiler runtimes travel with the bundle.
  while :; do
    count_before="$(find "${stage}/${package}/lib" -maxdepth 1 -type f | wc -l)"

    while IFS= read -r library; do
      is_system_library "${library}" && continue
      target="${stage}/${package}/lib/$(basename "${library}")"
      if [[ -e "${target}" ]] && ! cmp -s "${library}" "${target}"; then
        echo "Conflicting libraries named $(basename "${library}")" >&2
        exit 1
      fi
      if [[ ! -e "${target}" ]]; then
        cp -L "${library}" "${target}"
        if command -v dpkg-query >/dev/null; then
          library_path="$(readlink -f "${library}")"
          package_name="$(dpkg-query -S "${library_path}" 2>/dev/null |
            sed -n '1{s/: .*//;p}' || true)"
          package_name="${package_name%%:*}"
          copyright_file="/usr/share/doc/${package_name}/copyright"
          if [[ -n "${package_name}" && -f "${copyright_file}" ]]; then
            cp -L "${copyright_file}" \
              "${stage}/${package}/licenses/${package_name}.copyright"
          fi
        fi
      fi
    done < <(
      find "${stage}/${package}/bin" "${stage}/${package}/lib" \
        -maxdepth 1 -type f -exec ldd {} \; 2>/dev/null |
        awk '$2 == "=>" && $3 ~ /^\// { print $3 }' | sort -u
    )

    count_after="$(find "${stage}/${package}/lib" -maxdepth 1 -type f | wc -l)"
    [[ "${count_before}" == "${count_after}" ]] && break
  done

  while IFS= read -r -d '' target; do
    patchelf --set-rpath '$ORIGIN/../lib' "${target}"
  done < <(find "${stage}/${package}/bin" -maxdepth 1 -type f -print0)

  while IFS= read -r -d '' target; do
    patchelf --set-rpath '$ORIGIN' "${target}"
  done < <(find "${stage}/${package}/lib" -maxdepth 1 -type f -print0)

  {
    echo "MPTRAC working binaries"
    echo "Commit: $(git -C "${root_dir}" rev-parse HEAD)"
    echo "Source: https://github.com/slcs-jsc/mptrac/commit/$(git -C "${root_dir}" rev-parse HEAD)"
    echo "Built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Platform: Linux x86_64"
    echo "Compiler: $(gcc --version | head -n 1)"
    echo "Configuration: MPI=0 GPU=0 ZFP=0 ZSTD=0 LZ4=0 SZ3=0 ECCODES=0 KPP=0 CMS=0"
  } >"${stage}/${package}/BUILD-INFO.txt"

  {
    echo "MPTRAC precompiled development binaries"
    echo
    echo "This package contains Linux x86_64 CPU/OpenMP binaries."
    echo "MPI, GPU, compression extensions, ecCodes, KPP, and CMS are not included."
    echo
    echo "Run a program directly from the bin directory, for example:"
    echo
    echo "  ./bin/trac --help"
    echo
    echo "The bundled libraries are found relative to the executables; no"
    echo "compiler or separate GSL, netCDF, or HDF5 installation is needed."
    echo "See BUILD-INFO.txt for the exact source commit and build details."
  } >"${stage}/${package}/README-BINARY.txt"

  find "${stage}/${package}/lib" -maxdepth 1 -type f -printf '%f\n' |
    sort >"${stage}/${package}/THIRD-PARTY-LIBRARIES.txt"

  archive="${output_dir}/${package}.tar.gz"
  tar -C "${stage}" -czf "${archive}" "${package}"
  (
    cd "${output_dir}"
    sha256sum "${package}.tar.gz" >"${package}.tar.gz.sha256"
  )

  echo "Created ${archive}"
  if [[ -n "${GITHUB_OUTPUT:-}" ]]; then
    echo "package=${package}" >>"${GITHUB_OUTPUT}"
    echo "archive=${archive}" >>"${GITHUB_OUTPUT}"
  fi
}

verify_bundle() {
  local archive="$1"
  local stage package executable library

  test -f "${archive}"
  (
    cd "$(dirname "${archive}")"
    sha256sum --check "$(basename "${archive}").sha256"
  )

  stage="$(mktemp -d)"
  trap 'rm -rf "${stage}"' RETURN
  tar -C "${stage}" -xzf "${archive}"
  package="$(tar -tzf "${archive}" | sed -n '1{s|/.*||;p}')"

  test -f "${stage}/${package}/BUILD-INFO.txt"
  test -f "${stage}/${package}/README-BINARY.txt"
  test -s "${stage}/${package}/THIRD-PARTY-LIBRARIES.txt"
  test -f "${stage}/${package}/README.md"
  test -f "${stage}/${package}/COPYING"

  while IFS= read -r executable; do
    test -x "${stage}/${package}/bin/${executable}"
    if ldd "${stage}/${package}/bin/${executable}" | grep -q 'not found'; then
      echo "Unresolved dependency for ${executable}:" >&2
      ldd "${stage}/${package}/bin/${executable}" >&2
      exit 1
    fi

    while IFS= read -r library; do
      library="$(readlink -f "${library}")"
      is_system_library "${library}" && continue
      case "${library}" in
        "${stage}/${package}/lib/"*) ;;
        *)
          echo "Dependency outside bundle for ${executable}: ${library}" >&2
          exit 1
          ;;
      esac
    done < <(dependencies "${stage}/${package}/bin/${executable}")
  done < <(make -s -C "${src_dir}" print-exc)

  echo "Verified ${archive}"
}

case "${1:-}" in
  create)
    create_bundle "${2:-}"
    ;;
  verify)
    [[ $# -eq 2 ]] || usage
    verify_bundle "$2"
    ;;
  *)
    usage
    ;;
esac
