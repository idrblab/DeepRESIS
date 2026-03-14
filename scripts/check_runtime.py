from __future__ import annotations

import importlib
import shutil
import subprocess
import sys


PYTHON_MODULES = [
    "deepresis",
    "torch",
    "rdkit",
    "molmap",
    "mordred",
    "rpy2",
    "Bio",
]


def _check_python_modules() -> list[str]:
    errors: list[str] = []
    for module_name in PYTHON_MODULES:
        try:
            importlib.import_module(module_name)
        except Exception as exc:  # pragma: no cover - runtime check
            errors.append(f"Python module check failed for {module_name}: {exc}")
    return errors


def _check_binary(name: str) -> str | None:
    if shutil.which(name):
        return None
    return f"Required executable not found on PATH: {name}"


def _check_r_package(package_name: str) -> str | None:
    cmd = [
        "R",
        "-q",
        "-e",
        f'cat(requireNamespace("{package_name}", quietly = TRUE))',
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        return "Required executable not found on PATH: R"
    combined = (result.stdout + result.stderr).strip()
    if result.returncode != 0 or "TRUE" not in combined:
        return f"Required R package check failed for {package_name}: {combined}"
    return None


def main() -> int:
    errors = _check_python_modules()

    binary_error = _check_binary("RNAfold")
    if binary_error:
        errors.append(binary_error)

    r_error = _check_r_package("LncFinder")
    if r_error:
        errors.append(r_error)

    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1

    print("DeepRESIS runtime check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
