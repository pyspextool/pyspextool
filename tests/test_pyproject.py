from pathlib import Path

from packaging.requirements import Requirement

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib


def test_numpy_not_pinned_below_2():
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))

    dependencies = pyproject["project"]["dependencies"]
    requirements = [Requirement(dependency) for dependency in dependencies]
    numpy_requirements = [
        requirement for requirement in requirements if requirement.name == "numpy"
    ]

    assert numpy_requirements, "Expected numpy in project dependencies."
    assert all(
        specifier.operator not in {"<", "<="}
        for requirement in numpy_requirements
        for specifier in requirement.specifier
    )
