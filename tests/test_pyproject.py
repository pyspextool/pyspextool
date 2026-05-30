import re
from pathlib import Path


def test_numpy_not_pinned_below_2():
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject = pyproject_path.read_text(encoding="utf-8")

    numpy_requirements = re.findall(r'["\']numpy[^"\']*["\']', pyproject)

    assert numpy_requirements
    assert all("<" not in requirement for requirement in numpy_requirements)
