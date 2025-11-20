"""
Tests for example notebooks to ensure they are syntactically valid
and can be parsed correctly.
"""
import pytest
import json
import os
from pathlib import Path

# Get the root directory of the repository
REPO_ROOT = Path(__file__).parent.parent
NOTEBOOKS_DIR = REPO_ROOT / "notebooks"

# List of notebooks to test
NOTEBOOKS = [
    "spex-prism-example.ipynb",
    "spex-sxd-example.ipynb",
    "spex-lxd-example.ipynb",
    "uspex-prism-example.ipynb",
    "uspex-sxd-example.ipynb",
    "uspex-lxd-example.ipynb",
    "merge-orders-example.ipynb",
]


@pytest.mark.parametrize("notebook_name", NOTEBOOKS)
def test_notebook_is_valid_json(notebook_name):
    """Test that each notebook is valid JSON."""
    notebook_path = NOTEBOOKS_DIR / notebook_name
    assert notebook_path.exists(), f"Notebook {notebook_name} not found"
    
    with open(notebook_path, 'r') as f:
        try:
            nb = json.load(f)
        except json.JSONDecodeError as e:
            pytest.fail(f"Notebook {notebook_name} is not valid JSON: {e}")
    
    # Check basic notebook structure
    assert "cells" in nb, f"Notebook {notebook_name} missing 'cells' key"
    assert "metadata" in nb, f"Notebook {notebook_name} missing 'metadata' key"
    assert isinstance(nb["cells"], list), f"Notebook {notebook_name} 'cells' is not a list"


@pytest.mark.parametrize("notebook_name", NOTEBOOKS)
def test_notebook_has_code_cells(notebook_name):
    """Test that each notebook has at least one code cell."""
    notebook_path = NOTEBOOKS_DIR / notebook_name
    
    with open(notebook_path, 'r') as f:
        nb = json.load(f)
    
    code_cells = [cell for cell in nb["cells"] if cell.get("cell_type") == "code"]
    assert len(code_cells) > 0, f"Notebook {notebook_name} has no code cells"


@pytest.mark.parametrize("notebook_name", NOTEBOOKS)
def test_notebook_has_markdown_cells(notebook_name):
    """Test that each notebook has at least one markdown cell for documentation."""
    notebook_path = NOTEBOOKS_DIR / notebook_name
    
    with open(notebook_path, 'r') as f:
        nb = json.load(f)
    
    markdown_cells = [cell for cell in nb["cells"] if cell.get("cell_type") == "markdown"]
    assert len(markdown_cells) > 0, f"Notebook {notebook_name} has no markdown cells"


@pytest.mark.parametrize("notebook_name", NOTEBOOKS)
def test_notebook_code_syntax(notebook_name):
    """Test that all code cells have valid Python syntax."""
    notebook_path = NOTEBOOKS_DIR / notebook_name
    
    with open(notebook_path, 'r') as f:
        nb = json.load(f)
    
    code_cells = [cell for cell in nb["cells"] if cell.get("cell_type") == "code"]
    
    for i, cell in enumerate(code_cells):
        source = cell.get("source", [])
        if isinstance(source, list):
            code = "".join(source)
        else:
            code = source
        
        # Skip empty cells
        if not code.strip():
            continue
        
        # Try to compile the code to check syntax
        try:
            compile(code, f"{notebook_name}:cell_{i}", "exec")
        except SyntaxError as e:
            pytest.fail(
                f"Syntax error in {notebook_name} code cell {i}: {e}\n"
                f"Code:\n{code}"
            )


@pytest.mark.parametrize("notebook_name", NOTEBOOKS)
def test_notebook_imports_pyspextool(notebook_name):
    """Test that each notebook imports pyspextool."""
    notebook_path = NOTEBOOKS_DIR / notebook_name
    
    with open(notebook_path, 'r') as f:
        nb = json.load(f)
    
    code_cells = [cell for cell in nb["cells"] if cell.get("cell_type") == "code"]
    
    # Check if any code cell imports pyspextool
    has_import = False
    for cell in code_cells:
        source = cell.get("source", [])
        if isinstance(source, list):
            code = "".join(source)
        else:
            code = source
        
        if "import pyspextool" in code or "from pyspextool" in code:
            has_import = True
            break
    
    assert has_import, f"Notebook {notebook_name} does not import pyspextool"


def test_all_notebooks_present():
    """Test that all expected notebooks are present in the notebooks directory."""
    for notebook_name in NOTEBOOKS:
        notebook_path = NOTEBOOKS_DIR / notebook_name
        assert notebook_path.exists(), f"Expected notebook {notebook_name} not found"


def test_notebooks_directory_exists():
    """Test that the notebooks directory exists."""
    assert NOTEBOOKS_DIR.exists(), "Notebooks directory not found"
    assert NOTEBOOKS_DIR.is_dir(), "Notebooks path is not a directory"
