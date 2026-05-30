def test_numpy_not_pinned_below_2():
    with open("pyproject.toml", encoding="utf-8") as file:
        pyproject = file.read()

    assert '"numpy <2.0"' not in pyproject
