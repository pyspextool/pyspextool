# Telluric Correction Tests

Test push to a copilot branch

## Overview
The `test_telluric.py` module now includes comprehensive tests for telluric correction across multiple instrument modes using data from the PostExtractionTests folder.

## Test Coverage

### 1. test_telluric (existing test)
- **Mode**: uspex SXD
- **Data**: Uses pre-combined spectra from processed/uspex-SXD/proc/
- **Purpose**: Tests the original telluric correction workflow with combined spectra

### 2. test_telluric_postextraction (new parametrized test)
Tests telluric correction for 9 different instrument modes:
- **spex_prism** (LowRes15) - *Minimum requirement - crucial*
- **spex_sxd** (ShortXD) - *Optional*
- **spex_lxd_1_9** (LongXD1.9) - *Optional*
- **spex_lxd_2_1** (LongXD2.1) - *Preferred*
- **spex_lxd_2_3** (LongXD2.3) - *Optional*
- **uspex_prism** - *Optional*
- **uspex_lxd_short** (LXD_Short) - *Optional*
- **uspex_lxd_long** (LXD_Long) - *Optional*
- **uspex_sxd** (ShortXD) - *Optional*

**Workflow**:
1. Combines individual object spectra using the `combine` function
2. Combines individual standard spectra
3. Runs telluric correction on the combined spectra
4. Verifies the telluric correction file is created
5. Cleans up all generated files

### 3. test_telluric_precomputed (new test)
- **Mode**: spex_prism
- **Purpose**: Tests the workflow of applying a pre-computed telluric correction to individual spectra

**Workflow**:
1. Creates a telluric correction file from combined spectra
2. Applies that pre-computed correction to individual uncorrected spectra
3. Verifies corrected individual spectra are created
4. Cleans up all generated files

This tests the use case where a telluric correction is computed once and then applied to multiple individual observations.

## Data Sources

All test data comes from the `tests/test_data/PostExtractionTests/` directory, which contains extracted spectra ready for combining and telluric correction.

### Data Organization (from PostExtractionTests/README)

| Mode | Object Files | Standard Files | Standard Name |
|------|--------------|----------------|---------------|
| spex-LowRes15 (PRISM) | 17-28 | 29-40 | HD 101060 |
| spex-ShortXD (SXD) | 626-635 | 636-645 | HD 165029 |
| spex-LongXD1.9 (LXD) | 200-205 | 176-189 | HD 7927 |
| spex-LongXD2.1 (LXD) | 585-594 | 595-614 | HD 165029 |
| spex-LongXD2.3 (LXD) | 555-564 | 565-574 | BS3314 |
| uspex-prism | 1-2 | 7-8 | HD 223352 |
| uspex-LXD_Short | 11-20 | 1-10 | HD 17778 |
| uspex-LXD_Long | 23-32 | 33-42 | HD 223352 |
| uspex-ShortXD (SXD) | 1-6 | 9-16 | HD 222332 |

## Running the Tests

### Prerequisites
1. Install pyspextool with test dependencies:
   ```bash
   pip install -e ".[test]"
   ```

2. Ensure the test_data submodule is initialized:
   ```bash
   git submodule init
   git submodule update
   ```

### Run All Telluric Tests
```bash
pytest tests/test_telluric.py -v
```

### Run Specific Test
```bash
# Run only the postextraction tests
pytest tests/test_telluric.py::test_telluric_postextraction -v

# Run test for a specific mode
pytest tests/test_telluric.py::test_telluric_postextraction[spex_prism] -v
pytest tests/test_telluric.py::test_telluric_postextraction[spex_lxd_2_1] -v
pytest tests/test_telluric.py::test_telluric_postextraction[uspex_sxd] -v

# Run the pre-computed telluric test
pytest tests/test_telluric.py::test_telluric_precomputed -v
```

### Run with More Verbosity
```bash
pytest tests/test_telluric.py -vv -s
```

## Test Configuration

Test configurations are defined in `tests/conftest.py`:

- **proc_setup**: Original fixture for processed data with pre-combined spectra
- **postextraction_setup**: New fixture for PostExtractionTests data with individual spectra

## Notes

- All tests disable QA plot displays (`qa_showblock=False`) to prevent blocking during automated testing
- Tests clean up all generated files (telluric correction files, combined spectra, corrected spectra) after execution
- The PRISM mode test is marked as crucial because PRISM data is handled differently in the telluric correction pipeline
- Tests will create temporary files in the PostExtractionTests directories during execution but will clean them up afterward

## Troubleshooting

If tests fail:
1. Verify all dependencies are installed (astropy, scipy, matplotlib, numpy, pandas, astroquery, specutils, dust_extinction, jinja2, pooch)
2. Check that the test_data submodule is properly cloned and contains the PostExtractionTests folder
3. Ensure you have write permissions in the test_data directories
4. Check that pyspextool is properly installed and importable

## Future Enhancements

Potential additions to consider:
- Tests for different correction types ('A0 V', 'reflectance', 'basic')
- Tests for different output units
- Tests with user-defined shift ranges
- Tests for writing model spectra
- Additional instrument modes (spex-LongXD1.9, spex-LongXD2.3, uspex-LXD_Long)
