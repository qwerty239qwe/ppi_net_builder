# Testing Guide for ppi_net_builder

This directory contains comprehensive test suites for the ppi_net_builder package.

## Test Structure

- `test_data.py`: Tests for the Species enum and DataManager class
- `test_graph.py`: Tests for the NetworkBuilder class and require_attribute decorator
- `test_fetch.py`: Tests for fetch functions (mocked HTTP requests)
- `conftest.py`: Shared pytest fixtures and configuration

## Running Tests

### Prerequisites

Install test dependencies:
```bash
pip install -e .[dev]
```

### Run all tests
```bash
pytest
```

### Run specific test file
```bash
pytest tests/test_data.py
```

### Run with coverage
```bash
pytest --cov=ppi_net_builder
```

### Run with coverage HTML report
```bash
pytest --cov=ppi_net_builder --cov-report=html
open htmlcov/index.html
```

## Test Coverage

The tests cover:
- Species enum values and common_species list
- DataManager initialization and species handling
- NetworkBuilder initialization and network construction
- Community detection and subnetwork extraction
- Enrichment analysis functions
- HTTP request mocking for external API calls
- Error handling and edge cases

## Mocking Strategy

Since the package makes HTTP requests to STRING database APIs, all external requests are mocked using `unittest.mock`. This allows for:
- Fast test execution
- Reliable test results (no network dependency)
- Testing error conditions
- Controlled test data

## Fixtures

Common test fixtures are available in `conftest.py`:
- `sample_gene_list`: List of common gene names
- `sample_dataframe`: DataFrame with gene names
- `mock_interaction_data`: Mock PPI interaction data
- `mock_enrichment_data`: Mock enrichment analysis results

## Continuous Integration

The tests are configured to run with:
- Coverage reporting (minimum 80% coverage required)
- Strict test discovery
- Warning filtering
