import pytest
import pandas as pd
from unittest.mock import Mock, patch
from ppi_net_builder.src.fetch import _format_req_url, fetch_string_ids, fetch_stringdb, fetch_enrichment_figure, fetch_enrichment


class TestFormatReqUrl:
    """Test cases for URL formatting function."""

    def test_format_req_url_with_version(self):
        """Test URL formatting with specific version."""
        url = _format_req_url(version="12.0", method="network")
        expected = "https://version-12-0.string-db.org/api/tsv-no-header/network"
        assert url == expected

    def test_format_req_url_without_version(self):
        """Test URL formatting without version (uses latest)."""
        url = _format_req_url(version=None, method="enrichment")
        expected = "https://string-db.org/api/tsv-no-header/enrichment"
        assert url == expected

    def test_format_req_url_with_custom_output_format(self):
        """Test URL formatting with custom output format."""
        url = _format_req_url(version="11.5", method="ppi_enrichment", output_format="json")
        expected = "https://version-11-5.string-db.org/api/json/ppi_enrichment"
        assert url == expected


class TestFetchStringIds:
    """Test cases for fetch_string_ids function."""

    @patch('ppi_net_builder.src.fetch.requests.post')
    @patch('ppi_net_builder.src.fetch.tqdm')
    def test_fetch_string_ids_with_existing_file(self, mock_tqdm, mock_post):
        """Test fetch_string_ids with existing annotation file."""
        import tempfile
        import os

        # Create temporary file with existing data
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(",gene_name,STRING_ID,species,species_name,preferred_name,annotation\n")
            f.write("0,BRCA1,9606.ENSP00000471181,9606,Homo sapiens,BRCA1,existing\n")
            temp_file = f.name

        try:
            # Mock tqdm
            mock_tqdm.return_value.__enter__ = Mock()
            mock_tqdm.return_value.__exit__ = Mock(return_value=None)

            # Mock response for new gene
            mock_response = Mock()
            mock_response.text = "inputIdentifier\tqueryIndex\tstringId\tspeciesId\tspeciesName\tpreferredName\tannotation\nTP53\t0\t9606.ENSP00000269305\t9606\tHomo sapiens\tTP53\tnew_annotation"
            mock_post.return_value = mock_response

            genes = ["BRCA1", "TP53"]  # BRCA1 exists, TP53 is new
            result = fetch_string_ids(genes, species_id=9606, file_name=temp_file)

            assert len(result) == 2
            assert "BRCA1" in result["gene_name"].values
            assert "TP53" in result["gene_name"].values

        finally:
            os.unlink(temp_file)


class TestFetchStringdb:
    """Test cases for fetch_stringdb function."""

    def test_invalid_method_raises_assertion_error(self):
        """Test that invalid method raises AssertionError."""
        with pytest.raises(AssertionError):
            fetch_stringdb(["TP53"], method="invalid_method")

    @patch('ppi_net_builder.src.fetch.requests.post')
    def test_fetch_stringdb_network_method(self, mock_post):
        """Test fetch_stringdb with network method."""
        mock_response = Mock()
        mock_response.ok = True
        mock_response.text = "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tscore\n9606.ENSP00000269305\t9606.ENSP00000471181\tTP53\tBRCA1\t0.9"
        mock_post.return_value = mock_response

        result = fetch_stringdb(["TP53", "BRCA1"], method="network", species_id=9606)

        assert len(result) == 1
        assert result.iloc[0]["stringId_A"] == "9606.ENSP00000269305"
        assert result.iloc[0]["preferredName_A"] == "TP53"
        assert result.iloc[0]["score"] == "0.9"

    @patch('ppi_net_builder.src.fetch.requests.post')
    def test_fetch_stringdb_image_output(self, mock_post):
        """Test fetch_stringdb with image output."""
        mock_response = Mock()
        mock_response.ok = True
        mock_response.content = b"fake_image_data"
        mock_post.return_value = mock_response

        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img_file = f.name

        try:
            result = fetch_stringdb(["TP53"], output="image", img_file_name=img_file)

            assert result is None  # Image output returns None

            # Check that file was written
            with open(img_file, 'rb') as f:
                assert f.read() == b"fake_image_data"

        finally:
            import os
            if os.path.exists(img_file):
                os.unlink(img_file)

    @patch('ppi_net_builder.src.fetch.requests.post')
    def test_fetch_stringdb_failed_request(self, mock_post, capsys):
        """Test fetch_stringdb with failed HTTP request."""
        mock_response = Mock()
        mock_response.ok = False
        mock_response.status_code = 500
        mock_post.return_value = mock_response

        result = fetch_stringdb(["TP53"], method="network")

        captured = capsys.readouterr()
        assert "Request is not successful. Status Code: 500" in captured.out

        # Should return empty DataFrame
        assert len(result) == 0


class TestFetchEnrichment:
    """Test cases for enrichment functions."""

    @patch('ppi_net_builder.src.fetch.requests.post')
    def test_fetch_enrichment(self, mock_post):
        """Test fetch_enrichment function."""
        mock_response = Mock()
        mock_response.text = "term\tdescription\tFDR\tsignal\tstrength\tgene_count\nGO:0006915\tapoptotic process\t1.23e-15\t45.2\t78.9\t15"
        mock_post.return_value = mock_response

        result = fetch_enrichment(["TP53", "BRCA1"], species_id=9606)

        assert len(result) == 1
        assert result.iloc[0]["term"] == "GO:0006915"
        assert result.iloc[0]["description"] == "apoptotic process"
        assert result.iloc[0]["FDR"] == "1.23e-15"

    @patch('ppi_net_builder.src.fetch.requests.post')
    def test_fetch_enrichment_figure(self, mock_post):
        """Test fetch_enrichment_figure function."""
        mock_response = Mock()
        mock_response.ok = True
        mock_response.content = b"fake_plot_data"
        mock_post.return_value = mock_response

        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img_file = f.name

        try:
            result = fetch_enrichment_figure(
                ["TP53"], img_file_name=img_file, category="Process", n_terms=10
            )

            assert result is None  # Figure function returns None

            # Check that file was written
            with open(img_file, 'rb') as f:
                assert f.read() == b"fake_plot_data"

        finally:
            import os
            if os.path.exists(img_file):
                os.unlink(img_file)
