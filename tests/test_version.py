"""The version is set in scribonia/__init__.py; image tags in docs and the
version in CITATION.cff must match."""
import re
from pathlib import Path

from scribonia import __version__

ROOT = Path(__file__).resolve().parents[1]


def test_version_is_semver():
    assert re.fullmatch(r"\d+\.\d+\.\d+", __version__)


def test_image_tags_match_version():
    for name in ("README.md", "Dockerfile", "docs/usage.md"):
        text = (ROOT / name).read_text()
        tags = set(re.findall(r"gaiusaugustus/scribonia:([\w.\-]+)", text))
        assert tags == {__version__}, (name, tags)
    label = re.search(r'image\.version="([^"]+)"',
                      (ROOT / "Dockerfile").read_text())
    assert label and label.group(1) == __version__


def test_citation_file_matches_version():
    cff = re.search(r"^version: *(\S+)$", (ROOT / "CITATION.cff").read_text(), re.M)
    assert cff and cff.group(1).strip('"') == __version__
