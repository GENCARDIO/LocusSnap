from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class DocumentationParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids = set()
        self.links = []
        self.code_text = []
        self._inside_code = False

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if attributes.get("id"):
            assert attributes["id"] not in self.ids
            self.ids.add(attributes["id"])
        for key in ("href", "src"):
            if attributes.get(key):
                self.links.append((key, attributes[key]))
        if tag == "code":
            self._inside_code = True

    def handle_endtag(self, tag):
        if tag == "code":
            self._inside_code = False

    def handle_data(self, data):
        if self._inside_code:
            self.code_text.append(data)


def _parse_pages():
    parsed = {}
    for page in sorted(DOCS.glob("*.html")):
        parser = DocumentationParser()
        parser.feed(page.read_text(encoding="utf-8"))
        parser.close()
        parsed[page.resolve()] = parser
    return parsed


def test_documentation_pages_have_consistent_structure_and_resolved_links():
    parsed = _parse_pages()
    assert len(parsed) >= 8
    for page, parser in parsed.items():
        text = page.read_text(encoding="utf-8")
        assert text.count("<h1>") == 1
        assert '<main id="main">' in text
        assert "assets/site.css" in text
        assert "assets/site.js" in text
        for attribute, value in parser.links:
            if value.startswith(("http://", "https://", "mailto:", "data:")):
                continue
            parts = urlsplit(value)
            target = (page.parent / (parts.path or page.name)).resolve()
            assert target.exists(), f"{page.name}: broken {attribute}={value!r}"
            if attribute == "href" and parts.fragment and target.suffix == ".html":
                assert parts.fragment in parsed[target].ids, (
                    f"{page.name}: missing fragment {value!r}"
                )


def test_documented_locus_snap_flags_exist_in_the_current_parser():
    cli_source = (ROOT / "locus_snap" / "cli.py").read_text(encoding="utf-8")
    valid_options = set(re.findall(r'["\'](--[a-zA-Z][a-zA-Z0-9_-]*)["\']', cli_source))
    valid_options.add("--help")
    documented_options = set()
    for documentation_page in _parse_pages().values():
        documented_options.update(
            re.findall(
                r"--[a-zA-Z][a-zA-Z0-9_-]*",
                "\n".join(documentation_page.code_text),
            )
        )
    # These belong to pip and Python examples, rather than LocusSnap.
    documented_options.difference_update({"--upgrade", "--version"})
    unknown_options = documented_options - valid_options
    assert not unknown_options, sorted(unknown_options)
