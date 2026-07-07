"""Standalone maintenance script to keep a local pool of IZA zeolite
framework CIFs up to date.

Scrapes the IZA Database of Zeolite Structures' framework type code table
(https://america.iza-structure.org/IZA-SC/ftc_table.php), diffs the result
against a locally saved code list (iza_codes.txt next to this script), and
downloads the CIF for each newly discovered code into cif/.

Re-running only ever appends newly listed codes to iza_codes.txt -- existing
entries are never removed or re-downloaded.

It then also adds any downloaded framework not yet present in zse's
frameworks.db, computing the T-site/O-site/ring metadata the database needs
from the CIF via zse.cif_tools + zse.rings. This is expected to run rarely
(new IZA framework types are added a few times a year), so no attempt is made
to optimize it -- it's a straight-line, few-dozen-row batch job.

Run directly:
    python update_frameworks.py
"""

import argparse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import HTTPError

from ase import Atoms  # noqa: E402
from ase.db import connect  # noqa: E402

import zse.collections  # noqa: E402
from zse.cif_tools import read_cif  # noqa: E402
from zse.rings import get_unique_rings  # noqa: E402

IZA_TABLE_URL = "https://america.iza-structure.org/IZA-SC/ftc_table.php"
IZA_CIF_URL = "https://europe.iza-structure.org/IZA-SC/cif/"

HERE = Path(__file__).resolve().parent
DEFAULT_CODES_PATH = HERE / "iza_codes.txt"
DEFAULT_CIF_DIR = HERE / "cif"
DEFAULT_DB_PATH = Path(zse.collections.__file__).parent / "frameworks.db"


class _FrameworkCodeTableParser(HTMLParser):
    """Extracts framework type codes from the IZA 'All Codes' table page.

    Only '<td class="CodeTable">' cells are collected. This skips the
    discontinued disordered-structure "star codes" table and the duplicate
    entries in the Intergrowths table further down the page, both of which
    use different CSS classes.
    """

    def __init__(self) -> None:
        super().__init__()
        self.codes: list[str] = []
        self._in_code_cell = False
        self._buffer: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "td" and dict(attrs).get("class") == "CodeTable":
            self._in_code_cell = True
            self._buffer = []

    def handle_data(self, data: str) -> None:
        if self._in_code_cell:
            self._buffer.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag == "td" and self._in_code_cell:
            text = "".join(self._buffer).strip()
            if text:
                # a leading "-" marks an "interrupted" framework on the page;
                # it isn't part of the actual code used by the CIF endpoint.
                self.codes.append(text.lstrip("-").strip())
            self._in_code_cell = False


def scrape_framework_codes(url: str = IZA_TABLE_URL) -> list[str]:
    """Fetch the IZA 'All Codes' page and return every framework type code listed.

    Args:
        url: The IZA framework code table page to scrape.

    Returns:
        Sorted list of unique three-letter framework type codes.
    """
    with urllib.request.urlopen(url) as response:
        html = response.read().decode("utf-8")

    parser = _FrameworkCodeTableParser()
    parser.feed(html)
    return sorted(set(parser.codes))


def load_codes(path: Path) -> list[str]:
    """Load a previously saved list of framework codes, one per line."""
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def save_codes(path: Path, codes: list[str]) -> None:
    """Save a sorted, de-duplicated list of framework codes, one per line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sorted(set(codes))) + "\n")


def update_code_list(path: Path = DEFAULT_CODES_PATH, url: str = IZA_TABLE_URL) -> list[str]:
    """Scrape the current IZA code table and merge any newly listed codes
    into the saved code list. Existing entries are never removed.

    Args:
        path: Where the persisted code list is stored.
        url: The IZA framework code table page to scrape.

    Returns:
        The codes that were newly discovered (not already present in 'path').
    """
    known = set(load_codes(path))
    scraped = scrape_framework_codes(url)
    new_codes = sorted(set(scraped) - known)

    if new_codes:
        save_codes(path, known | set(new_codes))

    return new_codes


def download_cif(code: str, data_dir: Path = DEFAULT_CIF_DIR) -> Path | None:
    """Download a single framework's CIF from the IZA database.

    Args:
        code: Three-letter IZA framework type code (e.g. "CHA").
        data_dir: Directory the CIF is saved into.

    Returns:
        The path the CIF was saved to, or None if the framework code doesn't
        exist on the server (HTTP 404).
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / f"{code}.cif"
    url = f"{IZA_CIF_URL}{code}.cif"
    try:
        urllib.request.urlretrieve(url, output_path)
    except HTTPError as err:
        if err.code == 404:
            print(f"error code 404: Specified Zeolite Framework '{code}' not found")
            return None
        raise
    return output_path


def build_db_entry(cif_path: Path, max_ring: int = 12) -> tuple[Atoms, dict]:
    """Compute the metadata a frameworks.db row needs from a downloaded CIF.

    'tsites'/'tmult'/'osites'/'omult' come straight from zse.cif_tools.read_cif.
    'rings' (the unique ring sizes present in the framework) is computed with
    get_unique_rings(validation="vertex"), which was verified against a
    random sample of the existing database to reproduce its stored 'rings'
    field in the large majority of cases -- it can undercount for a handful
    of extreme large-pore frameworks (e.g. ETR's 18-ring) where the true
    largest ring isn't the shortest ring at any T-site's vertex angles; raise
    'max_ring' and check by hand for those.

    Args:
        cif_path: Path to the downloaded CIF.
        max_ring: Maximum ring size (in T-atoms) to search for.

    Returns:
        (atoms, data) ready to pass to 'ase.db.connect(...).write(atoms, fw=code, data=data)'.
    """
    atoms, tsites, tmult, tinds, osites, omult, _oinds = read_cif(str(cif_path))
    ring_list, _paths, _ring_atoms, _a = get_unique_rings(
        atoms, tinds, validation="vertex", max_ring=max_ring
    )
    rings = sorted({int(r) for r in ring_list}, reverse=True)
    if rings and rings[0] == max_ring:
        print(
            f"  warning: {cif_path.stem}'s largest found ring ({rings[0]}) sits at the "
            f"max_ring search ceiling -- it may be truncated, consider a larger max_ring"
        )

    data = {"rings": rings, "tsites": tsites, "tmult": tmult, "osites": osites, "omult": omult}
    return atoms, data


def add_new_frameworks_to_db(
    codes: list[str],
    cif_dir: Path = DEFAULT_CIF_DIR,
    db_path: Path = DEFAULT_DB_PATH,
    max_ring: int = 12,
) -> list[str]:
    """Add each of 'codes' to frameworks.db, computing its metadata from the
    already-downloaded CIF. Codes already present in the database, or whose
    CIF fails to parse/analyze, are skipped and reported rather than aborting
    the whole batch.

    Args:
        codes: Framework codes to add (their CIFs must already exist in 'cif_dir').
        cif_dir: Directory the CIFs were downloaded into.
        db_path: Path to the frameworks.db to update.
        max_ring: Maximum ring size (in T-atoms) to search for.

    Returns:
        The codes that were successfully added.
    """
    db = connect(db_path)
    existing = {row.fw for row in db.select()}

    added = []
    for code in codes:
        if code in existing:
            continue
        cif_path = Path(cif_dir) / f"{code}.cif"
        if not cif_path.exists():
            print(f"  skipping {code}: no CIF found at {cif_path}")
            continue
        try:
            atoms, data = build_db_entry(cif_path, max_ring=max_ring)
        except Exception as e:  # noqa: BLE001
            print(f"  skipping {code}: failed to analyze ({type(e).__name__}: {e})")
            continue
        db.write(atoms, fw=code, data=data)
        added.append(code)
        print(f"  added {code} to database: rings={data['rings']}")

    return added


def update_framework_database(
    codes_path: Path = DEFAULT_CODES_PATH,
    cif_dir: Path = DEFAULT_CIF_DIR,
    db_path: Path = DEFAULT_DB_PATH,
    url: str = IZA_TABLE_URL,
    max_ring: int = 16,
) -> list[str]:
    """Full update: refresh the saved code list, download the CIF for every
    newly discovered framework code, then add each of those to frameworks.db.

    Args:
        codes_path: Where the persisted code list is stored.
        cif_dir: Directory newly downloaded CIFs are saved into.
        db_path: Path to the frameworks.db to update.
        url: The IZA framework code table page to scrape.
        max_ring: Maximum ring size (in T-atoms) to search for when
            computing new entries' ring metadata.

    Returns:
        The codes that were newly discovered, downloaded, and added to the database.
    """
    new_codes = update_code_list(codes_path, url)
    if not new_codes:
        print("No new framework codes found.")
        return []

    print(f"Found {len(new_codes)} new framework code(s): {', '.join(new_codes)}")
    downloaded = []
    for code in new_codes:
        path = download_cif(code, cif_dir)
        if path is not None:
            print(f"  downloaded {code} -> {path}")
            downloaded.append(code)

    if not downloaded:
        return []

    print(f"Adding {len(downloaded)} new framework(s) to {db_path}")
    return add_new_frameworks_to_db(downloaded, cif_dir, db_path, max_ring)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument("--codes-path", type=Path, default=DEFAULT_CODES_PATH)
    arg_parser.add_argument("--cif-dir", type=Path, default=DEFAULT_CIF_DIR)
    arg_parser.add_argument("--db-path", type=Path, default=DEFAULT_DB_PATH)
    arg_parser.add_argument("--url", default=IZA_TABLE_URL)
    arg_parser.add_argument("--max-ring", type=int, default=12)
    args = arg_parser.parse_args()

    update_framework_database(args.codes_path, args.cif_dir, args.db_path, args.url, args.max_ring)
