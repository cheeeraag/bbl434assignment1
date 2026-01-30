import os
from plasmid_designer import UniversalPlasmidMaker


def test_construct_creates_output(tmp_path):
    markers = "markers.tab"
    genome = "pUC19.fa"
    design = "Design_pUC19.txt"
    out = tmp_path / "out.fa"

    upm = UniversalPlasmidMaker(markers)
    upm.construct(genome, design, str(out))

    assert out.exists()
    text = out.read_text()
    # EcoRI should have been removed
    assert "GAATTC" not in text
    # output contains fasta header
    assert text.startswith(">Universal_Plasmid")
    # length is reasonable
    seq = "".join(line.strip() for line in text.splitlines()[1:])
    assert len(seq) > 100
