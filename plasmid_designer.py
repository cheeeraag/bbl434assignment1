import re

class UniversalPlasmidMaker:

    def __init__(self, markers_file):
        self.restriction_sites = {}
        self.load_markers(markers_file)

    def _clean_name(self, name):
        # remove parenthetical aliases and extra whitespace
        return re.sub(r"\s+|\(.*?\)", "", name).strip()

    def load_markers(self, markers_file):
        """Load restriction enzymes from markers.tab (TSV or Markdown table).

        The function is resilient: it will extract the first DNA recognition
        sequence it can find in each line and normalize enzyme names.
        """
        seq_re = re.compile(r"[ACGTacgt]{4,12}")
        with open(markers_file) as f:
            for lineno, line in enumerate(f, start=1):
                s = line.strip()
                if not s or s.startswith("#"):
                    continue

                # TSV: name\tseq
                if "\t" in s:
                    parts = [p.strip() for p in s.split("\t") if p.strip()]
                    if len(parts) >= 2:
                        name = self._clean_name(parts[0])
                        m = seq_re.search(parts[1])
                        if m:
                            self.restriction_sites[name] = m.group(0).upper()
                        continue

                # Markdown table row: | Category | Name | Recognizes ... |
                if s.startswith("|"):
                    cols = [c.strip() for c in s.strip().strip("|").split("|")]
                    # only take rows that explicitly mention a recognition sequence
                    if len(cols) >= 3 and "recogn" in cols[2].lower():
                        name = self._clean_name(cols[1])
                        m = seq_re.search(cols[2])
                        if m:
                            self.restriction_sites[name] = m.group(0).upper()
                        continue

                # Fallback: only accept lines that explicitly say 'Recognizes' (safer)
                if "recogn" in s.lower():
                    m = seq_re.search(s)
                    if m:
                        name_candidate = re.split(r"\s|,|;|:|\||\\", s)[0]
                        name = self._clean_name(name_candidate)
                        self.restriction_sites[name] = m.group(0).upper()
                    continue

                # if we got here, we couldn't parse this line
                # don't raise — some rows (headers) are expected — but log for debug
                # print a warning to stderr to help users understand skipped lines
                # (kept minimal to avoid noisy output)
                # print(f"Skipped markers line {lineno}: {line.strip()}")

    def load_fasta(self, path):
        with open(path) as f:
            # ignore header lines that start with '>' (possibly with leading whitespace)
            seq = "".join(line.strip() for line in f if not line.lstrip().startswith(">"))
            return seq.upper()

    def load_design(self, path):
        design = []
        with open(path) as f:
            for lineno, line in enumerate(f, start=1):
                if not line.strip():
                    continue
                if "," not in line:
                    # skip or warn about malformed lines
                    # print(f"Skipping malformed design line {lineno}: {line.strip()}")
                    continue
                a, b = [p.strip() for p in line.strip().split(",", 1)]
                design.append((a, b))
        return design

    # ---------- ORI DETECTION ----------
    def find_ori_by_gc_skew(self, genome, window=500):
        """Simplified GC-skew based ORI detection.

        If the genome is shorter than the window, return the whole genome.
        Genome is normalized to uppercase for consistent counting.
        """
        genome = genome.upper()
        if len(genome) <= window:
            return genome

        skew = []
        for i in range(len(genome) - window + 1):
            segment = genome[i:i+window]
            g = segment.count("G")
            c = segment.count("C")
            skew.append(g - c)

        ori_index = skew.index(max(skew))
        return genome[ori_index:ori_index+window]

    # ---------- CORE GENES ----------
    def extract_gene(self, genome, gene_name):
        # generalized ORF finder: ATG ... (TAA|TAG|TGA), non-greedy, codon-aligned
        orf_re = re.compile(r"ATG(?:[ATGC]{3})+?(?:TAA|TAG|TGA)", re.IGNORECASE)
        # For some features (e.g., Ampicillin / lacZ alpha) we may want the longest ORF
        matches = orf_re.findall(genome.upper())
        if not matches:
            return ""
        # heuristics: Ampicillin (bla) is usually a larger ORF than lacZ alpha fragment
        if gene_name == "Ampicillin":
            return max(matches, key=len)
        if gene_name == "Blue_White_Selection":
            # return a shorter ORF that still meets a minimum length
            short_orfs = [m for m in matches if len(m) < 1500]
            return max(short_orfs, key=len) if short_orfs else max(matches, key=len)
        # default: return the longest match
        return max(matches, key=len)

    # ---------- MCS ----------
    def build_mcs(self, design):
        mcs = ""
        for name, enzyme in design:
            key = self._clean_name(enzyme)
            if key in self.restriction_sites:
                if key.upper() != "ECORI":   # deletion constraint (normalized)
                    mcs += self.restriction_sites[key]
            else:
                # unknown enzyme — skip but don't crash
                # print(f"Warning: enzyme '{enzyme}' not found in markers")
                continue
        return mcs

    # ---------- ASSEMBLY ----------
    def construct(self, genome_fa, design_txt, output_fa):
        genome = self.load_fasta(genome_fa)
        design = self.load_design(design_txt)

        ori = self.find_ori_by_gc_skew(genome)
        mcs = self.build_mcs(design)

        ampR = ""
        lacZ = ""

        for _, feature in design:
            if feature == "Ampicillin":
                ampR = self.extract_gene(genome, feature)
            if feature == "Blue_White_Selection":
                lacZ = self.extract_gene(genome, feature)

        plasmid = ori + ampR + lacZ + mcs

        # Final EcoRI check
        plasmid = plasmid.replace("GAATTC", "")

        with open(output_fa, "w") as f:
            f.write(">Universal_Plasmid\n")
            for i in range(0, len(plasmid), 70):
                f.write(plasmid[i:i+70] + "\n")

        print("Plasmid constructed successfully.")
        print("EcoRI removed:", "GAATTC" not in plasmid)
        print("Length:", len(plasmid))


# ---------- TEST ----------
if __name__ == "__main__":
    upm = UniversalPlasmidMaker("markers.tab")
    upm.construct(
        "pUC19.fa",
        "Design_pUC19.txt",
        "Output.Fa"
    )