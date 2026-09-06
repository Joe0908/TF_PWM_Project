from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = REPO / "scripts" / name
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


build = load_script("07_build_homology_table.py")


class HomologyTests(unittest.TestCase):
    def test_domain_aliases_and_c2h2_count(self) -> None:
        self.assertEqual(build.canonical_domains("zf-C2H2_4,RPA_C"), {"c2h2"})
        self.assertEqual(build.canonical_domains("Homeobox"), {"homeodomain"})
        mapping = {
            "Q": ({"c2h2"}, 4),
            "SAME": ({"c2h2"}, 4),
            "OTHER": ({"c2h2"}, 7),
            "HOME": ({"homeodomain"}, 1),
        }
        self.assertEqual(build.domain_compatibility("Q", "SAME", mapping), "yes")
        self.assertEqual(build.domain_compatibility("Q", "OTHER", mapping), "no")
        self.assertEqual(build.domain_compatibility("Q", "HOME", mapping), "no")
        self.assertEqual(build.domain_compatibility("Q", "MISSING", mapping), "unknown")

    def test_direct_mapping_filters_and_nearest_neighbour(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            queries = temp / "queries.fasta"
            comparison = temp / "matches.tsv"
            direct = temp / "direct.tsv"
            domains = temp / "domains.tsv"
            output = temp / "output.tsv"
            provenance = temp / "provenance.tsv"

            queries.write_text(
                ">sp|Q1|one\n" + "A" * 100 + "\n"
                ">sp|Q2|two\n" + "C" * 100 + "\n"
                ">sp|Q3|three\n" + "G" * 100 + "\n"
                ">sp|Q4|four\n" + "T" * 100 + "\n",
                encoding="utf-8",
            )
            rows = [
                ["Q1", "CODEBOOK|ENSG1|PWM1", 98, 80, 1, 0, 1, 80, 1, 80, 1e-20, 120, .8, 1, 100, 80],
                ["Q2", "JASPAR|R1|PWM2", 75, 80, 5, 0, 1, 80, 1, 80, 1e-20, 100, .8, 1, 100, 80],
                ["Q2", "HOCOMOCO|R1|PWM2B", 75, 80, 5, 0, 1, 80, 1, 80, 1e-20, 100, .8, 1, 100, 80],
                ["Q2", "JASPAR|R2|PWM3", 74, 80, 5, 0, 1, 80, 1, 80, 1e-20, 200, .8, 1, 100, 80],
                ["Q2", "CODEBOOK|SHORT|PWM_SHORT", 99, 20, 0, 0, 1, 20, 1, 20, 1e-20, 300, .2, .2, 100, 100],
                ["Q2", "CODEBOOK|BAD_DBD|PWM_BAD", 90, 80, 3, 0, 1, 80, 1, 80, 1e-20, 250, .8, 1, 100, 80],
                ["Q3", "CISBP|R3|PWM4", 65, 60, 8, 0, 1, 60, 1, 60, 1e-10, 80, .6, 1, 100, 60],
                ["Q4", "JASPAR|OTHER_TF|PWM5", 100, 100, 0, 0, 1, 100, 1, 100, 0, 300, 1, 1, 100, 100],
            ]
            pd.DataFrame(rows).to_csv(comparison, sep="\t", header=False, index=False)
            pd.DataFrame(
                [["Q1", "CODEBOOK", "ENSG1", "PWM1", "exact_construct"]],
                columns=build.DIRECT_COLUMNS,
            ).to_csv(direct, sep="\t", index=False)
            pd.DataFrame(
                [
                    ["Q2", "Homeobox", 1],
                    ["R1", "Homeodomain", 1],
                    ["R2", "Homeodomain", 1],
                    ["BAD_DBD", "C2H2 ZF", 4],
                ],
                columns=build.DOMAIN_COLUMNS,
            ).to_csv(domains, sep="\t", index=False)

            argv = sys.argv
            try:
                sys.argv = [
                    "07_build_homology_table.py",
                    "--queries", str(queries),
                    "--comparison", str(comparison),
                    "--direct-assignments", str(direct),
                    "--domain-map", str(domains),
                    "--output", str(output),
                    "--provenance-output", str(provenance),
                ]
                build.main()
            finally:
                sys.argv = argv

            result = pd.read_csv(output, sep="\t", dtype=str).fillna("").set_index("TF_name")
            self.assertEqual(result.at["Q1", "Direct_PWM"], "PWM1")
            self.assertEqual(set(result.at["Q2", "Homologous_PWM"].split(";")), {"PWM2", "PWM2B"})
            self.assertEqual(result.at["Q3", "Relatively_Homologous_PWM"], "PWM4")
            self.assertEqual(result.at["Q4", "Direct_PWM"], "")
            self.assertEqual(result.at["Q4", "Homologous_PWM"], "PWM5")
            self.assertNotIn("PWM_SHORT", result.at["Q2", "Homologous_PWM"])
            self.assertNotIn("PWM_BAD", result.at["Q2", "Homologous_PWM"])


if __name__ == "__main__":
    unittest.main()
