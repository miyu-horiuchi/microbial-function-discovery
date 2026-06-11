import unittest

from microbial_function_discovery.fasta import FastaRecord, parse_fasta


class FastaParsingTests(unittest.TestCase):
    def test_parse_fasta_preserves_id_description_and_sequence(self):
        records = parse_fasta(
            """>protein_1 GH5 glycoside hydrolase
MKT- aa*
>protein_2 nitrogenase iron protein
GGccTT
"""
        )

        self.assertEqual(
            records,
            [
                FastaRecord(
                    identifier="protein_1",
                    description="GH5 glycoside hydrolase",
                    sequence="MKTAA",
                ),
                FastaRecord(
                    identifier="protein_2",
                    description="nitrogenase iron protein",
                    sequence="GGCCTT",
                ),
            ],
        )

    def test_parse_fasta_accepts_sequence_without_header(self):
        records = parse_fasta("MKTAA\nGGG")

        self.assertEqual(records, [FastaRecord(identifier="sequence_1", description="", sequence="MKTAAGGG")])

    def test_parse_fasta_rejects_empty_input(self):
        with self.assertRaisesRegex(ValueError, "no FASTA records"):
            parse_fasta("\n\n")


if __name__ == "__main__":
    unittest.main()
