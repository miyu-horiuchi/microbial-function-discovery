import unittest

from microbial_function_discovery.annotations import AnnotationHit
from microbial_function_discovery.importers import parse_eggnog_mapper, parse_hmmer_domtblout


EGGNOG_TEXT = """#query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs
p1	1234.ENOG	1e-60	250	COG123	2	B	glycoside hydrolase family 5 cellulase	celA	-	3.2.1.4	ko:K01179	-	-	-	-	-	-	GH5	-	PF00150
p2	1235.ENOG	1e-55	240	COG124	2	E	nitrogenase iron protein	nifH	-	1.18.6.1	ko:K02588	-	-	-	-	-	-	-	-	PF00142
"""


HMMER_DOMTBLOUT = """# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description of target
GH5.hmm PF00150.20 300 p1 - 320 1e-40 180.0 0.0 1 1 1e-42 1e-40 180.0 0.0 5 290 10 300 8 305 0.98 glycoside hydrolase family 5 cellulase
NifH.hmm PF00142.23 280 p2 - 290 1e-30 160.0 0.0 1 1 1e-32 1e-30 160.0 0.0 2 260 8 280 7 282 0.97 nitrogenase iron protein
"""


class ImporterTests(unittest.TestCase):
    def test_parse_eggnog_mapper_returns_database_hits(self):
        hits = parse_eggnog_mapper(EGGNOG_TEXT)

        self.assertIn(
            AnnotationHit(
                protein_id="p1",
                database="CAZy",
                accession="GH5",
                name="glycoside hydrolase family 5 cellulase",
                evalue=1e-60,
            ),
            hits,
        )
        self.assertIn(
            AnnotationHit(
                protein_id="p2",
                database="KEGG",
                accession="K02588",
                name="nitrogenase iron protein",
                evalue=1e-55,
            ),
            hits,
        )

    def test_parse_hmmer_domtblout_returns_domain_hits(self):
        hits = parse_hmmer_domtblout(HMMER_DOMTBLOUT, database="Pfam")

        self.assertEqual(
            hits[0],
            AnnotationHit(
                protein_id="p1",
                database="Pfam",
                accession="PF00150.20",
                name="glycoside hydrolase family 5 cellulase",
                evalue=1e-40,
            ),
        )


if __name__ == "__main__":
    unittest.main()
