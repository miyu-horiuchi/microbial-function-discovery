import tempfile
import unittest
from pathlib import Path

from spaces.lead_browser import filter_leads, lead_detail, load_leads, option_values, summarize_leads


LEADS_TSV = """panel\ttarget_key\tlabel\tsource\ttarget_precision\trank\tgenome_id\tspecies\tgenus\tfamily\taccession\tscore\trisk_level\tbiosafety_flags\tevidence
biofuels_industrial\tbiofuels_industrial:carbon_utilization__trehalose\tcarbon_utilization__trehalose\tfusion:hybrid_esm2_bacformer\t1.0\t1\t1268\tMesobacillus foraminis\tMesobacillus\tBacillaceae\tGCA_004340465\t1.0\tmoderate\t\teggNOG:1RM9G;eggNOG:1RMCX
biosafety\tbiosafety:biosafety_level__bsl_1\tbiosafety_level__bsl_1\thybrid_esm2_bacformer\t1.0\t3\t11158\tActinopolymorpha singaporensis\tActinopolymorpha\tActinopolymorphaceae\tGCA_900104745\t1.0\tlow\t\tcandidate_score
food_fermentation_agriculture\tfood_fermentation_agriculture:metabolite_production__acetate\tmetabolite_production__acetate\tannotation\t1.0\t1\t132618\tAgathobaculum butyriciproducens\tAgathobaculum\tOscillospiraceae\tGCA_003096535\t1.0\tmoderate\t\teggNOG:1RJ03
"""


class LeadBrowserTests(unittest.TestCase):
    def test_load_leads_parses_numeric_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "leads.tsv"
            path.write_text(LEADS_TSV)

            leads = load_leads(path)

        self.assertEqual(len(leads), 3)
        self.assertEqual(leads[0]["genome_id"], "1268")
        self.assertEqual(leads[0]["score"], 1.0)
        self.assertEqual(leads[0]["rank"], 1)

    def test_filter_leads_supports_panel_risk_source_and_search(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "leads.tsv"
            path.write_text(LEADS_TSV)
            leads = load_leads(path)

        filtered = filter_leads(
            leads,
            panel="biosafety",
            risk="low",
            source="hybrid_esm2_bacformer",
            query="singaporensis",
        )

        self.assertEqual([lead["genome_id"] for lead in filtered], ["11158"])

    def test_summary_options_and_detail_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "leads.tsv"
            path.write_text(LEADS_TSV)
            leads = load_leads(path)

        self.assertEqual(summarize_leads(leads)["n_leads"], 3)
        self.assertIn("biofuels_industrial", option_values(leads, "panel"))
        detail = lead_detail(leads[0])

        self.assertIn("Mesobacillus foraminis", detail)
        self.assertIn("eggNOG:1RM9G", detail)
        self.assertIn("GCA_004340465", detail)


if __name__ == "__main__":
    unittest.main()
