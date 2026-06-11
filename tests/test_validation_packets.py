import unittest

from microbial_function_discovery.discovery import (
    parse_lead_table,
    render_validation_packets,
    select_validation_packet_leads,
)


LEADS_TSV = """panel\ttarget_key\tlabel\tsource\ttarget_precision\trank\tgenome_id\tspecies\tgenus\tfamily\taccession\tscore\trisk_level\tbiosafety_flags\tevidence
food_fermentation_agriculture\tfood_fermentation_agriculture:metabolite_production__butyrate\tmetabolite_production__butyrate\tannotation\t1.0\t1\tG1\tAgathobaculum butyriciproducens\tAgathobaculum\tOscillospiraceae\tGCA_001\t0.99\tmoderate\tpredicted_amr_signal\teggNOG:1RJ5B;eggNOG:1RJ6A
biosafety\tbiosafety:gram_stain__positive\tgram_stain__positive\tfusion\t1.0\t1\tG2\tNakamurella multipartita\tNakamurella\tNakamurellaceae\tGCA_002\t0.98\tlow\t\tcandidate_score
biofuels_industrial\tbiofuels_industrial:carbon_utilization__trehalose\tcarbon_utilization__trehalose\tfusion\t0.8\t2\tG3\tMesobacillus foraminis\tMesobacillus\tCytobacillaceae\tGCA_003\t0.92\tunknown\t\teggNOG:1RM9G
"""


class ValidationPacketTests(unittest.TestCase):
    def test_parse_lead_table_coerces_scores_and_keeps_ids(self):
        leads = parse_lead_table(LEADS_TSV)

        self.assertEqual(leads[0]["genome_id"], "G1")
        self.assertEqual(leads[0]["target_precision"], 1.0)
        self.assertEqual(leads[0]["rank"], 1)

    def test_select_validation_packet_leads_filters_panels_risk_and_limit(self):
        leads = parse_lead_table(LEADS_TSV)

        selected = select_validation_packet_leads(
            leads,
            limit=1,
            exclude_panels=["biosafety"],
            allowed_risks=["moderate", "unknown"],
        )

        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["genome_id"], "G1")
        self.assertEqual(selected[0]["panel"], "food_fermentation_agriculture")

    def test_render_validation_packets_includes_evidence_safety_and_next_steps(self):
        leads = parse_lead_table(LEADS_TSV)
        selected = select_validation_packet_leads(leads, limit=2)

        markdown = render_validation_packets(selected, title="Test Validation Packets")

        self.assertIn("# Test Validation Packets", markdown)
        self.assertIn("Agathobaculum butyriciproducens", markdown)
        self.assertIn("eggNOG:1RJ5B", markdown)
        self.assertIn("Biosafety review", markdown)
        self.assertIn("First validation step", markdown)
        self.assertIn("not release or deployment recommendations", markdown)


if __name__ == "__main__":
    unittest.main()
