import unittest

from microbial_function_discovery.panels import APPLICATION_AREAS, get_panel, list_panels


class PanelRegistryTests(unittest.TestCase):
    def test_all_schema_application_areas_have_panels(self):
        self.assertEqual(
            set(APPLICATION_AREAS),
            {
                "environmental_terraforming",
                "therapeutics_antimicrobials",
                "biofuels_industrial",
                "food_fermentation_agriculture",
                "biosafety",
            },
        )
        self.assertEqual({panel.key for panel in list_panels()}, set(APPLICATION_AREAS))

    def test_get_panel_returns_targets_and_evidence_sources(self):
        panel = get_panel("biofuels_industrial")
        self.assertIn("cellulose degradation", panel.example_targets)
        self.assertIn("CAZy", panel.evidence_sources)

    def test_unknown_panel_key_raises_key_error(self):
        with self.assertRaises(KeyError):
            get_panel("unknown")


if __name__ == "__main__":
    unittest.main()
