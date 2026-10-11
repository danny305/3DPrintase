import tempfile
import unittest
from pathlib import Path

from printase3d import cli as gp


class PipelineTests(unittest.TestCase):
    def test_nucleotide_atom_names_get_element_materials(self):
        self.assertEqual(gp.semantic_atom_material("P_Atom")[0], "atom_P_phosphorus")
        self.assertEqual(gp.semantic_atom_material("OP1_Atom")[0], "atom_O_oxygen")
        self.assertEqual(gp.semantic_atom_material("O5'_Atom")[0], "atom_O_oxygen")

    def test_cif_tokenizer_handles_quotes(self):
        singles, loops = gp.parse_cif("data_x\n_struct.title 'A glycan complex'\n")
        self.assertEqual(singles["_struct.title"], "A glycan complex")

    def test_render_template_is_packaged(self):
        template = Path(gp.__file__).parent / "jmol" / "3DP-Jmol.pipeline.spt"
        text = template.read_text()
        self.assertIn("__INPUT_PDB__", text)
        self.assertIn("__OUTPUT_OBJ__", text)

    def test_presets_have_search_terms(self):
        self.assertTrue(all(isinstance(p["keywords"], list) and p["keywords"] for p in gp.PRESETS.values()))

    def test_pubchem_urls_accept_cid(self):
        self.assertIn("444041", gp.PUBCHEM_SDF.format(cid=444041))

    def test_default_output_is_out_misc(self):
        self.assertEqual(gp.DEFAULT_OUTPUT, "out/misc")

    def test_pubchem_readable_stem(self):
        self.assertEqual(gp.safe_name("Alpha-Cyclodextrin"), "alpha-cyclodextrin")
        self.assertEqual(gp.safe_name("γ Cyclodextrin"), "gamma-cyclodextrin")

    def test_pubchem_short_output_pattern(self):
        self.assertEqual(f"{gp.safe_name('Alpha-Cyclodextrin')}_444913.obj", "alpha-cyclodextrin_444913.obj")

    def test_hydrogens_are_visible_by_default(self):
        args = gp.parser().parse_args(["pubchem", "5793"])
        self.assertTrue(args.hydrogens)
        args = gp.parser().parse_args(["pubchem", "5793", "--no-hydrogens"])
        self.assertFalse(args.hydrogens)

    def test_print_geometry_cli_defaults_and_overrides(self):
        args = gp.resolve_render_preset(gp.parser().parse_args(["pubchem", "5793"]))
        self.assertEqual((args.render_preset, args.spacefill, args.wireframe), ("ball-stick", 35.0, 0.3))
        args = gp.resolve_render_preset(gp.parser().parse_args(["pubchem", "5793", "--spacefill", "55", "--wireframe", "0.5"]))
        self.assertEqual((args.spacefill, args.wireframe), (55.0, 0.5))

    def test_ball_stick_preset_and_numeric_override(self):
        args = gp.resolve_render_preset(gp.parser().parse_args(["pubchem", "71485", "--ball-stick"]))
        self.assertEqual((args.render_preset, args.spacefill, args.wireframe), ("ball-stick", 35.0, 0.3))
        args = gp.resolve_render_preset(gp.parser().parse_args(["pubchem", "71485", "--ball-stick", "--wireframe", "0.35"]))
        self.assertEqual((args.spacefill, args.wireframe), (35.0, 0.35))

    def test_fructan_cli_uses_ball_stick_by_default(self):
        args = gp.resolve_render_preset(gp.parser().parse_args(["fructan", "--type", "levan", "--dp", "5", "8"]))
        self.assertEqual(args.dp, [5, 8])
        self.assertEqual((args.render_preset, args.spacefill, args.wireframe), ("ball-stick", 35.0, 0.3))

    def test_fructan_builder_formula_and_linkage(self):
        with tempfile.TemporaryDirectory() as folder:
            built = gp.build_fructan("inulin", 2, Path(folder) / "inulin_dp2.sdf", seed=61453)
        self.assertEqual(built["molecular_formula"], "C12H22O11")
        self.assertEqual(built["iupac_condensed"], "Fruf(b2-1)Fruf")

    def test_glucan_cli_and_builder(self):
        args = gp.resolve_render_preset(gp.parser().parse_args(["glucan", "--type", "cellulose", "--dp", "3", "10"]))
        self.assertEqual((args.dp, args.render_preset, args.spacefill, args.wireframe), ([3, 10], "ball-stick", 35.0, 0.3))
        with tempfile.TemporaryDirectory() as folder:
            built = gp.build_glucan("maltose", 2, Path(folder) / "maltose_dp2.sdf", seed=61453)
        self.assertEqual(built["molecular_formula"], "C12H22O11")
        self.assertEqual(built["iupac_condensed"], "Glc(a1-4)Glc")

    def test_natural_conformation_is_default_and_generic_is_separate(self):
        natural = gp.resolve_render_preset(gp.parser().parse_args([
            "glucan", "--type", "cellulose", "--dp", "4",
        ]))
        generic = gp.resolve_render_preset(gp.parser().parse_args([
            "glucan", "--type", "cellulose", "--dp", "4", "--conformation", "generic",
        ]))
        self.assertEqual(natural.conformation, "natural")
        self.assertEqual(gp.variant_name(natural, natural.seed), "ball-stick")
        self.assertEqual(gp.variant_name(generic, generic.seed), "ball-stick-generic-conformer")

    def test_natural_builder_records_restrained_glycosidic_torsions(self):
        with tempfile.TemporaryDirectory() as folder:
            built = gp.build_glucan("cellulose", 2, Path(folder) / "cellobiose.sdf")
        details = built["conformation_details"]
        self.assertEqual(details["profile"], "cellulose-beta-1-4")
        self.assertEqual(details["link_count"], 1)
        self.assertLessEqual(abs(details["measured_degrees"][0]["phi"] + 80.0), 3.1)
        self.assertLessEqual(abs(details["measured_degrees"][0]["psi"] - 100.0), 3.1)

    def test_galactan_and_mannan_use_general_pyranose_torsion_engine(self):
        self.assertIn("galactan", gp.GLUCAN_TYPES)
        self.assertIn("mannan", gp.GLUCAN_TYPES)
        with tempfile.TemporaryDirectory() as folder:
            gal = gp.build_glucan("galactan", 2, Path(folder) / "galactobiose.sdf")
            man = gp.build_glucan("mannan", 2, Path(folder) / "mannobiose.sdf")
        self.assertEqual(gal["molecular_formula"], "C12H22O11")
        self.assertEqual(man["molecular_formula"], "C12H22O11")
        self.assertEqual(gal["conformation_details"]["link_count"], 1)
        self.assertEqual(man["conformation_details"]["link_count"], 1)

    def test_fructan_omega_uses_carbon_reference_not_ring_oxygen(self):
        from glyles import convert
        from rdkit import Chem

        mol = Chem.AddHs(Chem.MolFromSmiles(convert(glycan="Fruf(b2-1)Fruf")[0][1]))
        links = gp.infer_glycosidic_torsions(mol)
        self.assertEqual(len(links), 1)
        omega = links[0]["omega"]
        self.assertTrue(all(mol.GetAtomWithIdx(i).GetAtomicNum() in (6, 8) for i in omega))
        self.assertEqual(mol.GetAtomWithIdx(omega[-1]).GetAtomicNum(), 6)

    def test_render_settings_get_distinct_variant_folders(self):
        default = gp.resolve_render_preset(gp.parser().parse_args(["pubchem", "5793"]))
        larger = gp.resolve_render_preset(gp.parser().parse_args(["pubchem", "5793", "--spacefill", "55"]))
        mono = gp.resolve_render_preset(gp.parser().parse_args(["pubchem", "5793", "--no-color"]))
        self.assertEqual(gp.variant_name(default, default.seed), "ball-stick")
        self.assertEqual(len({gp.variant_name(x, x.seed) for x in (default, larger, mono)}), 3)

    def test_bonds_are_uniform_single_cylinders_by_default(self):
        args = gp.parser().parse_args(["pubchem", "71485"])
        self.assertEqual(args.bond_colors, "uniform")
        split = gp.parser().parse_args(["pubchem", "71485", "--bond-colors", "split"])
        self.assertNotEqual(gp.variant_name(args, args.seed), gp.variant_name(split, split.seed))

    def test_general_molecule_sources_are_mutually_exclusive(self):
        args = gp.parser().parse_args(["molecule", "--smiles", "CCO", "--label", "ethanol"])
        self.assertEqual(args.smiles, "CCO")
        self.assertEqual(args.optimize, "auto")
        self.assertFalse(args.keep_all_fragments)

    def test_short_variant_names_describe_only_overrides(self):
        ball = gp.resolve_render_preset(gp.parser().parse_args(["molecule", "--cid", "5793", "--ball-stick"]))
        custom = gp.resolve_render_preset(gp.parser().parse_args(["molecule", "--cid", "5793", "--ball-stick", "--wireframe", "0.4", "--no-hydrogens"]))
        self.assertEqual(gp.variant_name(ball, ball.seed), "ball-stick")
        self.assertEqual(gp.variant_name(custom, custom.seed), "ball-stick-wf0p4-no-h")

    def test_multiple_bond_cylinders_are_off_by_default(self):
        default = gp.resolve_render_preset(gp.parser().parse_args(["molecule", "--cid", "11597698", "--ball-stick"]))
        scientific = gp.resolve_render_preset(gp.parser().parse_args(["molecule", "--cid", "11597698", "--ball-stick", "--multiple-bond-cylinders"]))
        self.assertFalse(default.multiple_bond_cylinders)
        self.assertEqual(gp.variant_name(default, default.seed), "ball-stick")
        self.assertEqual(gp.variant_name(scientific, scientific.seed), "ball-stick-multiple-bonds")

    def test_ph_variant_has_separate_short_folder(self):
        args = gp.resolve_render_preset(gp.parser().parse_args(["molecule", "--cid", "10624", "--ball-stick", "--ph", "7"]))
        self.assertEqual(gp.variant_name(args, args.seed), "ball-stick-ph7")

    def test_semantic_element_materials_distinguish_phosphorus_and_oxygen(self):
        self.assertEqual(gp.semantic_atom_material("P1_Atom")[0], "atom_P_phosphorus")
        self.assertEqual(gp.semantic_atom_material("O2_Atom")[0], "atom_O_oxygen")
        self.assertNotEqual(gp.semantic_atom_material("P1_Atom")[1], gp.semantic_atom_material("O2_Atom")[1])

    def test_mtl_preserves_element_palette_on_bambu_initial_import(self):
        # Bambu's mtl_parseline requires this exact comment without whitespace;
        # deal_algo uses each MTL entry instead of clustering nearby colors.
        materials = dict(gp.semantic_atom_material(f"{element}1_Atom")
                         for element in gp.ELEMENT_MATERIALS)
        materials["bond_uniform"] = gp.BOND_COLOR
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "palette.mtl"
            gp.write_print_mtl(path, materials)
            lines = path.read_text().splitlines()
        self.assertIn("#FirstTimeUsingMakerLab", lines)
        palette = {}
        for line in lines:
            if line.startswith("newmtl "):
                name = line.split()[1]
            elif line.startswith("Kd "):
                palette[name] = tuple(map(float, line.split()[1:]))
        self.assertEqual(palette, materials)
        self.assertEqual(palette["atom_O_oxygen"], (1.0, 0.055, 0.055))
        self.assertEqual(palette["atom_P_phosphorus"], (1.0, 0.504, 0.0))

    def test_bambu_palette_matches_bonds_to_carbon_not_hydrogen(self):
        carbon = gp.ELEMENT_MATERIALS["C"][1]
        hydrogen = gp.ELEMENT_MATERIALS["H"][1]
        bond = gp.BOND_COLOR
        self.assertGreater(sum(abs(a - b) for a, b in zip(carbon, hydrogen)), 1.5)
        self.assertEqual(carbon, bond)
        self.assertNotEqual(bond, hydrogen)

    def test_requested_elements_have_distinct_semantic_parts_and_colors(self):
        expected = {
            "N1_Atom": "N", "S1_Atom": "S", "F1_Atom": "F",
            "Cl1_Atom": "Cl", "CL1_Atom": "Cl", "Br1_Atom": "Br",
            "BR1_Atom": "Br", "I1_Atom": "I", "P1_Atom": "P",
        }
        colors = []
        for atom_name, element in expected.items():
            material, color = gp.semantic_atom_material(atom_name)
            self.assertTrue(material.startswith(f"atom_{element}_"))
            colors.append((element, color))
        unique_by_element = {element: color for element, color in colors}
        self.assertEqual(len(set(unique_by_element.values())), len(unique_by_element))

    def test_digit_prefixed_pdb_hydrogen_is_not_unknown(self):
        self.assertEqual(gp.semantic_atom_material("1HG1_Atom")[0], "atom_H_hydrogen")

    def test_false_tertiary_amide_protonation_is_corrected(self):
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CC[NH+](CC)C(=O)C[NH3+]")
        self.assertEqual(gp.correct_false_amide_protonation(mol), 1)
        self.assertEqual(Chem.GetFormalCharge(mol), 1)

    def test_protonation_transfers_all_deposited_heavy_atom_coordinates(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign

        neutral = Chem.AddHs(Chem.MolFromSmiles("N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"))
        self.assertEqual(AllChem.EmbedMolecule(neutral, randomSeed=7), 0)
        protonated = Chem.MolFromSmiles("[NH3+][C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O")
        transferred, detail = gp.transfer_heavy_atom_coordinates(neutral, protonated)
        self.assertTrue(transferred, detail)
        self.assertEqual(rdMolAlign.GetBestRMS(Chem.RemoveHs(neutral), protonated), 0.0)

    def test_auto_optimization_preserves_supplied_heavy_atom_geometry(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign

        source = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(source, randomSeed=9), 0)
        args = gp.resolve_render_preset(gp.parser().parse_args([
            "molecule", "--smiles", "CCO", "--label", "ethanol", "--optimize", "auto",
        ]))
        prepared, metadata = gp.prepare_molecule(Chem.Mol(source), args, "ethanol")
        self.assertIn("fixed-heavy-atom", metadata["coordinate_method"])
        self.assertAlmostEqual(
            rdMolAlign.GetBestRMS(Chem.RemoveHs(source), Chem.RemoveHs(prepared)), 0.0, places=6
        )

    def test_pubchem_geometry_and_ensemble_are_explicit_distinct_flags(self):
        pubchem = gp.resolve_render_preset(gp.parser().parse_args([
            "molecule", "--cid", "441477", "--ph", "7", "--pubchem-geometry",
        ]))
        ensemble = gp.resolve_render_preset(gp.parser().parse_args([
            "molecule", "--cid", "441477", "--ph", "7", "--ensemble",
        ]))
        self.assertEqual(pubchem.geometry_mode, "pubchem")
        self.assertEqual(ensemble.geometry_mode, "ensemble")
        self.assertIn("pubchem-geometry", gp.variant_name(pubchem, pubchem.seed))
        self.assertIn("ensemble", gp.variant_name(ensemble, ensemble.seed))

    def test_ensemble_ranks_charge_specific_conformers(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        source = Chem.AddHs(Chem.MolFromSmiles(
            "N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"
        ))
        self.assertEqual(AllChem.EmbedMolecule(source, randomSeed=7), 0)
        args = gp.resolve_render_preset(gp.parser().parse_args([
            "molecule", "--cid", "441477", "--ph", "7", "--ensemble",
            "--conformer-count", "8",
        ]))
        prepared, metadata = gp.prepare_molecule(source, args, "beta-D-glucosamine")
        search = metadata["conformer_search"]
        self.assertIsNotNone(search)
        self.assertGreaterEqual(search["total_candidates_including_source"], 2)
        self.assertEqual(search["forcefield"], "MMFF94s")
        self.assertEqual(Chem.GetFormalCharge(prepared), 1)


if __name__ == "__main__":
    unittest.main()
